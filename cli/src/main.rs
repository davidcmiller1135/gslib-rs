//! `gslib-rs`: command-line front end for the `gslib_rs` library.
//!
//! This is the imperative shell: it parses arguments, reads and writes files,
//! and seeds the random number generator. All numerical work is done by pure
//! library functions.

mod missing;
mod spec;

use std::error::Error;
use std::io::Write as _;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::{Args, Parser, Subcommand, ValueEnum};
use std::num::NonZeroUsize;

use gslib_rs::geometry::{Angles, Point3, Ranges};
use gslib_rs::io::{Extraction, GslibData, MissingCode, PresentSamples, read_gslib, write_gslib};
use gslib_rs::kriging::{Estimate, ordinary_kriging, simple_kriging};
use gslib_rs::model::{CovarianceModel, Structure, VariogramModel};
use gslib_rs::parallel::{experimental_variogram, try_map};
use gslib_rs::plot::{VariogramPlot, save};
use gslib_rs::random::{Seed, Stream};
use gslib_rs::sample::{Sample, mean_value};
use gslib_rs::search::{
    DataSearch, NodeSearch, OnSharedNode, OnTooFewData, OutsideGrid, SearchEllipsoid, SuperBlock,
    TableExtent,
};
use gslib_rs::simulation::{Axis, DataConditioning, Grid, Searches, sgsim, sisim};
use gslib_rs::transform::{BackTransform, Tail, normal_scores};
use gslib_rs::variogram::LagClasses;
use missing::MissingArgs;

type AppResult<T = ()> = Result<T, Box<dyn Error>>;

/// Variograms, kriging and simulation (a Rust port of GSLIB).
#[derive(Debug, Parser)]
#[command(name = "gslib-rs", version, about)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Experimental variogram, optionally plotted against a model (GSLIB `gamv`, `vargplt`).
    Variogram(VariogramArgs),
    /// Simple or ordinary kriging at target locations (GSLIB `kt3d`).
    Krige(KrigeArgs),
    /// Sequential Gaussian simulation with normal-score transform (GSLIB `sgsim`).
    Sgsim(SgsimArgs),
    /// Sequential indicator simulation at one threshold (GSLIB `sisim`).
    Sisim(SisimArgs),
}

/// Which columns hold coordinates and values (1-based, as in GSLIB).
#[derive(Debug, Clone, Copy, Args)]
struct Columns {
    /// Column of the x coordinate.
    #[arg(long = "x-col", default_value_t = 1, value_parser = clap::value_parser!(u16).range(1..))]
    x: u16,
    /// Column of the y coordinate.
    #[arg(long = "y-col", default_value_t = 2, value_parser = clap::value_parser!(u16).range(1..))]
    y: u16,
    /// Column of the z coordinate.
    #[arg(long = "z-col", default_value_t = 3, value_parser = clap::value_parser!(u16).range(1..))]
    z: u16,
    /// Column of the variable.
    #[arg(long = "value-col", default_value_t = 4, value_parser = clap::value_parser!(u16).range(1..))]
    value: u16,
}

impl Columns {
    fn zero_based(c: u16) -> usize {
        usize::from(c.saturating_sub(1))
    }

    fn present(self, data: &GslibData, policy: &Extraction) -> gslib_rs::Result<PresentSamples> {
        data.samples(
            Self::zero_based(self.x),
            Self::zero_based(self.y),
            Self::zero_based(self.z),
            Self::zero_based(self.value),
            policy,
        )
    }

    fn points(self, data: &GslibData) -> AppResult<Vec<Point3>> {
        Ok(data.points(
            Self::zero_based(self.x),
            Self::zero_based(self.y),
            Self::zero_based(self.z),
        )?)
    }
}

/// A variogram model: nugget plus nested structures.
#[derive(Debug, Clone, Args)]
struct ModelArgs {
    /// Nugget effect c0.
    #[arg(long, default_value_t = 0.0)]
    nugget: f64,
    /// Nested structure, repeatable:
    /// `TYPE:CC:A_HMAX[:A_HMIN:A_VERT[:AZIMUTH:DIP:PLUNGE]]` with `TYPE` one
    /// of `sph`, `exp`, `gau`, `hol`; or `pow:CC:OMEGA`.
    #[arg(long = "structure", value_parser = spec::parse_structure)]
    structures: Vec<Structure>,
}

impl ModelArgs {
    fn is_given(&self) -> bool {
        self.nugget > 0.0 || !self.structures.is_empty()
    }

    fn variogram(&self) -> AppResult<VariogramModel> {
        Ok(VariogramModel::new(self.nugget, self.structures.clone())?)
    }

    fn covariance(&self) -> AppResult<CovarianceModel> {
        Ok(CovarianceModel::try_from(self.variogram()?)?)
    }
}

/// The search neighbourhood (GSLIB `ndmin ndmax`, `noct`, `radius radius1
/// radius2`, `sang1 sang2 sang3`).
#[derive(Debug, Clone, Copy, Args)]
struct SearchArgs {
    /// Fewest conditioning data per estimate (GSLIB `ndmin`).
    #[arg(long, default_value_t = 0)]
    min_samples: usize,
    /// Most conditioning data per estimate (GSLIB `ndmax`).
    #[arg(long, default_value_t = 24)]
    max_samples: usize,
    /// Most data (and simulated nodes) per octant (GSLIB `noct`; default: no octant search).
    #[arg(long)]
    per_octant: Option<NonZeroUsize>,
    /// Search radius along the major axis (default: unlimited).
    #[arg(long, default_value_t = f64::INFINITY)]
    radius: f64,
    /// Radius along the minor horizontal axis (GSLIB `radius1`; default: --radius).
    #[arg(long)]
    radius_minor: Option<f64>,
    /// Radius along the vertical axis (GSLIB `radius2`; default: --radius).
    #[arg(long)]
    radius_vertical: Option<f64>,
    /// Orientation `AZIMUTH:DIP:PLUNGE` in degrees (GSLIB `sang1..3`).
    #[arg(long, value_parser = spec::parse_angles, allow_hyphen_values = true)]
    search_angles: Option<Angles>,
}

impl SearchArgs {
    /// The search ellipsoid: a sphere unless minor or vertical radii or
    /// angles are given, in which case every radius must be finite.
    fn ellipsoid(self) -> AppResult<SearchEllipsoid> {
        if self.radius_minor.is_none()
            && self.radius_vertical.is_none()
            && self.search_angles.is_none()
        {
            return Ok(SearchEllipsoid::isotropic(self.radius)?);
        }
        if !self.radius.is_finite() {
            return Err("an anisotropic search needs a finite --radius".into());
        }
        Ok(SearchEllipsoid::new(
            self.search_angles.unwrap_or_default(),
            Ranges {
                major: self.radius,
                minor: self.radius_minor.unwrap_or(self.radius),
                vertical: self.radius_vertical.unwrap_or(self.radius),
            },
        )?)
    }

    fn search(self) -> AppResult<DataSearch> {
        let search = DataSearch::new(self.ellipsoid()?, self.min_samples, self.max_samples)?;
        Ok(match self.per_octant {
            Some(n) => search.with_octants(n)?,
            None => search,
        })
    }
}

/// What `krige` writes where a target has too few data (GSLIB `ndmin`).
#[derive(Debug, Clone, Copy, ValueEnum)]
enum TooFew {
    /// Fail, naming the target.
    Error,
    /// Write the --unest code, as GSLIB kt3d does.
    Unestimated,
}

/// How simulation uses the data (GSLIB `sstrat`) and the node table.
#[derive(Debug, Clone, Copy, Args)]
struct ConditioningArgs {
    /// Move each datum to its nearest node and search only nodes (GSLIB sstrat = 1).
    #[arg(long)]
    assign_data: bool,
    /// With --assign-data: two data on one node.
    #[arg(long, value_enum, default_value_t = SharedNode::Error, requires = "assign_data")]
    on_shared_node: SharedNode,
    /// With --assign-data: a datum outside the grid.
    #[arg(long, value_enum, default_value_t = Outside::Error, requires = "assign_data")]
    outside_grid: Outside,
    /// Extent of the previously-simulated-node table.
    #[arg(long, value_enum, default_value_t = NodeTable::Ellipsoid)]
    node_table: NodeTable,
}

/// `--on-shared-node`.
#[derive(Debug, Clone, Copy, ValueEnum)]
enum SharedNode {
    /// Fail, naming both data and the node.
    Error,
    /// Keep the datum closer to the node (GSLIB; the later on a tie).
    KeepClosest,
}

/// `--outside-grid`.
#[derive(Debug, Clone, Copy, ValueEnum)]
enum Outside {
    /// Fail, naming the datum.
    Error,
    /// Move it to the nearest edge node (GSLIB).
    Clamp,
}

/// `--node-table`.
#[derive(Debug, Clone, Copy, ValueEnum)]
enum NodeTable {
    /// The whole search ellipsoid.
    Ellipsoid,
    /// GSLIB sgsim's compile-time limits: +-35 nodes in x and y, none in z.
    Gslib,
}

impl ConditioningArgs {
    const fn data(self, search: &DataSearch) -> DataConditioning<'_> {
        if self.assign_data {
            DataConditioning::AssignToNodes {
                shared: match self.on_shared_node {
                    SharedNode::Error => OnSharedNode::Error,
                    SharedNode::KeepClosest => OnSharedNode::KeepClosest,
                },
                outside: match self.outside_grid {
                    Outside::Error => OutsideGrid::Error,
                    Outside::Clamp => OutsideGrid::ClampToEdge,
                },
            }
        } else {
            DataConditioning::Search(search)
        }
    }

    const fn extent(self) -> TableExtent {
        match self.node_table {
            NodeTable::Ellipsoid => TableExtent::Ellipsoid,
            // SGSIM.INC: MAXCTX = MAXCTY = 71, MAXCTZ = 1.
            NodeTable::Gslib => TableExtent::Gslib { max: [71, 71, 1] },
        }
    }
}

/// The simulation grid, as GSLIB's `nx xmn xsiz` lines.
#[derive(Debug, Clone, Copy, Args)]
struct GridArgs {
    /// X axis `NX:XMN:XSIZ`: NX nodes, the first at XMN, XSIZ apart.
    #[arg(long = "grid-x", value_parser = spec::parse_axis, allow_hyphen_values = true)]
    grid_x: Axis,
    /// Y axis `NY:YMN:YSIZ` (default: one node at 0).
    #[arg(long = "grid-y", value_parser = spec::parse_axis, allow_hyphen_values = true, default_value = "1:0:1")]
    grid_y: Axis,
    /// Z axis `NZ:ZMN:ZSIZ` (default: one node at 0).
    #[arg(long = "grid-z", value_parser = spec::parse_axis, allow_hyphen_values = true, default_value = "1:0:1")]
    grid_z: Axis,
    /// Most previously simulated nodes per estimate (GSLIB `nodmax`).
    #[arg(long, default_value_t = 12)]
    max_nodes: usize,
}

impl GridArgs {
    fn grid(self) -> AppResult<Grid> {
        Ok(Grid::new(self.grid_x, self.grid_y, self.grid_z)?)
    }

    /// The node search over `grid`, with the data search's ellipsoid and
    /// octant limit, as in GSLIB.
    fn node_search(
        self,
        grid: &Grid,
        search: &SearchArgs,
        conditioning: ConditioningArgs,
        model: &CovarianceModel,
    ) -> AppResult<NodeSearch> {
        Ok(NodeSearch::new(
            grid,
            &search.ellipsoid()?,
            model,
            self.max_nodes,
            search.per_octant,
            conditioning.extent(),
        )?)
    }
}

/// The grid's node coordinates as a GSLIB table with columns X, Y, Z.
fn grid_table(title: &str, grid: &Grid) -> AppResult<GslibData> {
    let values = (0..grid.len())
        .flat_map(|i| {
            let p = grid.node(i);
            [p.x, p.y, p.z]
        })
        .collect();
    Ok(GslibData::new(
        title,
        vec!["X".to_owned(), "Y".to_owned(), "Z".to_owned()],
        values,
    )?)
}

impl RealizationArgs {
    /// The stream for realization `i`: stream `i` of the master seed, so a
    /// realization does not depend on the others (D-037).
    fn stream(&self, i: usize) -> AppResult<Stream> {
        Ok(Seed(self.seed).realization(u64::try_from(i)?))
    }
}

/// Realizations and seeding.
#[derive(Debug, Clone, Args)]
struct RealizationArgs {
    /// Number of realizations.
    #[arg(long, default_value_t = 1)]
    nsim: usize,
    /// Master seed: realization i uses stream i of it (D-037). 69069 is GSLIB's
    /// example seed in SGSIM.PAR and SISIM.PAR (F-016).
    #[arg(long, default_value_t = 69069)]
    seed: u64,
    /// Output prefix: writes `PREFIX_0.gslib`, `PREFIX_1.gslib`, ...
    #[arg(long)]
    out_prefix: PathBuf,
}

#[derive(Debug, Args)]
struct VariogramArgs {
    /// Input GSLIB data file.
    #[arg(short, long)]
    input: PathBuf,
    #[command(flatten)]
    columns: Columns,
    #[command(flatten)]
    missing: MissingArgs,
    /// Maximum separation distance.
    #[arg(short = 'd', long)]
    max_dist: f64,
    /// Number of lag classes.
    #[arg(short = 'n', long)]
    n_lags: usize,
    /// Write lags as CSV here instead of standard output.
    #[arg(short = 'o', long)]
    out_csv: Option<PathBuf>,
    /// Plot to this PNG or SVG file.
    #[arg(long)]
    plot: Option<PathBuf>,
    /// Model to overlay on the plot (optional).
    #[command(flatten)]
    model: ModelArgs,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum KrigingType {
    /// Simple kriging with a known mean (`--mean`).
    Simple,
    /// Ordinary kriging (unknown constant mean).
    Ordinary,
}

#[derive(Debug, Args)]
struct KrigeArgs {
    /// Conditioning data file.
    #[arg(long)]
    data: PathBuf,
    #[command(flatten)]
    columns: Columns,
    #[command(flatten)]
    missing: MissingArgs,
    /// Target locations file (x, y, z in columns 1-3).
    #[arg(long)]
    targets: PathBuf,
    /// Kriging type.
    #[arg(long, value_enum, default_value_t = KrigingType::Ordinary)]
    kind: KrigingType,
    /// Known mean for simple kriging (default: the data mean).
    #[arg(long)]
    mean: Option<f64>,
    #[command(flatten)]
    model: ModelArgs,
    #[command(flatten)]
    search: SearchArgs,
    /// What to do where a target has fewer than --min-samples data.
    #[arg(long, value_enum, default_value_t = TooFew::Error)]
    on_too_few: TooFew,
    /// Code written for unestimated targets (GSLIB kt3d's UNEST).
    #[arg(long, default_value_t = -999.0, allow_negative_numbers = true)]
    unest: f64,
    /// Output file: target columns plus `estimate` and `variance`.
    #[arg(long)]
    out: PathBuf,
}

#[derive(Debug, Args)]
struct SgsimArgs {
    /// Conditioning data file (original units; normal-scored internally).
    #[arg(long)]
    data: PathBuf,
    #[command(flatten)]
    columns: Columns,
    #[command(flatten)]
    missing: MissingArgs,
    #[command(flatten)]
    grid: GridArgs,
    #[command(flatten)]
    conditioning: ConditioningArgs,
    /// Normal-score variogram model (its sill should be 1).
    #[command(flatten)]
    model: ModelArgs,
    #[command(flatten)]
    search: SearchArgs,
    #[command(flatten)]
    realizations: RealizationArgs,
    /// Lower tail bound for the back-transform (default: clamp to data minimum).
    #[arg(long)]
    zmin: Option<f64>,
    /// Upper tail bound for the back-transform (default: clamp to data maximum).
    #[arg(long)]
    zmax: Option<f64>,
}

#[derive(Debug, Args)]
struct SisimArgs {
    /// Conditioning data file.
    #[arg(long)]
    data: PathBuf,
    #[command(flatten)]
    columns: Columns,
    #[command(flatten)]
    missing: MissingArgs,
    #[command(flatten)]
    grid: GridArgs,
    #[command(flatten)]
    conditioning: ConditioningArgs,
    /// Threshold: indicator is 1 where the value is below it.
    #[arg(long)]
    threshold: f64,
    /// Indicator variogram model.
    #[command(flatten)]
    model: ModelArgs,
    #[command(flatten)]
    search: SearchArgs,
    #[command(flatten)]
    realizations: RealizationArgs,
}

/// The x, y, z columns of a target file are always 1-3.
const TARGET_COLUMNS: Columns = Columns {
    x: 1,
    y: 2,
    z: 3,
    value: 4,
};

fn main() -> ExitCode {
    match run(Cli::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            // Nothing more can be done if stderr itself is unwritable.
            writeln!(std::io::stderr(), "error: {e}").ok();
            ExitCode::FAILURE
        }
    }
}

fn run(cli: Cli) -> AppResult {
    match cli.command {
        Command::Variogram(args) => variogram(&args),
        Command::Krige(args) => krige(&args),
        Command::Sgsim(args) => simulate_gaussian(&args),
        Command::Sisim(args) => simulate_indicator(&args),
    }
}

/// The data file's samples whose values are present, after the sentinel
/// guard (M6) and the trimming and non-finite policies.
fn read_samples(path: &Path, columns: Columns, missing: &MissingArgs) -> AppResult<Vec<Sample>> {
    let at = |e: String| format!("{}: {e}", path.display());
    let data = read_gslib(path)?;
    missing
        .guard(&data, Columns::zero_based(columns.value))
        .map_err(at)?;
    let present = columns
        .present(&data, &missing.extraction()?)
        .map_err(|e| {
            let hint = if matches!(e, gslib_rs::Error::NonFinite { .. }) {
                "; use --nan missing or --inf missing|clamp to accept such values"
            } else {
                ""
            };
            at(format!("{e}{hint}"))
        })?;
    if present.samples().is_empty() {
        return Err(at(if present.total_rows() == 0 {
            "no data rows".to_owned()
        } else {
            format!(
                "all {} rows are missing after trimming",
                present.total_rows()
            )
        })
        .into());
    }
    Ok(present.into_samples())
}

fn variogram(args: &VariogramArgs) -> AppResult {
    let samples = read_samples(&args.input, args.columns, &args.missing)?;
    let lags = experimental_variogram(&samples, &LagClasses::new(args.max_dist, args.n_lags)?);

    let sink: Box<dyn std::io::Write> = match &args.out_csv {
        Some(path) => Box::new(std::fs::File::create(path)?),
        None => Box::new(std::io::stdout().lock()),
    };
    let mut csv = csv::Writer::from_writer(sink);
    csv.write_record(["lag", "mean_distance", "pairs", "gamma"])?;
    let optional = |v: Option<f64>| v.map(|x| x.to_string()).unwrap_or_default();
    for lag in &lags {
        csv.write_record([
            lag.center.to_string(),
            optional(lag.mean_distance),
            lag.pairs.to_string(),
            optional(lag.gamma),
        ])?;
    }
    csv.flush()?;

    if let Some(path) = &args.plot {
        let experimental = lags
            .iter()
            .filter_map(|l| Some((l.mean_distance?, l.gamma?)))
            .collect();
        let model = if args.model.is_given() {
            let m = args.model.variogram()?;
            let origin = Point3::default();
            let curve = (0..=200_u32)
                .map(|i| {
                    let h = args.max_dist * f64::from(i) / 200.0;
                    (h, m.gamma(origin, Point3::new(h, 0.0, 0.0)))
                })
                .collect();
            Some(curve)
        } else {
            None
        };
        save(
            &VariogramPlot {
                experimental,
                model,
            },
            path,
        )?;
    }
    Ok(())
}

fn krige(args: &KrigeArgs) -> AppResult {
    let samples = read_samples(&args.data, args.columns, &args.missing)?;
    let targets_file = read_gslib(&args.targets)?;
    let targets = TARGET_COLUMNS.points(&targets_file)?;
    let search = args.search.search()?;
    // The data are fixed, so index them once; targets are kriged in
    // parallel, with results independent of the thread count (D-030).
    let index = SuperBlock::new(&samples);

    let policy = match args.on_too_few {
        TooFew::Error => OnTooFewData::Error,
        TooFew::Unestimated => OnTooFewData::Unestimated,
    };
    let estimates: Vec<Option<Estimate>> = match args.kind {
        KrigingType::Ordinary => {
            let model = args.model.variogram()?;
            try_map(&targets, |&t| {
                policy
                    .resolve(search.select_samples(&index, t))?
                    .map(|near| ordinary_kriging(&near, t, &model))
                    .transpose()
            })?
        }
        KrigingType::Simple => {
            let model = args.model.covariance()?;
            let mean = match args.mean {
                Some(m) => m,
                None => mean_value(&samples).ok_or("no data to compute the mean from")?,
            };
            try_map(&targets, |&t| {
                policy
                    .resolve(search.select_samples(&index, t))?
                    .map(|near| simple_kriging(&near, t, &model, mean))
                    .transpose()
            })?
        }
    };
    let code = MissingCode::new(args.unest)?;
    let value: Vec<Option<f64>> = estimates.iter().map(|e| e.map(|e| e.value)).collect();
    let variance: Vec<Option<f64>> = estimates.iter().map(|e| e.map(|e| e.variance)).collect();
    let out = targets_file
        .with_optional_column("estimate", &value, code)?
        .with_optional_column("variance", &variance, code)?;
    write_gslib(&args.out, &out)?;
    Ok(())
}

fn realization_path(prefix: &Path, i: usize) -> PathBuf {
    PathBuf::from(format!("{}_{i}.gslib", prefix.display()))
}

fn simulate_gaussian(args: &SgsimArgs) -> AppResult {
    let data_file = read_gslib(&args.data)?;
    let samples = read_samples(&args.data, args.columns, &args.missing)?;
    let grid = args.grid.grid()?;
    let model = args.model.covariance()?;
    if (model.sill() - 1.0).abs() > 1e-6 {
        writeln!(
            std::io::stderr(),
            "warning: model sill is {} but sgsim works in normal-score units (sill 1)",
            model.sill()
        )
        .ok();
    }
    let search = args.search.search()?;
    let nodes = args
        .grid
        .node_search(&grid, &args.search, args.conditioning, &model)?;
    let searches = Searches {
        data: args.conditioning.data(&search),
        nodes: &nodes,
    };
    let table = grid_table(data_file.title(), &grid)?;

    // Normal-score transform; simulate with mean 0; back-transform.
    let values: Vec<f64> = samples.iter().map(Sample::value).collect();
    let ns = normal_scores(&values, None)?;
    let scored = samples
        .iter()
        .zip(&ns.scores)
        .map(|(s, &y)| Sample::new(s.location(), y))
        .collect::<Result<Vec<_>, _>>()?;
    let tail = |bound: Option<f64>| bound.map_or(Tail::Clamp, Tail::Linear);
    let back = BackTransform::new(ns.table, tail(args.zmin), tail(args.zmax))?;

    for i in 0..args.realizations.nsim {
        let mut rng = args.realizations.stream(i)?;
        let simulated = sgsim(&scored, &grid, &model, 0.0, searches, &mut rng)?;
        let original_units: Vec<f64> = simulated.iter().map(|&y| back.apply(y)).collect();
        let out = table.with_column(format!("sim_{i}"), &original_units)?;
        write_gslib(realization_path(&args.realizations.out_prefix, i), &out)?;
    }
    Ok(())
}

fn simulate_indicator(args: &SisimArgs) -> AppResult {
    let data_file = read_gslib(&args.data)?;
    let samples = read_samples(&args.data, args.columns, &args.missing)?;
    let grid = args.grid.grid()?;
    let model = args.model.covariance()?;
    let search = args.search.search()?;
    let nodes = args
        .grid
        .node_search(&grid, &args.search, args.conditioning, &model)?;
    let searches = Searches {
        data: args.conditioning.data(&search),
        nodes: &nodes,
    };
    let table = grid_table(data_file.title(), &grid)?;

    for i in 0..args.realizations.nsim {
        let mut rng = args.realizations.stream(i)?;
        let below = sisim(&samples, &grid, args.threshold, &model, searches, &mut rng)?;
        let values: Vec<f64> = below
            .into_iter()
            .map(|b| if b { 1.0 } else { 0.0 })
            .collect();
        let out = table.with_column(format!("ind_{i}"), &values)?;
        write_gslib(realization_path(&args.realizations.out_prefix, i), &out)?;
    }
    Ok(())
}
