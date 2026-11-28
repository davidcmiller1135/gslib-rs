//! The `variogram` subcommand: GSLIB `gamv` (`docs/spec/gamv.md`, D-052).
//!
//! This is an adapter. It reads the data file, extracts the variables with
//! the missing-value flags (D-018 to D-021), builds the indicators, the
//! directions and the request, calls `gslib_rs::variogram::gamv`, and writes
//! the result.
//!
//! `--gslib-compat` selects every legacy policy:
//! - the kernel's (`GamvPolicies::gslib`);
//! - `atol` or `dtol ≤ 0` become 45° (D-044);
//! - gamv's indicator conventions: `I = 1` if `v ≥ c`, categories by
//!   `int(v + 0.5)`, and a categorical sill of 0 (D-050);
//! - standardisation only where GSLIB applies it, silently skipped for a
//!   sill ≤ 0 (D-049, D-051);
//! - GSLIB's values for undefined lags: 0 for types 4 and 5, −999 for −3,
//!   and zeros for empty lags (D-049);
//! - in the gamv layout, GSLIB's `np` (the orderings summed, D-043).
//!
//! It does not bring back self-pairs (D-017), nor gamv's three trimming
//! intervals, which only the fixture adapter reproduces (D-048).

use std::io::Write;
use std::path::PathBuf;

use clap::{Args, ValueEnum};
use gslib_rs::geometry::Point3;
use gslib_rs::io::{GslibData, read_gslib};
use gslib_rs::plot::{VariogramPlot, VariogramSeries, save};
use gslib_rs::variogram::{CentredLags, LagScheme};
use gslib_rs::variogram::{
    Data, Direction, DirectionalVariogram, GamvLag, GamvPolicies, GamvRequest, VariogramKind,
    VariogramSpec, gamv, indicator_categorical, indicator_continuous,
};

use crate::missing::MissingArgs;
use crate::{AppResult, Columns, ModelArgs};

/// Output format.
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub(crate) enum Format {
    /// Long-form CSV: one row per variogram, direction and lag.
    Csv,
    /// GSLIB gamv's layout (title line, then `nlag + 2` rows per direction),
    /// with full-precision numbers.
    Gamv,
}

#[derive(Debug, Args)]
pub(crate) struct VariogramArgs {
    /// Input GSLIB data file.
    #[arg(short, long)]
    input: PathBuf,
    /// Coordinate columns; `--value-col` is the variable when no `--vario`
    /// is given.
    #[command(flatten)]
    columns: Columns,
    #[command(flatten)]
    missing: MissingArgs,
    /// Number of lags (GSLIB `nlag`); the output has `nlag + 2`.
    #[arg(long)]
    nlag: usize,
    /// Lag separation distance (GSLIB `xlag`).
    #[arg(long)]
    xlag: f64,
    /// Lag tolerance (GSLIB `xltol`); omitted or ≤ 0 means `xlag/2`.
    #[arg(long, allow_hyphen_values = true)]
    xltol: Option<f64>,
    /// Direction, repeatable: `AZM,ATOL,BANDH,DIP,DTOL,BANDV` in degrees and
    /// data units (azimuth clockwise from north, dip negative downward,
    /// `ATOL ≥ 90` omnidirectional). Default: one omnidirectional direction.
    #[arg(long = "dir", value_parser = parse_direction, allow_hyphen_values = true)]
    directions: Vec<[f64; 6]>,
    /// Variogram, repeatable: `TAIL,HEAD,TYPE[,CUTOFF]` with TAIL and HEAD
    /// 1-based file columns and TYPE a GSLIB type: 1 semivariogram, 2 cross,
    /// 3 covariance, -3 sill minus covariance, 4 correlogram, 5 general
    /// relative, 6 pairwise relative, 7 logarithmic, 8 madogram, 9 indicator
    /// (CUTOFF), 10 categorical indicator (CUTOFF is the code). Default: the
    /// semivariogram of `--value-col`.
    #[arg(long = "vario", value_parser = parse_vario, allow_hyphen_values = true)]
    variograms: Vec<VarioArg>,
    /// Standardise semivariograms (types 1, 9, 10 with TAIL = HEAD) by the
    /// variable's variance (GSLIB `isill = 1`). Without `--gslib-compat`,
    /// every variogram must be one of those.
    #[arg(long)]
    standardize: bool,
    /// Reproduce GSLIB gamv exactly: select every legacy policy.
    #[arg(long)]
    gslib_compat: bool,
    /// Output format.
    #[arg(long, value_enum, default_value_t = Format::Csv)]
    format: Format,
    /// Write here instead of standard output.
    #[arg(short = 'o', long)]
    out: Option<PathBuf>,
    /// Plot to this PNG or SVG file (one series per variogram and direction).
    #[arg(long)]
    plot: Option<PathBuf>,
    /// Model to overlay on the plot, along the first direction (optional).
    #[command(flatten)]
    model: ModelArgs,
}

/// One `--vario` value.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct VarioArg {
    tail: u16,
    head: u16,
    measure: Measure,
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum Measure {
    Kind(VariogramKind),
    Indicator { cutoff: f64 },
    Category { code: f64 },
}

impl VarioArg {
    const fn gslib_type(self) -> i32 {
        match self.measure {
            Measure::Kind(k) => k.gslib_type(),
            Measure::Indicator { .. } => 9,
            Measure::Category { .. } => 10,
        }
    }
}

fn number(field: &str, what: &str) -> Result<f64, String> {
    field
        .trim()
        .parse()
        .map_err(|_| format!("{what} {field:?} is not a number"))
}

fn parse_direction(spec: &str) -> Result<[f64; 6], String> {
    let numbers = spec
        .split(',')
        .map(|f| number(f, "direction field"))
        .collect::<Result<Vec<f64>, String>>()?;
    numbers
        .try_into()
        .map_err(|_| format!("direction {spec:?} must be AZM,ATOL,BANDH,DIP,DTOL,BANDV"))
}

fn parse_vario(spec: &str) -> Result<VarioArg, String> {
    let fields: Vec<&str> = spec.split(',').map(str::trim).collect();
    let column = |f: &str, what: &str| {
        f.parse::<u16>()
            .ok()
            .filter(|c| *c >= 1)
            .ok_or_else(|| format!("{what} column {f:?} must be a 1-based column number"))
    };
    let (tail, head, ty, cutoff) = match fields.as_slice() {
        [t, h, ty] => (column(t, "tail")?, column(h, "head")?, *ty, None),
        [t, h, ty, c] => (
            column(t, "tail")?,
            column(h, "head")?,
            *ty,
            Some(number(c, "cut-off")?),
        ),
        _ => {
            return Err(format!(
                "variogram {spec:?} must be TAIL,HEAD,TYPE[,CUTOFF]"
            ));
        }
    };
    let known = match ty {
        "1" => Some(VariogramKind::Semivariogram),
        "2" => Some(VariogramKind::CrossSemivariogram),
        "3" => Some(VariogramKind::Covariance),
        "-3" => Some(VariogramKind::CovarianceAsVariogram),
        "4" => Some(VariogramKind::Correlogram),
        "5" => Some(VariogramKind::GeneralRelative),
        "6" => Some(VariogramKind::PairwiseRelative),
        "7" => Some(VariogramKind::Logarithmic),
        "8" => Some(VariogramKind::Madogram),
        _ => None,
    };
    let measure = match (known, ty, cutoff) {
        (Some(kind), _, None) => Measure::Kind(kind),
        (Some(_), _, Some(_)) => {
            return Err(format!(
                "variogram {spec:?}: only types 9 and 10 take a CUTOFF"
            ));
        }
        (None, "9", Some(cutoff)) => Measure::Indicator { cutoff },
        (None, "10", Some(code)) => Measure::Category { code },
        (None, "9" | "10", None) => {
            return Err(format!("variogram {spec:?}: type {ty} needs a CUTOFF"));
        }
        (None, _, _) => {
            return Err(format!(
                "variogram {spec:?}: unknown type {ty:?} (expected 1-10 or -3)"
            ));
        }
    };
    if matches!(
        measure,
        Measure::Indicator { .. } | Measure::Category { .. }
    ) && tail != head
    {
        return Err(format!(
            "variogram {spec:?}: an indicator variogram uses one column (TAIL = HEAD)"
        ));
    }
    Ok(VarioArg {
        tail,
        head,
        measure,
    })
}

/// A variable handed to the kernel.
struct Variable {
    name: String,
    values: Vec<Option<f64>>,
    sill: Option<f64>,
}

/// GSLIB's sill: the population variance `Σv²/n − (Σv/n)²` of the present
/// values (`GAMV.FOR` ll. 293–303).
#[expect(
    clippy::suboptimal_flops,
    clippy::suspicious_operation_groupings,
    reason = "GSLIB's unfused formula: mean of squares minus squared mean"
)]
fn population_variance(values: &[Option<f64>]) -> Option<f64> {
    let (n, sum, sq) = values
        .iter()
        .flatten()
        .fold((0_u32, 0.0, 0.0), |(n, s, q), &v| (n + 1, s + v, q + v * v));
    (n > 0).then(|| {
        let n = f64::from(n);
        let avg = sum / n;
        (sq / n) - avg * avg
    })
}

/// GSLIB gamv's indicator conventions (D-050): `I = 1` if `v ≥ c`;
/// categories by `int(v + 0.5)`, truncating toward zero; categorical sill 0.
fn gslib_indicator(values: &[Option<f64>], measure: Measure) -> (Vec<Option<f64>>, Option<f64>) {
    #[expect(
        clippy::float_cmp,
        reason = "GSLIB compares the truncated integers exactly"
    )]
    let same_category = |v: f64, code: f64| (v + 0.5).trunc() == (code + 0.5).trunc();
    let is_one = |v: f64| match measure {
        Measure::Indicator { cutoff } => v >= cutoff,
        Measure::Category { code } => same_category(v, code),
        Measure::Kind(_) => false,
    };
    let column: Vec<Option<f64>> = values
        .iter()
        .map(|v| v.map(|v| if is_one(v) { 1.0 } else { 0.0 }))
        .collect();
    let sill = match measure {
        Measure::Category { .. } => Some(0.0),
        Measure::Indicator { .. } | Measure::Kind(_) => population_variance(&column),
    };
    (column, sill)
}

/// gslib-rs's indicator conventions (D-050).
fn library_indicator(
    values: &[Option<f64>],
    measure: Measure,
) -> AppResult<(Vec<Option<f64>>, Option<f64>)> {
    let ind = match measure {
        Measure::Indicator { cutoff } => indicator_continuous(values, cutoff)?,
        Measure::Category { code } => indicator_categorical(values, code)?,
        Measure::Kind(_) => return Err("internal: not an indicator variogram".into()),
    };
    Ok((ind.values, ind.sill))
}

fn short_name(name: &str) -> String {
    name.chars().take(12).collect()
}

/// Kernel variables, and each `--vario` with its kernel spec.
type Built = (Vec<Variable>, Vec<(VarioArg, VariogramSpec)>);

/// Build the kernel's variables and variogram specs from `--vario`.
fn variables(
    file: &GslibData,
    varios: &[VarioArg],
    missing: &MissingArgs,
    compat: bool,
) -> AppResult<Built> {
    let extraction = missing.extraction()?;
    let extract = |c: u16| -> AppResult<Vec<Option<f64>>> {
        let index = Columns::zero_based(c);
        missing.guard(file, index)?;
        file.extract(index, &extraction).map_err(|e| {
            let hint = if matches!(e, gslib_rs::Error::NonFinite { .. }) {
                "; use --nan missing or --inf missing|clamp to accept such values"
            } else {
                ""
            };
            format!("{e}{hint}").into()
        })
    };
    let mut vars: Vec<Variable> = Vec::new();
    let mut by_column: Vec<(u16, usize)> = Vec::new();
    let mut column_var = |c: u16, list: &mut Vec<Variable>| -> AppResult<usize> {
        if let Some(&(_, v)) = by_column.iter().find(|(col, _)| *col == c) {
            return Ok(v);
        }
        let index = Columns::zero_based(c);
        let values = extract(c)?;
        let name = file
            .names()
            .get(index)
            .map_or_else(String::new, |n| short_name(n));
        list.push(Variable {
            name,
            sill: population_variance(&values),
            values,
        });
        by_column.push((c, list.len() - 1));
        Ok(list.len() - 1)
    };
    let mut cuts = 0;
    let mut specs = Vec::new();
    for &vario in varios {
        let (tail, head, kind) = match vario.measure {
            Measure::Kind(kind) => (
                column_var(vario.tail, &mut vars)?,
                column_var(vario.head, &mut vars)?,
                kind,
            ),
            Measure::Indicator { .. } | Measure::Category { .. } => {
                let source = extract(vario.tail)?;
                let (values, sill) = if compat {
                    gslib_indicator(&source, vario.measure)
                } else {
                    library_indicator(&source, vario.measure)?
                };
                cuts += 1;
                vars.push(Variable {
                    name: format!("Indicator {cuts:>2}"),
                    values,
                    sill,
                });
                (vars.len() - 1, vars.len() - 1, VariogramKind::Semivariogram)
            }
        };
        specs.push((
            vario,
            VariogramSpec {
                tail,
                head,
                kind,
                standardize: false,
            },
        ));
    }
    Ok((vars, specs))
}

/// Apply `--standardize` (D-049, D-051).
fn standardise(
    specs: &mut [(VarioArg, VariogramSpec)],
    vars: &[Variable],
    compat: bool,
) -> AppResult {
    for (k, (_, spec)) in specs.iter_mut().enumerate() {
        let eligible = spec.kind == VariogramKind::Semivariogram && spec.tail == spec.head;
        if compat {
            let sill = vars.get(spec.tail).and_then(|v| v.sill);
            spec.standardize = eligible && sill.is_some_and(|s| s > 0.0);
        } else if eligible {
            spec.standardize = true;
        } else {
            return Err(format!(
                "--standardize: variogram {} is not a semivariogram with TAIL = HEAD \
                 (GSLIB would silently leave it; use --gslib-compat for that)",
                k + 1
            )
            .into());
        }
    }
    Ok(())
}

fn direction(fields: [f64; 6], compat: bool) -> AppResult<Direction> {
    let [
        azimuth,
        azimuth_tol,
        horizontal_band,
        dip,
        dip_tol,
        vertical_band,
    ] = fields;
    let legacy = |t: f64| if compat && t <= 0.0 { 45.0 } else { t };
    Ok(Direction::new(
        azimuth,
        legacy(azimuth_tol),
        horizontal_band,
        dip,
        legacy(dip_tol),
        vertical_band,
    )?)
}

/// Default direction: omnidirectional with unbounded bandwidths.
const OMNIDIRECTIONAL: [f64; 6] = [0.0, 90.0, 1.0e21, 0.0, 90.0, 1.0e21];

pub(crate) fn run(args: &VariogramArgs) -> AppResult {
    let compat = args.gslib_compat;
    let file = read_gslib(&args.input)?;
    if file.nrows() == 0 {
        return Err(format!("{}: no data rows", args.input.display()).into());
    }
    let points = args.columns.points(&file)?;
    let varios = if args.variograms.is_empty() {
        vec![VarioArg {
            tail: args.columns.value,
            head: args.columns.value,
            measure: Measure::Kind(VariogramKind::Semivariogram),
        }]
    } else {
        args.variograms.clone()
    };
    let at = |e: Box<dyn std::error::Error>| format!("{}: {e}", args.input.display());
    let (mut vars, mut specs) = variables(&file, &varios, &args.missing, compat).map_err(at)?;
    // As GSLIB's reader (D-048): a row with no present variable is dropped.
    // Its coordinates would otherwise still shape the super-block search,
    // and with it the summation order.
    let keep: Vec<bool> = (0..points.len())
        .map(|row| {
            vars.iter()
                .any(|v| v.values.get(row).copied().flatten().is_some())
        })
        .collect();
    let kept = keep.iter().filter(|k| **k).count();
    if kept == 0 {
        return Err(
            at(format!("all {} rows are missing after trimming", points.len()).into()).into(),
        );
    }
    let points: Vec<Point3> = points
        .into_iter()
        .zip(&keep)
        .filter_map(|(p, k)| k.then_some(p))
        .collect();
    for v in &mut vars {
        v.values = std::mem::take(&mut v.values)
            .into_iter()
            .zip(&keep)
            .filter_map(|(x, k)| k.then_some(x))
            .collect();
    }
    if args.standardize {
        standardise(&mut specs, &vars, compat)?;
    }
    let dirs = if args.directions.is_empty() {
        vec![OMNIDIRECTIONAL]
    } else {
        args.directions.clone()
    };
    let directions = dirs
        .iter()
        .map(|&d| direction(d, compat))
        .collect::<AppResult<Vec<_>>>()?;
    // D-044: an omitted or non-positive tolerance is GSLIB's "xlag/2". NaN is
    // kept, so that CentredLags rejects it.
    let xltol = args
        .xltol
        .filter(|t| t.is_nan() || *t > 0.0)
        .unwrap_or(0.5 * args.xlag);
    let request = GamvRequest {
        lags: CentredLags::new(args.nlag, args.xlag, xltol)?,
        directions,
        variograms: specs.iter().map(|(_, s)| *s).collect(),
        sills: vars.iter().map(|v| v.sill).collect(),
        policies: if compat {
            GamvPolicies::gslib()
        } else {
            GamvPolicies::default()
        },
    };
    let names: Vec<String> = vars.iter().map(|v| v.name.clone()).collect();
    let data = Data::new(points, vars.into_iter().map(|v| v.values).collect())?;
    let result = gamv(&data, &request)?;

    let out = Output {
        request: &request,
        varios: &specs,
        names: &names,
        compat,
    };
    let mut sink: Box<dyn Write> = match &args.out {
        Some(path) => Box::new(std::fs::File::create(path)?),
        None => Box::new(std::io::stdout().lock()),
    };
    match args.format {
        Format::Csv => out.csv(&result, &mut sink)?,
        Format::Gamv => out.gamv(&result, &mut sink)?,
    }
    sink.flush()?;

    if let Some(path) = &args.plot {
        save(&plot(args, &request, &result, specs.len())?, path)?;
    }
    Ok(())
}

/// A lag as GSLIB writes it: zeros for an empty lag, and GSLIB's (or, by
/// default, −999) values where the measure is undefined (D-049).
struct GslibRecord {
    distance: f64,
    gamma: f64,
    tail_mean: f64,
    head_mean: f64,
    tail_variance: f64,
    head_variance: f64,
}

/// The missing-value marker written in the gamv layout by default.
const UNDEFINED: f64 = -999.0;

/// γ as written: zero for an empty lag; for an undefined value, GSLIB's 0
/// (types 4, 5) or −999 (type −3) under `--gslib-compat`, else −999.
const fn gslib_gamma(gamma: Option<f64>, pairs: usize, gslib_type: i32, compat: bool) -> f64 {
    let undefined = if compat && gslib_type != -3 {
        0.0
    } else {
        UNDEFINED
    };
    match (gamma, pairs) {
        (Some(g), _) => g,
        (None, 0) => 0.0,
        (None, _) => undefined,
    }
}

fn gslib_record(lag: &GamvLag, gslib_type: i32, compat: bool) -> GslibRecord {
    GslibRecord {
        distance: lag.mean_distance.unwrap_or(0.0),
        gamma: gslib_gamma(lag.gamma, lag.pairs, gslib_type, compat),
        tail_mean: lag.tail_mean.unwrap_or(0.0),
        head_mean: lag.head_mean.unwrap_or(0.0),
        tail_variance: lag.tail_variance.unwrap_or(0.0),
        head_variance: lag.head_variance.unwrap_or(0.0),
    }
}

/// GSLIB's title for a variogram type (`GAMV.FOR`, `writeout`).
const fn title(gslib_type: i32) -> &'static str {
    match gslib_type.abs() {
        1 => "Semivariogram",
        2 => "Cross Semivariogram",
        3 => "Covariance",
        4 => "Correlogram",
        5 => "General Relative",
        6 => "Pairwise Relative",
        7 => "Variogram of Logarithms",
        8 => "Semimadogram",
        _ => "Indicator 1/2 Variogram",
    }
}

struct Output<'a> {
    request: &'a GamvRequest,
    varios: &'a [(VarioArg, VariogramSpec)],
    names: &'a [String],
    compat: bool,
}

impl Output<'_> {
    fn gslib_type(&self, v: usize) -> i32 {
        self.varios.get(v).map_or(1, |(a, _)| a.gslib_type())
    }

    fn csv(&self, result: &[DirectionalVariogram], sink: &mut dyn Write) -> AppResult {
        let mut csv = csv::Writer::from_writer(sink);
        csv.write_record([
            "variogram",
            "direction",
            "lag",
            "nominal_distance",
            "mean_distance",
            "gamma",
            "pairs",
            "orderings",
            "tail_mean",
            "head_mean",
            "tail_variance",
            "head_variance",
            "skipped",
        ])?;
        let text = |v: Option<f64>| v.map(|x| x.to_string()).unwrap_or_default();
        for dv in result {
            let ty = self.gslib_type(dv.variogram);
            for (i, lag) in dv.lags.iter().enumerate() {
                let nominal = self.request.lags.center(i).to_string();
                let cells: [Option<f64>; 6] = if self.compat {
                    let r = gslib_record(lag, ty, true);
                    let correlogram = ty == 4;
                    [
                        Some(r.distance),
                        Some(r.gamma),
                        Some(r.tail_mean),
                        Some(r.head_mean),
                        correlogram.then_some(r.tail_variance),
                        correlogram.then_some(r.head_variance),
                    ]
                } else {
                    [
                        lag.mean_distance,
                        lag.gamma,
                        lag.tail_mean,
                        lag.head_mean,
                        lag.tail_variance,
                        lag.head_variance,
                    ]
                };
                let [dist, gamma, tm, hm, tv, hv] = cells.map(text);
                csv.write_record([
                    (dv.variogram + 1).to_string(),
                    (dv.direction + 1).to_string(),
                    (i + 1).to_string(),
                    nominal,
                    dist,
                    gamma,
                    lag.pairs.to_string(),
                    lag.orderings.to_string(),
                    tm,
                    hm,
                    tv,
                    hv,
                    lag.skipped.to_string(),
                ])?;
            }
        }
        csv.flush()?;
        Ok(())
    }

    /// GSLIB's layout (`writeout`, format `(1x,i3,1x,f12.3,1x,f12.5,1x,i8,
    /// 4(1x,f14.5))`) with every number written in full precision.
    fn gamv(&self, result: &[DirectionalVariogram], sink: &mut dyn Write) -> AppResult {
        let name = |v: usize| self.names.get(v).map_or("", String::as_str);
        for dv in result {
            let ty = self.gslib_type(dv.variogram);
            let (tail, head) = self
                .varios
                .get(dv.variogram)
                .map_or((0, 0), |(_, s)| (s.tail, s.head));
            writeln!(
                sink,
                "{:<24}tail:{:<12} head:{:<12}  direction {:>2}",
                title(ty),
                name(tail),
                name(head),
                dv.direction + 1
            )?;
            for (i, lag) in dv.lags.iter().enumerate() {
                let r = gslib_record(lag, ty, self.compat);
                write!(
                    sink,
                    " {:>3} {:>24} {:>24} {:>8} {:>24} {:>24}",
                    i + 1,
                    r.distance,
                    r.gamma,
                    if self.compat {
                        lag.orderings
                    } else {
                        lag.pairs
                    },
                    r.tail_mean,
                    r.head_mean
                )?;
                if ty == 4 {
                    write!(sink, " {:>24} {:>24}", r.tail_variance, r.head_variance)?;
                }
                writeln!(sink)?;
            }
        }
        Ok(())
    }
}

/// Unit vector of a direction (azimuth clockwise from north, dip negative
/// downward: z = sin(dip), F-005).
fn unit_vector([azm, _, _, dip, _, _]: [f64; 6]) -> Point3 {
    let (sa, ca) = azm.to_radians().sin_cos();
    let (sd, cd) = dip.to_radians().sin_cos();
    Point3::new(sa * cd, ca * cd, sd)
}

fn plot(
    args: &VariogramArgs,
    request: &GamvRequest,
    result: &[DirectionalVariogram],
    nvario: usize,
) -> AppResult<VariogramPlot> {
    let experimental = result
        .iter()
        .map(|dv| {
            let label = if nvario > 1 {
                format!(
                    "variogram {}, direction {}",
                    dv.variogram + 1,
                    dv.direction + 1
                )
            } else {
                format!("direction {}", dv.direction + 1)
            };
            let points = dv
                .lags
                .iter()
                .filter_map(|l| Some((l.mean_distance?, l.gamma?)))
                .collect();
            VariogramSeries { label, points }
        })
        .collect();
    let model = if args.model.is_given() {
        let m = args.model.variogram()?;
        let u = unit_vector(args.directions.first().copied().unwrap_or(OMNIDIRECTIONAL));
        let reach = request.lags.max_distance();
        let origin = Point3::default();
        Some(
            (0..=200_u32)
                .map(|i| {
                    let h = reach * f64::from(i) / 200.0;
                    (h, m.gamma(origin, Point3::new(h * u.x, h * u.y, h * u.z)))
                })
                .collect(),
        )
    } else {
        None
    };
    Ok(VariogramPlot {
        experimental,
        model,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_directions() {
        assert_eq!(
            parse_direction("0, 22.5, 25, -10, 22.5, 25"),
            Ok([0.0, 22.5, 25.0, -10.0, 22.5, 25.0])
        );
        assert!(parse_direction("0,22.5,25").is_err(), "too few");
        assert!(parse_direction("0,a,25,0,1,1").is_err(), "not a number");
    }

    #[test]
    fn parses_variograms() {
        let v = parse_vario("4,5,2").unwrap();
        assert_eq!((v.tail, v.head, v.gslib_type()), (4, 5, 2));
        assert_eq!(parse_vario("4,4,-3").unwrap().gslib_type(), -3);
        assert_eq!(
            parse_vario("4,4,9,2.5").unwrap().measure,
            Measure::Indicator { cutoff: 2.5 }
        );
        assert!(parse_vario("4,4,9").is_err(), "indicator needs a cut-off");
        assert!(parse_vario("4,4,1,2").is_err(), "cut-off only for 9 and 10");
        assert!(parse_vario("4,5,9,1").is_err(), "indicator uses one column");
        assert!(parse_vario("4,4,11").is_err(), "unknown type");
        assert!(parse_vario("0,4,1").is_err(), "columns are 1-based");
    }

    #[test]
    fn gslib_indicator_conventions() {
        let v = [Some(1.0), Some(2.0), Some(3.0), Some(-0.7)];
        let (ind, _) = gslib_indicator(&v, Measure::Indicator { cutoff: 2.0 });
        assert_eq!(ind, [Some(0.0), Some(1.0), Some(1.0), Some(0.0)], "v ≥ c");
        let (cat, sill) = gslib_indicator(&v, Measure::Category { code: 0.0 });
        assert_eq!(
            cat,
            [Some(0.0), Some(0.0), Some(0.0), Some(1.0)],
            "−0.7 → 0"
        );
        assert_eq!(sill, Some(0.0), "categorical sill left at 0");
    }

    #[test]
    fn gslib_values_for_undefined_lags() {
        assert_eq!(
            gslib_gamma(None, 0, 4, true).to_bits(),
            0.0_f64.to_bits(),
            "empty"
        );
        assert_eq!(
            gslib_gamma(None, 3, 4, true).to_bits(),
            0.0_f64.to_bits(),
            "type 4"
        );
        assert_eq!(
            gslib_gamma(None, 3, 5, true).to_bits(),
            0.0_f64.to_bits(),
            "type 5"
        );
        assert_eq!(
            gslib_gamma(None, 3, -3, true).to_bits(),
            UNDEFINED.to_bits(),
            "type −3"
        );
        assert_eq!(
            gslib_gamma(None, 3, 4, false).to_bits(),
            UNDEFINED.to_bits(),
            "default"
        );
        assert_eq!(
            gslib_gamma(Some(0.5), 3, 4, true).to_bits(),
            0.5_f64.to_bits(),
            "defined"
        );
    }
}
