use clap::{Parser, Subcommand};
use std::{error::Error, path::{Path, PathBuf}};
use ndarray::Array1;


use gslib_rs::{
    io::{read_gslib, write_gslib, GslibData},
    variogram,
    plot::plot_variogram,
    simulation::{sgsim, sisim},
};

/// GSLIB-RS command-line interface
#[derive(Parser)]
#[command(name = "gslib-rs")]
#[command(version)]
#[command(about = "Variograms, kriging, simulation & plotting in Rust", long_about = None)]
struct Cli {
    #[command(subcommand)]
    cmd: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Compute an empirical variogram
    Variogram {
        /// Input GSLIB data file (ASCII)
        #[arg(short, long)]
        input: PathBuf,

        /// Maximum distance to consider
        #[arg(short='d', long)]
        max_dist: f64,

        /// Number of lag classes
        #[arg(short='n', long)]
        n_lags: usize,

        /// Optional CSV to write lag,semivar pairs
        #[arg(short='o', long)]
        out_csv: Option<PathBuf>,

        /// Fit and overlay a model ("spherical" or "exponential")
        #[arg(long)]
        model: Option<String>,

        /// Optional PNG/SVG to plot variogram
        #[arg(long)]
        plot: Option<PathBuf>,
    },

    /// Sequential Gaussian Simulation
    Sgsim {
        /// Conditioning data file (GSLIB)
        #[arg(long)]
        data: PathBuf,

        /// Simulation targets file (GSLIB; must list x,y,z only)
        #[arg(long)]
        sim: PathBuf,

        /// Number of realizations
        #[arg(long, default_value_t = 1)]
        nsim: usize,

        /// Seed for RNG
        #[arg(long)]
        seed: Option<u64>,

        /// Output prefix (files will be PREFIX_0.gslib, PREFIX_1.gslib, …)
        #[arg(long)]
        out_prefix: PathBuf,
    },

    /// Sequential Indicator Simulation
    Sisim {
        #[arg(long)] data: PathBuf,
        #[arg(long)] sim: PathBuf,
        #[arg(long)] threshold: f64,
        #[arg(long, default_value_t = 1)] nsim: usize,
        #[arg(long)] out_prefix: PathBuf,
    },

    /// SASIM (Simulated Annealing).  See `--help` for parameters.
    Sasim {
        #[arg(long)] data: PathBuf,
        #[arg(long)] out_prefix: PathBuf,
        // ... you can add flags here to capture SasimParams fields ...
    },
}

type Coords = Vec<(f64, f64, f64)>;

/// Read a GSLIB file and check it has at least `min_cols` columns.
fn read_checked(path: &Path, min_cols: usize, what: &str) -> Result<GslibData, Box<dyn Error>> {
    let ds = read_gslib(path)?;
    let ncols = ds.data.ncols();
    if ncols < min_cols {
        return Err(format!(
            "{what} file {} has {ncols} column(s); expected at least {min_cols}",
            path.display()
        )
        .into());
    }
    Ok(ds)
}

/// Columns 0..3 as (x, y, z).
fn coords_of(ds: &GslibData) -> Coords {
    ds.data.rows().into_iter().map(|r| (r[0], r[1], r[2])).collect()
}

/// Read a data file with columns x, y, z, value.
fn read_data(path: &Path) -> Result<(GslibData, Coords, Vec<f64>), Box<dyn Error>> {
    let ds = read_checked(path, 4, "data")?;
    let coords = coords_of(&ds);
    let values = ds.data.column(3).to_vec();
    Ok((ds, coords, values))
}

/// Write one realization: every column of the simulation-target file, plus
/// `values` appended as a new column called `name`.
fn write_realization(
    sim_ds: &GslibData,
    title: &str,
    name: String,
    values: Vec<f64>,
    path: PathBuf,
) -> Result<(), Box<dyn Error>> {
    let mut data = sim_ds.data.clone();
    data.push_column(Array1::from(values).view())?;
    let mut var_names = sim_ds.var_names.clone();
    var_names.push(name);
    let out = GslibData { title: title.to_string(), var_names, data };
    write_gslib(path, &out)?;
    Ok(())
}

fn realization_path(prefix: &Path, i: usize) -> PathBuf {
    PathBuf::from(format!("{}_{}.gslib", prefix.to_string_lossy(), i))
}

fn main() -> Result<(), Box<dyn Error>> {
    let cli = Cli::parse();

    match cli.cmd {
        Commands::Variogram { input, max_dist, n_lags, out_csv, model, plot } => {
            // 1) Read data (columns: x, y, z, value)
            let (_ds, coords, values) = read_data(&input)?;

            // 2) Compute empirical variogram
            let (lags, semis) = variogram::compute_empirical(&coords, &values, max_dist, n_lags);

            // 3) Write CSV if requested
            if let Some(path) = out_csv {
                let mut wtr = csv::Writer::from_path(path)?;
                wtr.write_record(["lag", "semivar"])?;
                for (&h,&g) in lags.iter().zip(&semis) {
                    wtr.write_record(&[h.to_string(), g.to_string()])?;
                }
                wtr.flush()?;
            } else {
                // fallback to stdout
                println!("lag,semivar");
                for (&h,&g) in lags.iter().zip(&semis) {
                    println!("{},{}", h, g);
                }
            }

            // 4) Optionally fit & plot
            if let Some(p) = plot {
                if let Some(m) = model {
                    let model_semis = variogram::fit_model(&lags, &semis, &m);
                    plot_variogram(&lags, &semis, Some(&lags), Some(&model_semis), &p)?;
                } else {
                    plot_variogram(&lags, &semis, None, None, &p)?;
                }
            }
        }

        Commands::Sgsim { data, sim, nsim, seed, out_prefix } => {
            if let Some(s) = seed {
                gslib_rs::utils::set_seed(s);
            }

            let (ds_data, data_coords, data_vals) = read_data(&data)?;
            let ds_sim = read_checked(&sim, 3, "simulation-target")?;
            let sim_coords = coords_of(&ds_sim);

            // Global mean of the conditioning values
            let mean = data_vals.iter().sum::<f64>() / (data_vals.len() as f64);

            // Default covariance function: exponential covariance with unit sill
            let cov_func = |h: f64| (-h).exp();

            let realizations = sgsim(&data_coords, &data_vals, &sim_coords, cov_func, mean, nsim);

            for (i, sim_vals) in realizations.into_iter().enumerate() {
                write_realization(
                    &ds_sim,
                    &ds_data.title,
                    format!("sim_{i}"),
                    sim_vals,
                    realization_path(&out_prefix, i),
                )?;
            }
        }

        Commands::Sisim { data, sim, threshold, nsim, out_prefix } => {
            let (ds_data, data_coords, data_vals) = read_data(&data)?;
            let ds_sim = read_checked(&sim, 3, "simulation-target")?;
            let sim_coords = coords_of(&ds_sim);

            // Default covariance model: exponential with unit sill
            let cov_func = |h: f64| (-h).exp();

            let realizations =
                sisim(&data_coords, &data_vals, &sim_coords, threshold, cov_func, nsim);

            for (i, sim_inds) in realizations.into_iter().enumerate() {
                write_realization(
                    &ds_sim,
                    &ds_data.title,
                    format!("ind_{i}"),
                    sim_inds.into_iter().map(f64::from).collect(),
                    realization_path(&out_prefix, i),
                )?;
            }
        }

        Commands::Sasim { .. } => {
            // TODO: parse SasimParams from flags, call `sasim`, write realizations.
            return Err("the sasim subcommand is not implemented yet".into());
        }
    }

    Ok(())
}
