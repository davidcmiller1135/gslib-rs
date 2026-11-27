use clap::{Parser, Subcommand};
use std::{error::Error, path::PathBuf};
use ndarray::Array2;


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

fn main() -> Result<(), Box<dyn Error>> {
    let cli = Cli::parse();

    match cli.cmd {
        Commands::Variogram { input, max_dist, n_lags, out_csv, model, plot } => {
            // 1) Read data
            let ds = read_gslib(&input)?;
            // assume columns: x,y,z,var
            let coords = ds.data.rows()
                .into_iter()
                .map(|r| (r[0], r[1], r[2]))
                .collect::<Vec<_>>();
            let values = ds.data.column(3).to_vec();

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
            // 1) Seed the RNG if requested
            if let Some(s) = seed {
                // fully qualified so you don't need to import utils
                gslib_rs::utils::set_seed(s);
            }

            // 2) Read conditioning data
            let ds_data = read_gslib(&data)?;
            let data_coords: Vec<(f64, f64, f64)> = ds_data
                .data
                .rows()
                .into_iter()
                .map(|r| (r[0], r[1], r[2]))
                .collect();
            let data_vals: Vec<f64> = ds_data.data.column(3).to_vec();

            // 3) Read simulation target locations
            let ds_sim = read_gslib(&sim)?;
            let sim_coords: Vec<(f64, f64, f64)> = ds_sim
                .data
                .rows()
                .into_iter()
                .map(|r| (r[0], r[1], r[2]))
                .collect();

            // 4) Compute global mean of the conditioning values
            let mean = data_vals.iter().sum::<f64>() / (data_vals.len() as f64);

            // 5) Default covariance function: exponential covariance with unit sill
            let cov_func = |h: f64| (-h).exp();

            // 6) Run the simulation
            let realizations = sgsim(
                &data_coords,
                &data_vals,
                &sim_coords,
                cov_func,
                mean,
                nsim,
            );

            // 7) Write out each realization as a GSLIB file
            for (i, sim_vals) in realizations.into_iter().enumerate() {
                // Build new var_names = original coords names + "sim_i"
                let mut var_names = ds_sim.var_names.clone();
                var_names.push(format!("sim_{}", i));

                let nrec = sim_coords.len();
                let nvar = var_names.len();

                // Build data array: [x,y,z,sim_val]
                let mut arr = Array2::<f64>::zeros((nrec, nvar));
                for (j, &(x, y, z)) in sim_coords.iter().enumerate() {
                    arr[(j, 0)] = x;
                    arr[(j, 1)] = y;
                    arr[(j, 2)] = z;
                    arr[(j, 3)] = sim_vals[j];
                }

                let out_ds = GslibData {
                    title: ds_data.title.clone(),
                    var_names: var_names.clone(),
                    data: arr,
                };

                let prefix = out_prefix.to_string_lossy();
                let out_path = std::path::PathBuf::from(format!("{}_{}.gslib", prefix, i));
                write_gslib(out_path, &out_ds)?;
            }
        }


        Commands::Sisim { data, sim, threshold, nsim, out_prefix } => {
            // 1) Read conditioning data
            let ds_data = read_gslib(&data)?;
            let data_coords: Vec<(f64, f64, f64)> = ds_data
                .data
                .rows()
                .into_iter()
                .map(|r| (r[0], r[1], r[2]))
                .collect();
            let data_vals: Vec<f64> = ds_data.data.column(3).to_vec();

            // 2) Read simulation targets (x,y,z only)
            let ds_sim = read_gslib(&sim)?;
            let sim_coords: Vec<(f64, f64, f64)> = ds_sim
                .data
                .rows()
                .into_iter()
                .map(|r| (r[0], r[1], r[2]))
                .collect();

            // 3) Default covariance model: exponential with unit sill
            let cov_func = |h: f64| (-h).exp();

            // 4) Run SISIM
            let realizations: Vec<Vec<u8>> = sisim(
                &data_coords,
                &data_vals,
                &sim_coords,
                threshold,
                cov_func,
                nsim,
            );

            // 5) Write out each realization as a GSLIB file
            for (i, sim_inds) in realizations.into_iter().enumerate() {
                // take original var_names (["X","Y","Z"]) and append indicator column
                let mut var_names = ds_sim.var_names.clone();
                var_names.push(format!("ind_{}", i));

                let nrec = sim_coords.len();
                let nvar = var_names.len();

                // Build an Array2<f64> with indicator values cast to f64
                let mut arr = Array2::<f64>::zeros((nrec, nvar));
                for (j, &(x, y, z)) in sim_coords.iter().enumerate() {
                    arr[(j, 0)] = x;
                    arr[(j, 1)] = y;
                    arr[(j, 2)] = z;
                    arr[(j, 3)] = sim_inds[j] as f64;
                }

                let out_ds = GslibData {
                    title: ds_data.title.clone(),
                    var_names: var_names.clone(),
                    data: arr,
                };

                let prefix = out_prefix.to_string_lossy();
                let out_path = std::path::PathBuf::from(format!("{}_{}.gslib", prefix, i));
                write_gslib(out_path, &out_ds)?;
            }
        }



        Commands::Sasim { .. } => {
            // TODO: parse SasimParams from flags, call `sasim`, write realizations.
            return Err("the sasim subcommand is not implemented yet".into());
        }
    }

    Ok(())
}
