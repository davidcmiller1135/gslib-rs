use clap::{Parser, Subcommand};

/// GSLIB-RS command-line interface
#[derive(Parser)]
#[command(name = "gslib-rs")]
#[command(author, version, about)]
struct Cli {
    #[command(subcommand)]
    cmd: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Compute an empirical variogram
    Variogram {
        #[arg(long)]
        input: String,
        #[arg(long)]
        output: String,
    },
    /// Run sequential Gaussian simulation
    Sgsim {
        #[arg(long)]
        input: String,
        #[arg(long)]
        nsim: usize,
        #[arg(long)]
        output: String,
    },
    // … add other subcommands …
}

fn main() {
    let cli = Cli::parse();
    match cli.cmd {
        Commands::Variogram { input, output } => {
            // call gslib_rs::variogram::compute_empirical + plot, etc.
            unimplemented!();
        }
        Commands::Sgsim { input, nsim, output } => {
            unimplemented!();
        }
    }
}
