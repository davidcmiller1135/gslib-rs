# gslib-rs

An early-stage Rust port of [GSLIB](http://www.gslib.com/) geostatistics routines, as a library
(`gslib_rs`) and a command-line tool (`cli`).

> **Status: prototype.** APIs will change, and results have not yet been validated against the
> original Fortran GSLIB. Do not use for production estimates.

## Workspace

| Crate | Contents |
|---|---|
| `gslib_rs` | GSLIB/Geo-EAS file I/O, empirical variograms and spherical/exponential models, ordinary/simple kriging and simple cokriging, normal-score and back transforms, coordinate rotation, SGS, SIS, SASIM (partial), plotting |
| `cli` | `gslib-rs` binary with `variogram`, `sgsim` and `sisim` subcommands |

Not implemented yet: declustering, annealing calibration, the `sasim` CLI subcommand,
anisotropy, nested/nugget variogram models, and neighbourhood search for simulation
(SGS/SIS currently krige with every datum and previously simulated node, so keep grids small).

## Usage

```sh
cargo build --release

# Empirical variogram to stdout (input columns: x y z value)
cargo run -p cli -- variogram --input data.dat --max-dist 50 --n-lags 10

# ...with a CSV and a plot with a model overlay
cargo run -p cli -- variogram -i data.dat -d 50 -n 10 -o vario.csv --model spherical --plot vario.png

# Sequential Gaussian / indicator simulation
cargo run -p cli -- sgsim --data data.dat --sim grid.dat --nsim 5 --seed 42 --out-prefix sgs
cargo run -p cli -- sisim --data data.dat --sim grid.dat --threshold 1.5 --nsim 5 --out-prefix sis
```

Data files use the simplified Geo-EAS format: a title line, the number of variables, one
variable name per line, then whitespace-separated rows. The CLI expects the first three
columns to be x, y, z and (for data files) the fourth to be the value.

## Development

```sh
cargo test --workspace
cargo clippy --workspace --all-targets -- -D warnings
```

## License

MIT, see [LICENSE](LICENSE).
