# gslib-rs

A Rust port of [GSLIB](http://www.gslib.com/), the Geostatistical Software
Library: variograms, kriging, and stochastic simulation. It is a library
(`gslib_rs`) plus a command-line tool (`gslib-rs`).

> **Status: early.** The core is built and tested, but it covers only part of
> what GSLIB does and has not yet been checked against the Fortran.
> [`docs/ROADMAP.md`](docs/ROADMAP.md) compares it with GSLIB program by program.

## What works

- **Variogram models:** nugget plus nested spherical, exponential, Gaussian,
  power and hole-effect structures, with GSLIB three-angle anisotropy.
- **Experimental variograms:** omnidirectional, reporting pairs and mean
  distance per lag.
- **Kriging:** simple, ordinary (also with the power model) and simple
  cokriging, each returning the kriging variance.
- **Search neighbourhood:** GSLIB's anisotropic search ellipsoid with `ndmin`/`ndmax` and octants (the CLI exposes the nearest `--max-samples` within `--radius`).
- **Transforms:** normal scores (weighted, ties handled) and back-transform
  with GSLIB tail options.
- **Simulation:** sequential Gaussian (SGS), sequential indicator (SIS) and
  simulated annealing (SASIM). Each takes an explicit, seedable random number
  generator.
- **GSLIB file I/O:** tolerant of blank lines, grid headers and Fortran `D`
  exponents; parse errors report line numbers.
- **Plots (PNG/SVG):** variogram, scatter, histogram and pixel map. These sit
  behind the default `plot` feature.

## Command line

```sh
cargo build --release

# Experimental variogram (columns default to x y z value = 1 2 3 4),
# plotted against a model: nugget 0.1 + spherical (sill 0.9, range 40).
gslib-rs variogram -i data.dat -d 100 -n 10 -o vario.csv --plot vario.png \
    --nugget 0.1 --structure sph:0.9:40

# Ordinary kriging at target locations, with an anisotropic model
# (major range 40 at azimuth 45, minor 20, vertical 5).
gslib-rs krige --data data.dat --targets grid.dat --out est.dat \
    --nugget 0.1 --structure sph:0.9:40:20:5:45:0:0 --max-samples 16

# Sequential Gaussian simulation. Data are normal-scored and results
# back-transformed, so give the model in normal-score units (sill 1).
gslib-rs sgsim --data data.dat --grid-x 50:0.5:1 --grid-y 50:0.5:1 \
    --nsim 10 --seed 69069 --max-nodes 12 \
    --out-prefix sgs --structure exp:1:50

# Sequential indicator simulation at one threshold.
gslib-rs sisim --data data.dat --grid-x 50:0.5:1 --grid-y 50:0.5:1 --threshold 2.5 \
    --out-prefix sis --structure sph:0.25:30
```

`--structure` takes `TYPE:CC:A_HMAX[:A_HMIN:A_VERT[:AZIMUTH:DIP:PLUNGE]]`,
where `TYPE` is `sph`, `exp`, `gau` or `hol`, or takes `pow:CC:OMEGA`. These
are the same numbers as a GSLIB parameter file. `--grid-x`, `--grid-y` and
`--grid-z` take `N:MIN:SIZ`, GSLIB's `nx xmn xsiz` (y and z default to one
node); simulation output lists X, Y, Z and one value column per node, x
fastest. Run `gslib-rs <command> --help` for every option.

**Missing values.** No trimming is applied unless you ask for it, but a
command refuses data whose value column contains −999, −99 or a value of
magnitude 10²¹ or more, since GSLIB files use these as missing-value codes.
Give `--tmin`/`--tmax` to trim them (a value is kept if `tmin ≤ v < tmax`,
so `--tmin -998` drops −999), or `--no-trimming` if they are real data. NaN
and infinite values are errors unless `--nan missing` or
`--inf missing|clamp` says otherwise. Coordinates are never trimmed.

**Search.** `--min-samples`/`--max-samples` (GSLIB `ndmin`/`ndmax`),
`--per-octant`, and an anisotropic ellipsoid from `--radius`,
`--radius-minor`, `--radius-vertical` and `--search-angles AZ:DIP:PLUNGE`.
`krige --on-too-few unestimated` writes `--unest` (default −999) where a
target has too few data, as `kt3d` does. `sgsim`/`sisim` take
`--max-nodes`, `--assign-data` (GSLIB's two-part search, with
`--on-shared-node` and `--outside-grid`) and `--node-table gslib` for
GSLIB's node-table limits.

**Seeds.** `--seed` is a master seed: realization i draws from its own
stream (ChaCha8, stream i), so it is the same whatever `--nsim` is and in
whatever order realizations run.

**Threads.** `krige`, `variogram`, `sgsim` and `sisim` use every core (the
simulations run realizations in parallel); set `RAYON_NUM_THREADS` to limit
them.
Results are identical for any number of threads.

## Library

```rust
use gslib_rs::geometry::Point3;
use gslib_rs::kriging::ordinary_kriging;
use gslib_rs::model::{Kind, Structure, VariogramModel};
use gslib_rs::sample::Sample;

let model = VariogramModel::new(0.1, vec![Structure::isotropic(Kind::Spherical, 0.9, 10.0)?])?;
let data = [
    Sample::new(Point3::new(0.0, 0.0, 0.0), 1.0)?,
    Sample::new(Point3::new(4.0, 0.0, 0.0), 3.0)?,
];
let estimate = ordinary_kriging(&data, Point3::new(2.0, 0.0, 0.0), &model)?;
println!("{} ± {}", estimate.value, estimate.variance.sqrt());
```

## Engineering standards

This project is deliberately stricter than most Rust code:

- **Lints:** every clippy group is denied (`all`, `pedantic`, `nursery`,
  `cargo`), plus restriction lints:
  - library code cannot `unwrap`, `expect`, `panic`, or index with `[]`;
  - there are no `as` casts outside [`num.rs`](gslib_rs/src/num.rs);
  - every exception is an `#[expect(lint, reason = "…")]`.

  There are currently two exceptions in the whole codebase.
- **Errors:** library functions return `Result` and never panic on bad input.
- **Toolchain:** pinned in `rust-toolchain.toml` (so new lints arrive only
  deliberately). The MSRV is 1.89 and is tested in CI.
- **Supply chain:** `cargo deny` checks advisories, licences, duplicate crate
  versions and crate sources.
- **Tests:** unit, property-based (`proptest`) and end-to-end CLI tests.

Before pushing, run:

```sh
cargo fmt --all --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
cargo deny check
```

## License

MIT. See [LICENSE](LICENSE).

The GSLIB 2.0 source is vendored unmodified in
[`vendor/gslib-2.0/`](vendor/gslib-2.0/README.md) as the parity reference.
Parts of gslib-rs are derived from it; its 1996 Stanford copyright notice is
reproduced in [NOTICE](NOTICE).
