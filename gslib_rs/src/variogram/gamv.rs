//! The `gamv` kernel (`docs/spec/gamv.md`, D-051).
//!
//! # Pair enumeration
//!
//! Pairs are enumerated exactly as for [`experimental_variogram`]
//! (D-031, D-032, D-034): rows in fixed blocks of 256, each row's later
//! partners from the super-block search within the lags' cut-off, in block
//! order, and the block sums combined in order. Each pair `i < j` is visited
//! once, tested against every direction, and added to every variogram and
//! every lag class it falls in. Pair acceptance does not depend on summation
//! order, so pair counts match GSLIB's exactly; the sums differ from GSLIB's
//! double loop only in the last bits.
//!
//! There are no self-pairs, and each pair is reported once (D-017).
//!
//! # Tail and head
//!
//! For an accepted pair the [`PairOrder`] decides which datum is the tail
//! (GSLIB ll. 591–613). The tail value comes from the variogram's `tail`
//! variable at the tail datum, and the head value from its `head` variable at
//! the head datum.
//!
//! # Omnidirectional directions (D-043)
//!
//! When `atol ≥ 90`, each accepted pair also adds its other ordering to the
//! sums (tail and head data swapped), except for the cross-semivariogram, as
//! GSLIB does. Values therefore equal GSLIB's and do not depend on the data's
//! row order. [`GamvLag::pairs`] still counts each pair once, and
//! [`GamvLag::orderings`] counts the orderings summed (GSLIB's `np`).
//!
//! # Missing values
//!
//! A pair contributes to a variogram only if every value it uses is present.
//! Trimming happens when values are extracted (D-019, D-048); here values are
//! already `Option`.
//!
//! [`experimental_variogram`]: super::experimental_variogram

use crate::error::{Error, Result, check_lengths};
use crate::geometry::{PiConstant, Point3};
use crate::num::to_f64;

use super::directional::{Direction, DirectionWindow, EPSLON, PairOrder, Separation};
use super::{CentredLags, LagScheme, PRUNE_MARGIN, ROWS_PER_BLOCK};
use crate::search::Blocks;

/// Coordinates and variables for [`gamv`].
#[derive(Debug, Clone, PartialEq)]
pub struct Data {
    points: Vec<Point3>,
    columns: Vec<Vec<Option<f64>>>,
}

impl Data {
    /// Data at `points`, with one column per variable, each holding one value
    /// (or `None`) per point.
    ///
    /// # Errors
    /// - [`Error::InvalidInput`] if a coordinate or a present value is not
    ///   finite.
    /// - [`Error::LengthMismatch`] if a column's length differs from the
    ///   number of points.
    pub fn new(points: Vec<Point3>, columns: Vec<Vec<Option<f64>>>) -> Result<Self> {
        if let Some(bad) = points.iter().position(|p| !p.is_finite()) {
            return Err(Error::invalid(
                "data coordinates",
                format!("point {bad} is not finite"),
            ));
        }
        for (k, column) in columns.iter().enumerate() {
            check_lengths("data column and points", column.len(), points.len())?;
            if let Some(bad) = column
                .iter()
                .position(|v| v.is_some_and(|v| !v.is_finite()))
            {
                return Err(Error::invalid(
                    "data values",
                    format!("variable {k}, point {bad} is not finite"),
                ));
            }
        }
        Ok(Self { points, columns })
    }

    /// Number of data points.
    #[must_use]
    pub const fn len(&self) -> usize {
        self.points.len()
    }

    /// True if there are no data points.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.points.is_empty()
    }

    /// Number of variables.
    #[must_use]
    pub const fn variables(&self) -> usize {
        self.columns.len()
    }
}

/// Which measure a variogram computes (GSLIB's `ivtype`, `GAMV.FOR` ll.
/// 404–414 and 628–871).
///
/// Below, `t` is a tail value and `h` a head value of a pair, `t'` and `h'`
/// are the tail and head variables at the other datum, `N` is the number of
/// pairs and `m_t`, `m_h` are the lag's tail and head means.
///
/// Indicator variograms (GSLIB types 9 and 10) are semivariograms of derived
/// indicator variables (D-051).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VariogramKind {
    /// Type 1: `Σ (t − h)² / (2N)`.
    Semivariogram,
    /// Type 2: `Σ (t' − t)(h − h') / (2N)`. Needs all four values. Never
    /// doubled for omnidirectional directions.
    CrossSemivariogram,
    /// Type 3: `Σ t·h / N − m_t·m_h`.
    Covariance,
    /// Type −3: `√(sill_tail · sill_head)` minus the type-3 covariance.
    /// `None` where a sill is unavailable.
    CovarianceAsVariogram,
    /// Type 4: the covariance divided by the tail and head standard
    /// deviations of the lag. `None` where their product is below 1e-20.
    Correlogram,
    /// Type 5: `Σ (t − h)² / N` divided by `((m_t + m_h)/2)²`. GSLIB applies
    /// no factor of 1/2. `None` where the squared mean is below 1e-20.
    GeneralRelative,
    /// Type 6: `Σ (2(h − t)/(h + t))² / (2N)`, over pairs with `|h + t| >
    /// 1e-20`. Other pairs are counted as skipped.
    PairwiseRelative,
    /// Type 7: `Σ (ln h − ln t)² / (2N)`, over pairs with both values `>
    /// 1e-20`. Other pairs are counted as skipped.
    Logarithmic,
    /// Type 8: `Σ |t − h| / (2N)`.
    Madogram,
}

impl VariogramKind {
    /// GSLIB's `ivtype` (1 to 8, or −3).
    #[must_use]
    pub const fn gslib_type(self) -> i32 {
        match self {
            Self::Semivariogram => 1,
            Self::CrossSemivariogram => 2,
            Self::Covariance => 3,
            Self::CovarianceAsVariogram => -3,
            Self::Correlogram => 4,
            Self::GeneralRelative => 5,
            Self::PairwiseRelative => 6,
            Self::Logarithmic => 7,
            Self::Madogram => 8,
        }
    }
}

/// One variogram to compute: tail and head variables (indices into
/// [`Data`]'s columns) and the measure.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct VariogramSpec {
    /// Variable for the tail values.
    pub tail: usize,
    /// Variable for the head values.
    pub head: usize,
    /// The measure.
    pub kind: VariogramKind,
    /// Divide by the variable's sill ([`GamvRequest::sills`]). GSLIB's
    /// `isill = 1`, which applies only to semivariograms with tail = head;
    /// requesting it for anything else is an error rather than ignored
    /// (D-049).
    pub standardize: bool,
}

/// The second, swapped contribution of omnidirectional pairwise-relative and
/// logarithmic variograms (D-047).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum OmniRelativeLog {
    /// Computed from the swapped values, like every other measure. The
    /// default.
    #[default]
    Swapped,
    /// GSLIB's defect (`GAMV.FOR` ll. 730, 758): guarded by the swapped
    /// values but computed from the unswapped ones. Legacy policy, for
    /// parity. A non-finite result is still an error.
    GslibDuplicate,
}

/// The legacy policies that act inside the `gamv` kernel.
///
/// `Default` is gslib-rs's behaviour; [`GamvPolicies::gslib`] reproduces
/// GSLIB. Policies that act on input or output (trimming, indicators,
/// undefined values, non-positive tolerances, which sills to standardise
/// with) belong to the adapters (D-051). There is no self-pair policy
/// (D-017).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct GamvPolicies {
    /// π in angle conversions (D-046).
    pub pi: PiConstant,
    /// Omnidirectional types 6 and 7 (D-047).
    pub omni_relative_log: OmniRelativeLog,
}

impl GamvPolicies {
    /// Every legacy policy: GSLIB's behaviour.
    #[must_use]
    pub const fn gslib() -> Self {
        Self {
            pi: PiConstant::Gslib,
            omni_relative_log: OmniRelativeLog::GslibDuplicate,
        }
    }
}

/// What [`gamv`] computes.
#[derive(Debug, Clone, PartialEq)]
pub struct GamvRequest {
    /// Lag classes (GSLIB's `nlag`, `xlag`, `xltol`).
    pub lags: CentredLags,
    /// Directions (at least one).
    pub directions: Vec<Direction>,
    /// Variograms (at least one).
    pub variograms: Vec<VariogramSpec>,
    /// One sill (variance) per variable, `None` where unavailable, or empty
    /// if no variogram needs one. Used for standardisation and by
    /// [`VariogramKind::CovarianceAsVariogram`]. GSLIB uses each variable's
    /// population variance over its valid data (D-048).
    pub sills: Vec<Option<f64>>,
    /// Kernel policies.
    pub policies: GamvPolicies,
}

/// One lag of one directional variogram.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub struct GamvLag {
    /// Mean separation of the pairs (`None` if there are none).
    pub mean_distance: Option<f64>,
    /// The variogram value: `None` if there are no pairs, or where the
    /// measure is undefined (see [`VariogramKind`]).
    pub gamma: Option<f64>,
    /// Number of pairs, each counted once (D-017).
    pub pairs: usize,
    /// Number of orderings summed: equal to `pairs`, except in an
    /// omnidirectional direction, where each pair adds both orderings that
    /// have the values they need and lie in the measure's domain (D-043).
    /// This is GSLIB's `np`.
    pub orderings: usize,
    /// Mean of the tail values (`None` if there are no pairs).
    pub tail_mean: Option<f64>,
    /// Mean of the head values (`None` if there are no pairs).
    pub head_mean: Option<f64>,
    /// Correlogram only: variance of the tail values.
    pub tail_variance: Option<f64>,
    /// Correlogram only: variance of the head values.
    pub head_variance: Option<f64>,
    /// Pairwise relative and logarithmic only: orderings outside the
    /// measure's domain, not counted in `orderings`.
    pub skipped: usize,
}

/// One variogram in one direction: `nlag + 2` lags.
#[derive(Debug, Clone, PartialEq)]
pub struct DirectionalVariogram {
    /// Index into [`GamvRequest::variograms`].
    pub variogram: usize,
    /// Index into [`GamvRequest::directions`].
    pub direction: usize,
    /// The lags, indexed as in [`CentredLags`].
    pub lags: Vec<GamvLag>,
}

/// One contribution to a lag: GSLIB's increments to `gam`, `hm`, `tm`, `hv`
/// and `tv`.
#[derive(Debug, Clone, Copy)]
enum Increment {
    Add {
        gamma: f64,
        tail: f64,
        head: f64,
        tail_sq: f64,
        head_sq: f64,
    },
    Skip,
}

/// Running sums for one lag (`np`, `dis`, `gam`, `hm`, `tm`, `hv`, `tv` in
/// GSLIB).
#[derive(Debug, Clone, Copy, Default)]
struct Sums {
    pairs: usize,
    orderings: usize,
    skipped: usize,
    distance: f64,
    gamma: f64,
    tail: f64,
    head: f64,
    tail_sq: f64,
    head_sq: f64,
}

/// Sills resolved for one variogram.
#[derive(Debug, Clone, Copy)]
struct VariogramSills {
    standardize: Option<f64>,
    tail: Option<f64>,
    head: Option<f64>,
}

impl Sums {
    /// Add another block's sums (D-031: blocks are combined in order).
    fn combine(&mut self, other: &Self) {
        self.pairs += other.pairs;
        self.orderings += other.orderings;
        self.skipped += other.skipped;
        self.distance += other.distance;
        self.gamma += other.gamma;
        self.tail += other.tail;
        self.head += other.head;
        self.tail_sq += other.tail_sq;
        self.head_sq += other.head_sq;
    }

    /// Add one pair's increments (one or two orderings).
    fn add_pair(&mut self, h: f64, incs: &[Option<Increment>; 2]) {
        let before = self.orderings;
        for inc in incs.iter().flatten() {
            self.add(h, *inc);
        }
        if self.orderings > before {
            self.pairs += 1;
        }
    }

    fn add(&mut self, h: f64, inc: Increment) {
        match inc {
            Increment::Add {
                gamma,
                tail,
                head,
                tail_sq,
                head_sq,
            } => {
                self.orderings += 1;
                self.distance += h;
                self.gamma += gamma;
                self.tail += tail;
                self.head += head;
                self.tail_sq += tail_sq;
                self.head_sq += head_sq;
            }
            Increment::Skip => self.skipped += 1,
        }
    }

    /// Averages and the measure (`GAMV.FOR` ll. 798–872), in GSLIB's order
    /// of operations.
    #[expect(
        clippy::suboptimal_flops,
        clippy::suspicious_operation_groupings,
        reason = "GSLIB's unfused formulas, term for term (`sq/n − mean²` is the variance)"
    )]
    fn finish(self, kind: VariogramKind, sills: VariogramSills) -> Result<GamvLag> {
        let mut lag = GamvLag {
            mean_distance: None,
            gamma: None,
            pairs: self.pairs,
            orderings: self.orderings,
            tail_mean: None,
            head_mean: None,
            tail_variance: None,
            head_variance: None,
            skipped: self.skipped,
        };
        if self.orderings == 0 {
            return Ok(lag);
        }
        let n = to_f64(self.orderings);
        let mut gam = self.gamma / n;
        let hm = self.tail / n;
        let tm = self.head / n;
        if let Some(sill) = sills.standardize {
            gam /= sill;
        }
        let gamma = match kind {
            VariogramKind::Semivariogram
            | VariogramKind::CrossSemivariogram
            | VariogramKind::PairwiseRelative
            | VariogramKind::Logarithmic
            | VariogramKind::Madogram => Some(0.5 * gam),
            VariogramKind::Covariance => Some(gam - hm * tm),
            VariogramKind::CovarianceAsVariogram => sills
                .tail
                .zip(sills.head)
                .map(|(st, sh)| st.sqrt() * sh.sqrt() - (gam - hm * tm)),
            VariogramKind::Correlogram => {
                let sd = |sq: f64, mean: f64| (sq / n - mean * mean).max(0.0).sqrt();
                let (hv, tv) = (sd(self.tail_sq, hm), sd(self.head_sq, tm));
                lag.tail_variance = Some(hv * hv);
                lag.head_variance = Some(tv * tv);
                (hv * tv >= EPSLON).then(|| (gam - hm * tm) / (hv * tv))
            }
            VariogramKind::GeneralRelative => {
                let mean = 0.5 * (hm + tm);
                let squared = mean * mean;
                (squared >= EPSLON).then(|| gam / squared)
            }
        };
        lag.mean_distance = Some(self.distance / n);
        lag.tail_mean = Some(hm);
        lag.head_mean = Some(tm);
        lag.gamma = gamma;
        let finite = [
            lag.gamma,
            lag.tail_mean,
            lag.head_mean,
            lag.tail_variance,
            lag.head_variance,
        ]
        .into_iter()
        .flatten()
        .all(f64::is_finite);
        if finite {
            Ok(lag)
        } else {
            Err(Error::invalid(
                "variogram data",
                format!(
                    "a {kind:?} lag is not finite (values too large, or the \
                     legacy types 6/7 policy met an undefined pair)"
                ),
            ))
        }
    }
}

/// Tail and head values for one variogram and one pair, in GSLIB's naming:
/// `vrh`/`vrt` for the pair as ordered, `vrhpr`/`vrtpr` with the roles of the
/// two data swapped.
#[derive(Debug, Clone, Copy)]
struct PairValues {
    tail: Option<f64>,
    head: Option<f64>,
    tail_swapped: Option<f64>,
    head_swapped: Option<f64>,
}

fn value(columns: &[Vec<Option<f64>>], variable: usize, datum: usize) -> Option<f64> {
    columns.get(variable)?.get(datum).copied().flatten()
}

fn pair_values(
    columns: &[Vec<Option<f64>>],
    spec: VariogramSpec,
    (first, second): (usize, usize),
    order: PairOrder,
) -> PairValues {
    let (t, h) = match order {
        PairOrder::AsGiven => (first, second),
        PairOrder::Swapped => (second, first),
    };
    PairValues {
        tail: value(columns, spec.tail, t),
        head: value(columns, spec.head, h),
        tail_swapped: value(columns, spec.tail, h),
        head_swapped: value(columns, spec.head, t),
    }
}

fn plain(gamma: f64, tail: f64, head: f64, with_squares: bool) -> Increment {
    let (tail_sq, head_sq) = if with_squares {
        (tail * tail, head * head)
    } else {
        (0.0, 0.0)
    };
    Increment::Add {
        gamma,
        tail,
        head,
        tail_sq,
        head_sq,
    }
}

/// The measure-specific increment for values `(t, h)`; for types 6 and 7,
/// `(ft, fh)` are the values the formula uses (they differ from `(t, h)` only
/// under the legacy omnidirectional policy).
fn increment(kind: VariogramKind, (t, h): (f64, f64), (ft, fh): (f64, f64)) -> Increment {
    match kind {
        VariogramKind::Semivariogram | VariogramKind::GeneralRelative => {
            let d = t - h;
            plain(d * d, t, h, false)
        }
        VariogramKind::Covariance | VariogramKind::CovarianceAsVariogram => {
            plain(t * h, t, h, false)
        }
        VariogramKind::Correlogram => plain(t * h, t, h, true),
        VariogramKind::Madogram => plain((t - h).abs(), t, h, false),
        VariogramKind::PairwiseRelative => {
            if (h + t).abs() > EPSLON {
                let g = 2.0 * (fh - ft) / (fh + ft);
                plain(g * g, t, h, false)
            } else {
                Increment::Skip
            }
        }
        VariogramKind::Logarithmic => {
            if h > EPSLON && t > EPSLON {
                let g = fh.ln() - ft.ln();
                plain(g * g, t, h, false)
            } else {
                Increment::Skip
            }
        }
        // Handled by `increments`: it needs all four values.
        VariogramKind::CrossSemivariogram => Increment::Skip,
    }
}

/// The (one or two) increments a pair makes to a variogram, or none if a
/// value it needs is missing (`GAMV.FOR` ll. 615–785).
fn increments(
    kind: VariogramKind,
    values: PairValues,
    omni: bool,
    policy: OmniRelativeLog,
) -> [Option<Increment>; 2] {
    let (Some(t), Some(h)) = (values.tail, values.head) else {
        return [None, None];
    };
    let swapped = values.tail_swapped.zip(values.head_swapped);
    if kind == VariogramKind::CrossSemivariogram {
        return [
            swapped.map(|(tp, hp)| {
                // vrh = t, vrt = h, vrhpr = tp, vrtpr = hp.
                plain((tp - t) * (h - hp), 0.5 * (t + tp), 0.5 * (h + hp), false)
            }),
            None,
        ];
    }
    let first = increment(kind, (t, h), (t, h));
    let second = if omni {
        swapped.map(|pair| {
            let formula = match policy {
                OmniRelativeLog::Swapped => pair,
                OmniRelativeLog::GslibDuplicate => (t, h),
            };
            increment(kind, pair, formula)
        })
    } else {
        None
    };
    [Some(first), second]
}

fn validate(data: &Data, request: &GamvRequest) -> Result<Vec<VariogramSills>> {
    if request.directions.is_empty() {
        return Err(Error::invalid("directions", "at least one is required"));
    }
    if request.variograms.is_empty() {
        return Err(Error::invalid("variograms", "at least one is required"));
    }
    if !request.sills.is_empty() {
        check_lengths("sills and variables", request.sills.len(), data.variables())?;
    }
    if let Some(bad) = request.sills.iter().flatten().find(|s| !s.is_finite()) {
        return Err(Error::invalid("sill", format!("must be finite, got {bad}")));
    }
    let sill = |v: usize| request.sills.get(v).copied().flatten();
    request
        .variograms
        .iter()
        .enumerate()
        .map(|(k, spec)| {
            for (what, v) in [("tail", spec.tail), ("head", spec.head)] {
                if v >= data.variables() {
                    return Err(Error::invalid(
                        "variogram variable",
                        format!(
                            "variogram {k}: {what} variable {v} does not exist ({} variables)",
                            data.variables()
                        ),
                    ));
                }
            }
            let standardize = if spec.standardize {
                if spec.kind != VariogramKind::Semivariogram || spec.tail != spec.head {
                    return Err(Error::invalid(
                        "standardisation",
                        format!(
                            "variogram {k}: only a semivariogram with tail = head can be \
                             standardised (GSLIB ignores the request)"
                        ),
                    ));
                }
                match sill(spec.tail) {
                    Some(s) if s > 0.0 => Some(s),
                    other => {
                        return Err(Error::invalid(
                            "standardisation",
                            format!(
                                "variogram {k}: the sill of variable {} must be positive, \
                                 got {other:?} (GSLIB silently skips it)",
                                spec.tail
                            ),
                        ));
                    }
                }
            } else {
                None
            };
            // GSLIB marks an unavailable sill as negative (−999).
            let usable = |v: usize| sill(v).filter(|s| *s >= 0.0);
            Ok(VariogramSills {
                standardize,
                tail: usable(spec.tail),
                head: usable(spec.head),
            })
        })
        .collect()
}

/// Experimental variograms of scattered data, as GSLIB `gamv` computes them.
///
/// Returns one [`DirectionalVariogram`] per variogram and direction, ordered
/// by variogram and then direction (GSLIB's output order), each with
/// `nlag + 2` lags.
///
/// # Errors
/// [`Error::InvalidInput`] if:
/// - there are no directions or no variograms;
/// - a variogram names a variable that `data` does not have;
/// - sills are given for the wrong number of variables, or are not finite;
/// - standardisation is requested for anything but a semivariogram with
///   tail = head, or the sill is not positive;
/// - a result is not finite.
pub fn gamv(data: &Data, request: &GamvRequest) -> Result<Vec<DirectionalVariogram>> {
    let sills = validate(data, request)?;
    let windows: Vec<DirectionWindow> = request
        .directions
        .iter()
        .map(|d| d.window(request.policies.pi))
        .collect();
    let size = request.variograms.len() * windows.len() * request.lags.len();
    let blocks = Blocks::new(&data.points);
    // As for the semivariogram (D-032): the cut-off, widened so that no pair
    // at it is lost to rounding; the lag classes still decide every pair.
    let reach = request.lags.max_distance() * (1.0 + PRUNE_MARGIN);
    let plan = Plan {
        data,
        request,
        windows: &windows,
        blocks: &blocks,
        reach,
        size,
    };
    let mut sums = vec![Sums::default(); size];
    for block in 0..data.points.len().div_ceil(ROWS_PER_BLOCK) {
        for (t, b) in sums.iter_mut().zip(&plan.block_sums(block)) {
            t.combine(b);
        }
    }

    let nlags = request.lags.len();
    let mut out = Vec::with_capacity(request.variograms.len() * windows.len());
    let mut chunks = sums.chunks_exact(nlags);
    for ((v, spec), &sill) in request.variograms.iter().enumerate().zip(&sills) {
        for d in 0..windows.len() {
            let lags = chunks
                .next()
                .map(|chunk| {
                    chunk
                        .iter()
                        .map(|s| s.finish(spec.kind, sill))
                        .collect::<Result<Vec<_>>>()
                })
                .transpose()?
                .unwrap_or_default();
            out.push(DirectionalVariogram {
                variogram: v,
                direction: d,
                lags,
            });
        }
    }
    Ok(out)
}

/// What every row block needs.
struct Plan<'a> {
    data: &'a Data,
    request: &'a GamvRequest,
    windows: &'a [DirectionWindow],
    blocks: &'a Blocks,
    reach: f64,
    size: usize,
}

impl Plan<'_> {
    /// The sums from the pairs whose first datum is in row block `block`.
    fn block_sums(&self, block: usize) -> Vec<Sums> {
        let mut acc = vec![Sums::default(); self.size];
        let start = block.saturating_mul(ROWS_PER_BLOCK);
        let end = start
            .saturating_add(ROWS_PER_BLOCK)
            .min(self.data.points.len());
        for i in start..end {
            let Some(&p) = self.data.points.get(i) else {
                continue;
            };
            self.blocks.visit_within(p, self.reach, |j, q| {
                if j > i {
                    self.accumulate(&Separation::new(p, q), (i, j), &mut acc);
                }
            });
        }
        acc
    }

    /// Add one pair to every direction and variogram that accepts it.
    fn accumulate(&self, pair: &Separation, data_indices: (usize, usize), sums: &mut [Sums]) {
        let lag_range = self.request.lags.classes(pair.distance());
        if lag_range.is_empty() {
            return;
        }
        let nlags = self.request.lags.len();
        let ndir = self.windows.len();
        let h = pair.distance();
        for (d, window) in self.windows.iter().enumerate() {
            let Some(order) = window.accept(pair) else {
                continue;
            };
            for (v, &spec) in self.request.variograms.iter().enumerate() {
                let values = pair_values(&self.data.columns, spec, data_indices, order);
                let incs = increments(
                    spec.kind,
                    values,
                    window.is_omnidirectional(),
                    self.request.policies.omni_relative_log,
                );
                if incs[0].is_none() {
                    continue;
                }
                let base = (v * ndir + d) * nlags;
                for lag in lag_range.clone() {
                    if let Some(s) = sums.get_mut(base + lag) {
                        s.add_pair(h, &incs);
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    use proptest::prelude::*;

    fn p(x: f64, y: f64) -> Point3 {
        Point3::new(x, y, 0.0)
    }

    fn semivariogram(variable: usize) -> VariogramSpec {
        VariogramSpec {
            tail: variable,
            head: variable,
            kind: VariogramKind::Semivariogram,
            standardize: false,
        }
    }

    /// Azimuth `azm` ± `atol`, narrow bandwidths off, horizontal.
    fn direction(azm: f64, atol: f64) -> Direction {
        Direction::new(azm, atol, 1e9, 0.0, 90.0, 1e9).unwrap()
    }

    fn request(lags: CentredLags, directions: Vec<Direction>) -> GamvRequest {
        GamvRequest {
            lags,
            directions,
            variograms: vec![semivariogram(0)],
            sills: vec![],
            policies: GamvPolicies::default(),
        }
    }

    fn one_variable(points: Vec<Point3>, values: &[f64]) -> Data {
        Data::new(points, vec![values.iter().copied().map(Some).collect()]).unwrap()
    }

    #[test]
    fn invalid_inputs_are_errors() {
        let pts = vec![p(0.0, 0.0), p(1.0, 0.0)];
        assert!(
            Data::new(pts.clone(), vec![vec![Some(1.0)]]).is_err(),
            "short column"
        );
        assert!(
            Data::new(pts.clone(), vec![vec![Some(1.0), Some(f64::NAN)]]).is_err(),
            "NaN value"
        );
        assert!(
            Data::new(vec![p(f64::INFINITY, 0.0)], vec![]).is_err(),
            "infinite coordinate"
        );
        let data = one_variable(pts, &[1.0, 2.0]);
        let lags = CentredLags::new(2, 1.0, 0.5).unwrap();
        assert!(
            gamv(&data, &request(lags, vec![])).is_err(),
            "no directions"
        );
        let mut r = request(lags, vec![direction(0.0, 90.0)]);
        r.variograms.clear();
        assert!(gamv(&data, &r).is_err(), "no variograms");
        r.variograms.push(VariogramSpec {
            tail: 0,
            head: 1,
            kind: VariogramKind::Semivariogram,
            standardize: false,
        });
        assert!(gamv(&data, &r).is_err(), "head variable out of range");
    }

    #[test]
    fn output_is_ordered_by_variogram_then_direction() {
        let data = Data::new(
            vec![p(0.0, 0.0), p(1.0, 0.0)],
            vec![vec![Some(0.0), Some(1.0)], vec![Some(0.0), Some(3.0)]],
        )
        .unwrap();
        let mut r = request(
            CentredLags::new(3, 1.0, 0.5).unwrap(),
            vec![direction(0.0, 10.0), direction(90.0, 10.0)],
        );
        r.variograms.push(semivariogram(1));
        let out = gamv(&data, &r).unwrap();
        let order: Vec<_> = out.iter().map(|v| (v.variogram, v.direction)).collect();
        assert_eq!(order, vec![(0, 0), (0, 1), (1, 0), (1, 1)]);
        assert!(out.iter().all(|v| v.lags.len() == 5), "nlag + 2 lags each");
        // Only the east–west direction sees the pair, at lag k = 1 (index 2).
        assert_eq!(out[1].lags[2].gamma, Some(0.5));
        assert_eq!(out[3].lags[2].gamma, Some(4.5));
        assert_eq!(out[0].lags[2].pairs, 0, "north sees nothing");
    }

    /// z = x along a line: every pair at distance k has squared difference
    /// k², so γ(k) = k²/2 exactly, with 11 − k pairs.
    #[test]
    fn linear_trend_along_the_direction() {
        let pts: Vec<Point3> = (0..=10).map(|i| p(f64::from(i), 0.0)).collect();
        let vals: Vec<f64> = (0..=10).map(f64::from).collect();
        let data = one_variable(pts, &vals);
        let lags = CentredLags::new(5, 1.0, 0.25).unwrap();
        let out = gamv(&data, &request(lags, vec![direction(90.0, 10.0)])).unwrap();
        for k in 1..=5_u32 {
            let lag = out[0].lags[usize::try_from(k).unwrap() + 1];
            let kf = f64::from(k);
            assert_eq!(lag.pairs, usize::try_from(11 - k).unwrap(), "pairs at {k}");
            assert_abs_diff_eq!(lag.gamma.unwrap(), kf * kf / 2.0, epsilon = 1e-12);
            assert_abs_diff_eq!(lag.mean_distance.unwrap(), kf, epsilon = 1e-12);
        }
    }

    /// On a 5 × 5 grid with value = x, pairs straight north have equal
    /// values. With atol = 10° no diagonal pair qualifies (the nearest is
    /// atan(1/4) ≈ 14°), so γ is exactly 0 to the north.
    #[test]
    fn constant_across_the_direction_gives_zero() {
        let mut pts = Vec::new();
        let mut vals = Vec::new();
        for y in 0..5 {
            for x in 0..5 {
                pts.push(p(f64::from(x), f64::from(y)));
                vals.push(f64::from(x));
            }
        }
        let data = one_variable(pts, &vals);
        let lags = CentredLags::new(4, 1.0, 0.25).unwrap();
        let out = gamv(&data, &request(lags, vec![direction(0.0, 10.0)])).unwrap();
        for (k, lag) in (1..=4_usize).zip(&out[0].lags[2..]) {
            assert_eq!(lag.pairs, 5 * (5 - k), "five columns of 5 − k pairs at {k}");
            assert_eq!(
                lag.gamma.map(f64::to_bits),
                Some(0.0_f64.to_bits()),
                "lag {k}"
            );
        }
    }

    #[test]
    fn tail_and_head_follow_the_pair_order() {
        // Datum 0 at the origin, datum 1 one unit north. Tail variable a,
        // head variable b.
        let data = Data::new(
            vec![p(0.0, 0.0), p(0.0, 1.0)],
            vec![vec![Some(1.0), Some(2.0)], vec![Some(10.0), Some(20.0)]],
        )
        .unwrap();
        let spec = VariogramSpec {
            tail: 0,
            head: 1,
            kind: VariogramKind::Semivariogram,
            standardize: false,
        };
        let lags = CentredLags::new(2, 1.0, 0.5).unwrap();
        let mut r = request(lags, vec![direction(0.0, 10.0), direction(180.0, 10.0)]);
        r.variograms = vec![spec];
        let out = gamv(&data, &r).unwrap();
        // North: the pair is as given; tail a at datum 0, head b at datum 1.
        let north = out[0].lags[2];
        assert_eq!(north.tail_mean, Some(1.0));
        assert_eq!(north.head_mean, Some(20.0));
        assert_eq!(north.gamma, Some(0.5 * 19.0 * 19.0));
        // South: swapped; tail a at datum 1, head b at datum 0.
        let south = out[1].lags[2];
        assert_eq!(south.tail_mean, Some(2.0));
        assert_eq!(south.head_mean, Some(10.0));
        assert_eq!(south.gamma, Some(0.5 * 8.0 * 8.0));
    }

    #[test]
    fn coincident_data_fill_lag_zero_without_self_pairs() {
        // Two data share a location; the third is apart. D-017: lag 0 holds
        // the one coincident pair, never a datum paired with itself.
        let pts = vec![p(0.0, 0.0), p(0.0, 0.0), p(0.0, 2.0)];
        let data = one_variable(pts, &[1.0, 4.0, 9.0]);
        let lags = CentredLags::new(3, 1.0, 0.5).unwrap();
        let out = gamv(
            &data,
            &request(lags, vec![direction(0.0, 45.0), direction(0.0, 90.0)]),
        )
        .unwrap();
        assert_eq!((out[0].lags[0].pairs, out[0].lags[0].orderings), (1, 1));
        assert_eq!(
            (out[1].lags[0].pairs, out[1].lags[0].orderings),
            (1, 2),
            "omni"
        );
        assert_eq!(out[0].lags[0].gamma, Some(4.5));
    }

    #[test]
    fn gslib_policies_select_every_legacy_option() {
        assert_eq!(
            GamvPolicies::gslib(),
            GamvPolicies {
                pi: PiConstant::Gslib,
                omni_relative_log: OmniRelativeLog::GslibDuplicate,
            }
        );
        assert_eq!(
            GamvPolicies::default(),
            GamvPolicies {
                pi: PiConstant::Exact,
                omni_relative_log: OmniRelativeLog::Swapped,
            }
        );
    }

    // --- Step 3: every measure, hand-computed ----------------------------
    //
    // Three data one unit apart going north; values a = [1, 2, 4] and
    // b = [3, 1, 5]. With xlag = 1, xltol = 0.25, lag index 2 (distance 1)
    // holds the pairs (0, 1) and (1, 2), both "as given" for a northward
    // direction, so the earlier datum is the tail.

    const A: usize = 0;
    const B: usize = 1;

    fn three_north(extra: &[[f64; 3]]) -> Data {
        let mut columns = vec![
            vec![Some(1.0), Some(2.0), Some(4.0)],
            vec![Some(3.0), Some(1.0), Some(5.0)],
        ];
        columns.extend(extra.iter().map(|c| c.iter().copied().map(Some).collect()));
        Data::new(vec![p(0.0, 0.0), p(0.0, 1.0), p(0.0, 2.0)], columns).unwrap()
    }

    fn spec(tail: usize, head: usize, kind: VariogramKind) -> VariogramSpec {
        VariogramSpec {
            tail,
            head,
            kind,
            standardize: false,
        }
    }

    /// Lag index 2 of each variogram, for a north direction.
    fn lag_one(
        data: &Data,
        variograms: Vec<VariogramSpec>,
        sills: Vec<Option<f64>>,
        atol: f64,
        policies: GamvPolicies,
    ) -> Result<Vec<GamvLag>> {
        let r = GamvRequest {
            lags: CentredLags::new(2, 1.0, 0.25)?,
            directions: vec![direction(0.0, atol)],
            variograms,
            sills,
            policies,
        };
        Ok(gamv(data, &r)?.iter().map(|v| v.lags[2]).collect())
    }

    fn gamma_of(data: &Data, kind: VariogramKind, sills: Vec<Option<f64>>) -> Option<f64> {
        lag_one(
            data,
            vec![spec(A, A, kind)],
            sills,
            10.0,
            GamvPolicies::default(),
        )
        .unwrap()[0]
            .gamma
    }

    #[test]
    fn every_measure_matches_hand_computation() {
        use VariogramKind as K;
        let data = three_north(&[]);
        // Pairs (t, h) of a: (1, 2) and (2, 4).
        let cases = [
            (K::Semivariogram, (1.0 + 4.0) / 4.0),
            (K::Covariance, 5.0 - 4.5),
            (K::Correlogram, 1.0),
            (K::GeneralRelative, 2.5 / (2.25 * 2.25)),
            (K::PairwiseRelative, (4.0 / 9.0 + 4.0 / 9.0) / 4.0),
            (K::Logarithmic, 2.0_f64.ln().powi(2) / 2.0),
            (K::Madogram, (1.0 + 2.0) / 4.0),
        ];
        for (kind, expected) in cases {
            let got = gamma_of(&data, kind, vec![]).unwrap();
            assert_abs_diff_eq!(got, expected, epsilon = 1e-12);
        }
        let lag = lag_one(
            &data,
            vec![spec(A, A, K::Correlogram)],
            vec![],
            10.0,
            GamvPolicies::default(),
        )
        .unwrap()[0];
        assert_abs_diff_eq!(lag.tail_variance.unwrap(), 0.25, epsilon = 1e-12);
        assert_abs_diff_eq!(lag.head_variance.unwrap(), 1.0, epsilon = 1e-12);
        let semi = lag_one(
            &data,
            vec![spec(A, A, K::Semivariogram)],
            vec![],
            10.0,
            GamvPolicies::default(),
        )
        .unwrap()[0];
        assert_eq!(
            (semi.tail_variance, semi.head_variance),
            (None, None),
            "correlogram only"
        );
    }

    #[test]
    fn covariance_as_variogram_is_sill_minus_covariance() {
        let data = three_north(&[]);
        let sills = vec![Some(2.0), Some(1.0)];
        let g = gamma_of(&data, VariogramKind::CovarianceAsVariogram, sills.clone()).unwrap();
        let c = gamma_of(&data, VariogramKind::Covariance, sills).unwrap();
        let root = 2.0_f64.sqrt();
        assert_abs_diff_eq!(g, root.mul_add(root, -c), epsilon = 1e-15);
        assert_eq!(
            gamma_of(&data, VariogramKind::CovarianceAsVariogram, vec![]),
            None,
            "no sill: undefined (GSLIB writes −999)"
        );
        assert_eq!(
            gamma_of(
                &data,
                VariogramKind::CovarianceAsVariogram,
                vec![Some(-999.0), None]
            ),
            None,
            "GSLIB's −999 marker is not a sill"
        );
    }

    #[test]
    fn cross_semivariogram_matches_hand_computation() {
        // γ_ab = Σ (a_j − a_i)(b_j − b_i) / 2N = ((2−1)(1−3) + (4−2)(5−1)) / 4.
        let data = three_north(&[]);
        let lag = lag_one(
            &data,
            vec![spec(A, B, VariogramKind::CrossSemivariogram)],
            vec![],
            10.0,
            GamvPolicies::default(),
        )
        .unwrap()[0];
        assert_abs_diff_eq!(lag.gamma.unwrap(), 1.5, epsilon = 1e-12);
        assert_abs_diff_eq!(lag.tail_mean.unwrap(), 2.25, epsilon = 1e-12);
        assert_abs_diff_eq!(lag.head_mean.unwrap(), 2.5, epsilon = 1e-12);
        // Not doubled when omnidirectional.
        let omni = lag_one(
            &data,
            vec![spec(A, B, VariogramKind::CrossSemivariogram)],
            vec![],
            90.0,
            GamvPolicies::default(),
        )
        .unwrap()[0];
        assert_eq!(omni.pairs, 2);
    }

    #[test]
    fn cross_semivariogram_of_a_variable_with_itself_is_its_semivariogram() {
        let data = three_north(&[]);
        for atol in [10.0, 90.0] {
            let lags = lag_one(
                &data,
                vec![
                    spec(B, B, VariogramKind::CrossSemivariogram),
                    spec(B, B, VariogramKind::Semivariogram),
                ],
                vec![],
                atol,
                GamvPolicies::default(),
            )
            .unwrap();
            assert_abs_diff_eq!(
                lags[0].gamma.unwrap(),
                lags[1].gamma.unwrap(),
                epsilon = 1e-12
            );
        }
    }

    #[test]
    fn standardisation_divides_by_the_sill() {
        let data = three_north(&[]);
        let mut s = spec(A, A, VariogramKind::Semivariogram);
        s.standardize = true;
        let run = |v: VariogramSpec, sills: Vec<Option<f64>>| {
            lag_one(&data, vec![v], sills, 10.0, GamvPolicies::default())
        };
        let lag = run(s, vec![Some(2.0), None]).unwrap()[0];
        assert_abs_diff_eq!(lag.gamma.unwrap(), 1.25 / 2.0, epsilon = 1e-12);
        // D-023: GSLIB silently skips these; here they are errors.
        assert!(run(s, vec![Some(0.0), None]).is_err(), "zero sill");
        assert!(run(s, vec![None, None]).is_err(), "missing sill");
        assert!(run(s, vec![]).is_err(), "no sills");
        let mut cov = spec(A, A, VariogramKind::Covariance);
        cov.standardize = true;
        assert!(
            run(cov, vec![Some(2.0), None]).is_err(),
            "not a semivariogram"
        );
        let mut cross = spec(A, B, VariogramKind::Semivariogram);
        cross.standardize = true;
        assert!(
            run(cross, vec![Some(2.0), Some(2.0)]).is_err(),
            "tail ≠ head"
        );
        assert!(run(s, vec![Some(2.0)]).is_err(), "wrong number of sills");
        assert!(run(s, vec![Some(f64::NAN), None]).is_err(), "NaN sill");
    }

    #[test]
    fn undefined_measures_are_none_with_pairs_reported() {
        // c is constant, so the lag's tail and head have no spread; z is all
        // zero, so the general relative mean vanishes.
        let data = three_north(&[[1.0, 1.0, 1.0], [0.0, 0.0, 0.0]]);
        let lags = lag_one(
            &data,
            vec![
                spec(2, 2, VariogramKind::Correlogram),
                spec(3, 3, VariogramKind::GeneralRelative),
            ],
            vec![],
            10.0,
            GamvPolicies::default(),
        )
        .unwrap();
        for lag in &lags {
            assert_eq!(lag.gamma, None, "{lag:?}");
            assert_eq!(lag.pairs, 2, "pairs still reported");
        }
        assert_eq!(lags[0].tail_variance, Some(0.0));
    }

    #[test]
    fn relative_and_log_skip_pairs_outside_their_domain() {
        // d sums to zero over each unit pair; e has a zero in the middle.
        let data = three_north(&[[1.0, -1.0, 1.0], [1.0, 0.0, 2.0]]);
        let lags = lag_one(
            &data,
            vec![
                spec(2, 2, VariogramKind::PairwiseRelative),
                spec(3, 3, VariogramKind::Logarithmic),
            ],
            vec![],
            10.0,
            GamvPolicies::default(),
        )
        .unwrap();
        for lag in &lags {
            assert_eq!((lag.pairs, lag.skipped, lag.gamma), (0, 2, None), "{lag:?}");
        }
    }

    /// D-021: tail ≠ head, omnidirectional. Unit pairs as (t, h | t', h'):
    /// (1, 1 | 2, 3) and (2, 5 | 4, 1).
    #[test]
    fn omnidirectional_pairwise_relative_uses_swapped_values() {
        let data = three_north(&[]);
        let pr = |t: f64, h: f64| (2.0 * (h - t) / (h + t)).powi(2);
        let run = |policy| {
            let policies = GamvPolicies {
                omni_relative_log: policy,
                ..GamvPolicies::default()
            };
            lag_one(
                &data,
                vec![spec(A, B, VariogramKind::PairwiseRelative)],
                vec![],
                90.0,
                policies,
            )
            .unwrap()[0]
                .gamma
                .unwrap()
        };
        let swapped = (pr(1.0, 1.0) + pr(2.0, 3.0) + pr(2.0, 5.0) + pr(4.0, 1.0)) / 8.0;
        let duplicate = (pr(1.0, 1.0) + pr(2.0, 5.0)) / 4.0;
        assert_abs_diff_eq!(run(OmniRelativeLog::Swapped), swapped, epsilon = 1e-12);
        assert_abs_diff_eq!(
            run(OmniRelativeLog::GslibDuplicate),
            duplicate,
            epsilon = 1e-12
        );
    }

    #[test]
    fn omnidirectional_policies_agree_when_tail_equals_head() {
        let data = three_north(&[]);
        for kind in [VariogramKind::PairwiseRelative, VariogramKind::Logarithmic] {
            let run = |policy| {
                let policies = GamvPolicies {
                    omni_relative_log: policy,
                    ..GamvPolicies::default()
                };
                lag_one(&data, vec![spec(B, B, kind)], vec![], 90.0, policies).unwrap()[0]
            };
            assert_eq!(
                run(OmniRelativeLog::Swapped),
                run(OmniRelativeLog::GslibDuplicate),
                "{kind:?}"
            );
        }
    }

    /// D-021: under the legacy policy the swapped log contribution is
    /// guarded by the swapped values (both 1) but computed from the
    /// unswapped head value 0: ln 0 = −∞. That is an error, never an
    /// infinite γ. The default computes it from the swapped values.
    #[test]
    fn legacy_log_defect_is_an_error_not_infinity() {
        let data = Data::new(
            vec![p(0.0, 0.0), p(0.0, 1.0)],
            vec![vec![Some(1.0), Some(1.0)], vec![Some(1.0), Some(0.0)]],
        )
        .unwrap();
        let run = |policy| {
            let policies = GamvPolicies {
                omni_relative_log: policy,
                ..GamvPolicies::default()
            };
            lag_one(
                &data,
                vec![spec(0, 1, VariogramKind::Logarithmic)],
                vec![],
                90.0,
                policies,
            )
        };
        assert!(run(OmniRelativeLog::GslibDuplicate).is_err());
        let lag = run(OmniRelativeLog::Swapped).unwrap()[0];
        assert_eq!((lag.pairs, lag.skipped), (1, 1));
        assert_eq!(lag.gamma, Some(0.0));
    }

    #[test]
    fn gslib_type_numbers() {
        use VariogramKind as K;
        let all = [
            K::Semivariogram,
            K::CrossSemivariogram,
            K::Covariance,
            K::CovarianceAsVariogram,
            K::Correlogram,
            K::GeneralRelative,
            K::PairwiseRelative,
            K::Logarithmic,
            K::Madogram,
        ];
        let types: Vec<i32> = all.iter().map(|k| k.gslib_type()).collect();
        assert_eq!(types, vec![1, 2, 3, -3, 4, 5, 6, 7, 8]);
    }

    /// D-043: an omnidirectional cross measure averages both orderings, so
    /// reversing the order of the data rows cannot change it, whereas a
    /// single ordering would.
    #[test]
    fn omnidirectional_cross_measures_ignore_row_order() {
        let pts: Vec<Point3> = (0..12_u32)
            .map(|i| p(f64::from(i % 4), f64::from(i.div_euclid(4))))
            .collect();
        let a: Vec<Option<f64>> = (0..12).map(|i| Some(f64::from(i * i % 7))).collect();
        let b: Vec<Option<f64>> = (0..12).map(|i| Some(f64::from(3 * i % 5) + 0.5)).collect();
        let run = |at: Vec<Point3>, tail: Vec<Option<f64>>, head: Vec<Option<f64>>| {
            let data = Data::new(at, vec![tail, head]).unwrap();
            let mut r = request(
                CentredLags::new(3, 1.0, 0.5).unwrap(),
                vec![direction(0.0, 90.0)],
            );
            r.variograms = [VariogramKind::Covariance, VariogramKind::Madogram]
                .map(|kind| VariogramSpec {
                    tail: 0,
                    head: 1,
                    kind,
                    standardize: false,
                })
                .to_vec();
            gamv(&data, &r).unwrap()
        };
        let forward = run(pts.clone(), a.clone(), b.clone());
        let rev = |v: Vec<Option<f64>>| v.into_iter().rev().collect::<Vec<_>>();
        let backward = run(pts.into_iter().rev().collect(), rev(a), rev(b));
        for (f, r) in forward.iter().zip(&backward) {
            for (lf, lr) in f.lags.iter().zip(&r.lags) {
                assert_eq!(lf.pairs, lr.pairs);
                match (lf.gamma, lr.gamma) {
                    (Some(x), Some(y)) => assert_abs_diff_eq!(x, y, epsilon = 1e-12),
                    (x, y) => assert_eq!(x, y),
                }
            }
        }
    }

    /// One omnidirectional semivariogram equals PR #4's
    /// `experimental_variogram` (same pairs, same enumeration), up to the
    /// rounding of its fused accumulation.
    #[test]
    fn omnidirectional_semivariogram_matches_experimental_variogram() {
        use crate::sample::Sample;
        use crate::variogram::experimental_variogram;
        let pts: Vec<Point3> = (0..40)
            .map(|i| p(f64::from(i * 7 % 13), f64::from(i * 5 % 11)))
            .collect();
        let vals: Vec<f64> = (0..40).map(|i| f64::from(i * i % 17) / 3.0).collect();
        let lags = CentredLags::new(6, 2.0, 1.5).unwrap();
        let samples: Vec<Sample> = pts
            .iter()
            .zip(&vals)
            .map(|(&q, &v)| Sample::new(q, v).unwrap())
            .collect();
        let reference = experimental_variogram(&samples, &lags);
        let out = gamv(
            &one_variable(pts, &vals),
            &request(lags, vec![direction(0.0, 90.0)]),
        )
        .unwrap();
        for (ours, theirs) in out[0].lags.iter().zip(&reference) {
            assert_eq!(ours.pairs, theirs.pairs);
            match (ours.gamma, theirs.gamma) {
                (Some(x), Some(y)) => assert_abs_diff_eq!(x, y, epsilon = 1e-12 * y.abs().max(1.0)),
                (x, y) => assert_eq!(x, y),
            }
            if let (Some(m), Some(pm)) = (ours.tail_mean, theirs.pair_mean) {
                assert_abs_diff_eq!(m, pm, epsilon = 1e-12 * pm.abs().max(1.0));
            }
        }
    }

    // --- Properties ------------------------------------------------------

    fn scattered() -> impl Strategy<Value = (Vec<Point3>, Vec<f64>)> {
        prop::collection::vec((0.0..20.0_f64, 0.0..20.0_f64, -5.0..5.0_f64), 2..25).prop_map(
            |rows| {
                let pts = rows.iter().map(|&(x, y, _)| p(x, y)).collect();
                let vals = rows.iter().map(|&(_, _, v)| v).collect();
                (pts, vals)
            },
        )
    }

    proptest! {
        /// D-017: omnidirectional pair counts are exactly twice the number
        /// of unordered pairs in each lag, counted independently here.
        #[test]
        fn omnidirectional_pairs_counted_once((pts, vals) in scattered()) {
            let lags = CentredLags::new(8, 2.0, 1.5).unwrap();
            let data = one_variable(pts.clone(), &vals);
            let out = gamv(&data, &request(lags, vec![direction(0.0, 90.0)])).unwrap();
            let reach = (8.0 + 0.5) * 2.0;
            for (index, lag) in out[0].lags.iter().enumerate().skip(1) {
                let centre = 2.0 * to_f64(index - 1);
                let mut expected = 0;
                for (i, a) in pts.iter().enumerate() {
                    for b in &pts[i + 1..] {
                        let h = a.distance(*b);
                        if h <= reach && h > 0.0 && (h - centre).abs() <= 1.5 {
                            expected += 1;
                        }
                    }
                }
                prop_assert_eq!(lag.pairs, expected, "lag {}", index);
                prop_assert_eq!(lag.orderings, 2 * expected, "lag {}", index);
            }
        }

        /// γ ≥ 0, and each lag's mean distance lies within its class.
        #[test]
        fn gamma_non_negative_and_mean_distance_in_class(
            (pts, vals) in scattered(), azm in 0.0..180.0_f64
        ) {
            let lags = CentredLags::new(6, 2.5, 1.25).unwrap();
            let data = one_variable(pts, &vals);
            let out = gamv(&data, &request(lags, vec![direction(azm, 30.0)])).unwrap();
            for (index, lag) in out[0].lags.iter().enumerate().skip(1) {
                if let (Some(g), Some(m)) = (lag.gamma, lag.mean_distance) {
                    prop_assert!(g >= 0.0);
                    let centre = 2.5 * to_f64(index - 1);
                    prop_assert!((m - centre).abs() <= 1.25 + 1e-12, "lag {}: {}", index, m);
                }
            }
        }

        /// The correlogram lies in [−1, 1] (up to rounding).
        #[test]
        fn correlogram_is_bounded((pts, vals) in scattered(), atol in 10.0..90.0_f64) {
            let lags = CentredLags::new(6, 2.5, 1.25).unwrap();
            let data = one_variable(pts, &vals);
            let mut r = request(lags, vec![direction(0.0, atol)]);
            r.variograms = vec![VariogramSpec {
                tail: 0,
                head: 0,
                kind: VariogramKind::Correlogram,
                standardize: false,
            }];
            for lag in &gamv(&data, &r).unwrap()[0].lags {
                if let Some(g) = lag.gamma {
                    prop_assert!(g.abs() <= 1.0 + 1e-9, "{}", g);
                }
            }
        }

        /// A missing value removes exactly the pairs of its datum: the
        /// result equals the one computed without that datum.
        #[test]
        fn missing_values_never_pair((pts, vals) in scattered(), drop in 0..25_usize) {
            let drop = drop % pts.len();
            let lags = CentredLags::new(6, 2.5, 1.25).unwrap();
            let dirs = vec![direction(30.0, 40.0), direction(0.0, 90.0)];
            let mut column: Vec<Option<f64>> = vals.iter().copied().map(Some).collect();
            column[drop] = None;
            let with_none = Data::new(pts.clone(), vec![column]).unwrap();
            let mut kept_pts = pts;
            let mut kept_vals = vals;
            kept_pts.remove(drop);
            kept_vals.remove(drop);
            let removed = one_variable(kept_pts, &kept_vals);
            let a = gamv(&with_none, &request(lags, dirs.clone())).unwrap();
            let b = gamv(&removed, &request(lags, dirs)).unwrap();
            for (x, y) in a.iter().zip(&b) {
                for (lx, ly) in x.lags.iter().zip(&y.lags) {
                    prop_assert_eq!(lx.pairs, ly.pairs);
                    match (lx.gamma, ly.gamma) {
                        (Some(gx), Some(gy)) => prop_assert!((gx - gy).abs() <= 1e-12 * gy.abs().max(1.0)),
                        (gx, gy) => prop_assert_eq!(gx, gy),
                    }
                }
            }
        }
    }
}
