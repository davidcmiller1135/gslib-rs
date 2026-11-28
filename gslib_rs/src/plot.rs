//! Plots rendered with `plotters` to PNG or SVG.
//!
//! Each plot is a value implementing [`Figure`]; [`save`] picks the backend
//! from the file extension, so drawing code is written once.

use std::ops::Range;
use std::path::Path;

use plotters::coord::Shift;
use plotters::prelude::*;

use crate::error::{Error, Result};
use crate::num::{floor_index, to_f64};

/// Result of drawing on a backend.
pub type DrawResult<DB> =
    std::result::Result<(), DrawingAreaErrorKind<<DB as DrawingBackend>::ErrorType>>;

/// Something that can be drawn onto any plotters backend.
pub trait Figure {
    /// Image size in pixels.
    fn size(&self) -> (u32, u32);

    /// Draw onto `root`.
    ///
    /// # Errors
    /// Any backend error.
    fn draw<DB: DrawingBackend>(&self, root: &DrawingArea<DB, Shift>) -> DrawResult<DB>;
}

/// Render `figure` to `path`: SVG for a `.svg` extension, PNG otherwise.
///
/// # Errors
/// [`Error::Plot`] if rendering or writing fails.
pub fn save(figure: &impl Figure, path: &Path) -> Result<()> {
    let is_svg = path
        .extension()
        .and_then(|e| e.to_str())
        .is_some_and(|e| e.eq_ignore_ascii_case("svg"));
    let plot_err = |e: &dyn std::fmt::Display| Error::Plot(format!("{}: {e}", path.display()));
    if is_svg {
        let root = SVGBackend::new(path, figure.size()).into_drawing_area();
        figure.draw(&root).map_err(|e| plot_err(&e))?;
        root.present().map_err(|e| plot_err(&e))
    } else {
        let root = BitMapBackend::new(path, figure.size()).into_drawing_area();
        figure.draw(&root).map_err(|e| plot_err(&e))?;
        root.present().map_err(|e| plot_err(&e))
    }
}

/// `[min, max]` of the finite values, widened if degenerate; `None` if there
/// are no finite values.
fn extent(values: impl IntoIterator<Item = f64>) -> Option<Range<f64>> {
    let (lo, hi) = values
        .into_iter()
        .filter(|v| v.is_finite())
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), v| {
            (lo.min(v), hi.max(v))
        });
    if lo > hi {
        None
    } else if hi - lo < f64::EPSILON * lo.abs().max(1.0) {
        Some((lo - 1.0)..(hi + 1.0))
    } else {
        Some(lo..hi)
    }
}

/// Upper end of `[0, max]` for non-negative quantities.
fn zero_based(values: impl IntoIterator<Item = f64>) -> Range<f64> {
    let top = values
        .into_iter()
        .filter(|v| v.is_finite())
        .fold(0.0, f64::max);
    0.0..(if top > 0.0 { top * 1.05 } else { 1.0 })
}

/// One labelled series of experimental variogram points.
#[derive(Debug, Clone, PartialEq)]
pub struct VariogramSeries {
    /// Legend label (e.g. the direction).
    pub label: String,
    /// `(lag distance, γ)` points; lags without a value should be omitted.
    pub points: Vec<(f64, f64)>,
}

/// Experimental variograms (one series per direction or variogram) with an
/// optional model curve.
#[derive(Debug, Clone, PartialEq)]
pub struct VariogramPlot {
    /// The experimental series.
    pub experimental: Vec<VariogramSeries>,
    /// Optional `(distance, γ)` model curve.
    pub model: Option<Vec<(f64, f64)>>,
}

/// Series colours, cycled.
const SERIES_COLOURS: [RGBColor; 6] = [RED, GREEN, MAGENTA, CYAN, BLACK, YELLOW];

impl Figure for VariogramPlot {
    fn size(&self) -> (u32, u32) {
        (800, 600)
    }

    fn draw<DB: DrawingBackend>(&self, root: &DrawingArea<DB, Shift>) -> DrawResult<DB> {
        let model = self.model.as_deref().unwrap_or_default();
        let all = || {
            self.experimental
                .iter()
                .flat_map(|s| &s.points)
                .chain(model)
        };
        root.fill(&WHITE)?;
        let mut chart = ChartBuilder::on(root)
            .caption("Variogram", ("sans-serif", 24))
            .margin(20)
            .x_label_area_size(40)
            .y_label_area_size(50)
            .build_cartesian_2d(
                zero_based(all().map(|p| p.0)),
                zero_based(all().map(|p| p.1)),
            )?;
        chart
            .configure_mesh()
            .x_desc("distance")
            .y_desc("γ")
            .draw()?;
        for (series, colour) in self.experimental.iter().zip(SERIES_COLOURS.iter().cycle()) {
            let colour = *colour;
            chart
                .draw_series(
                    series
                        .points
                        .iter()
                        .map(move |&p| Circle::new(p, 4, colour.filled())),
                )?
                .label(series.label.clone())
                .legend(move |(x, y)| Circle::new((x + 10, y), 4, colour.filled()));
        }
        if !model.is_empty() {
            chart
                .draw_series(LineSeries::new(model.iter().copied(), BLUE.stroke_width(2)))?
                .label("model")
                .legend(|(x, y)| PathElement::new([(x, y), (x + 20, y)], BLUE.stroke_width(2)));
        }
        chart.configure_series_labels().border_style(BLACK).draw()
    }
}

/// A scatter plot, optionally with the least-squares line.
#[derive(Debug, Clone, PartialEq)]
pub struct ScatterPlot {
    /// `(x, y)` points.
    pub points: Vec<(f64, f64)>,
    /// Draw the ordinary least-squares regression of y on x.
    pub regression: bool,
}

impl ScatterPlot {
    /// Slope and intercept of the least-squares line, if x has variance.
    #[must_use]
    pub fn regression_line(&self) -> Option<(f64, f64)> {
        let n = to_f64(self.points.len());
        let mx = self.points.iter().map(|p| p.0).sum::<f64>() / n;
        let my = self.points.iter().map(|p| p.1).sum::<f64>() / n;
        let (sxy, sxx) = self.points.iter().fold((0.0, 0.0), |(sxy, sxx), &(x, y)| {
            ((x - mx).mul_add(y - my, sxy), (x - mx).mul_add(x - mx, sxx))
        });
        (sxx > 0.0).then(|| {
            let slope = sxy / sxx;
            (slope, slope.mul_add(-mx, my))
        })
    }
}

impl Figure for ScatterPlot {
    fn size(&self) -> (u32, u32) {
        (600, 600)
    }

    fn draw<DB: DrawingBackend>(&self, root: &DrawingArea<DB, Shift>) -> DrawResult<DB> {
        let xr = extent(self.points.iter().map(|p| p.0)).unwrap_or(0.0..1.0);
        let yr = extent(self.points.iter().map(|p| p.1)).unwrap_or(0.0..1.0);
        root.fill(&WHITE)?;
        let mut chart = ChartBuilder::on(root)
            .caption("Scatter plot", ("sans-serif", 24))
            .margin(20)
            .x_label_area_size(40)
            .y_label_area_size(50)
            .build_cartesian_2d(xr.clone(), yr)?;
        chart.configure_mesh().draw()?;
        chart
            .draw_series(self.points.iter().map(|&p| Circle::new(p, 3, RED.filled())))?
            .label("data")
            .legend(|(x, y)| Circle::new((x + 10, y), 3, RED.filled()));
        if let Some((slope, intercept)) = self.regression.then(|| self.regression_line()).flatten()
        {
            let line = [xr.start, xr.end].map(|x| (x, slope.mul_add(x, intercept)));
            chart
                .draw_series(LineSeries::new(line, BLUE.stroke_width(2)))?
                .label(format!("y = {slope:.3}x + {intercept:.3}"))
                .legend(|(x, y)| PathElement::new([(x, y), (x + 20, y)], BLUE.stroke_width(2)));
        }
        chart.configure_series_labels().border_style(BLACK).draw()
    }
}

/// A histogram of values in equal-width classes.
#[derive(Debug, Clone, PartialEq)]
pub struct HistogramPlot {
    /// Class range; values outside fall in the end classes.
    pub range: Range<f64>,
    /// Count in each class.
    pub counts: Vec<usize>,
}

impl HistogramPlot {
    /// Bin `values` into `bins` classes over `range` (or the data's range).
    ///
    /// # Errors
    /// [`Error::InvalidInput`] if `bins` is zero, there are no finite values,
    /// or the range is not finite and increasing.
    pub fn new(values: &[f64], bins: usize, range: Option<Range<f64>>) -> Result<Self> {
        if bins == 0 {
            return Err(Error::invalid("histogram", "needs at least one bin"));
        }
        let range = match range {
            Some(r) => r,
            None => extent(values.iter().copied())
                .ok_or_else(|| Error::invalid("histogram", "no finite values"))?,
        };
        if !(range.start.is_finite() && range.end.is_finite() && range.start < range.end) {
            return Err(Error::invalid(
                "histogram range",
                "must be finite and increasing",
            ));
        }
        let width = (range.end - range.start) / to_f64(bins);
        let mut counts = vec![0; bins];
        for v in values.iter().filter(|v| v.is_finite()) {
            let class = floor_index((v - range.start) / width).map_or(0, |c| c.min(bins - 1));
            if let Some(c) = counts.get_mut(class) {
                *c += 1;
            }
        }
        Ok(Self { range, counts })
    }
}

impl Figure for HistogramPlot {
    fn size(&self) -> (u32, u32) {
        (600, 400)
    }

    fn draw<DB: DrawingBackend>(&self, root: &DrawingArea<DB, Shift>) -> DrawResult<DB> {
        let width = (self.range.end - self.range.start) / to_f64(self.counts.len().max(1));
        root.fill(&WHITE)?;
        let mut chart = ChartBuilder::on(root)
            .caption("Histogram", ("sans-serif", 24))
            .margin(20)
            .x_label_area_size(40)
            .y_label_area_size(50)
            .build_cartesian_2d(
                self.range.clone(),
                zero_based(self.counts.iter().map(|&c| to_f64(c))),
            )?;
        chart.configure_mesh().draw()?;
        chart.draw_series(self.counts.iter().enumerate().map(|(i, &c)| {
            let x0 = width.mul_add(to_f64(i), self.range.start);
            Rectangle::new([(x0, 0.0), (x0 + width, to_f64(c))], BLUE.filled())
        }))?;
        Ok(())
    }
}

/// A 2-D grid of values drawn as coloured cells, x fastest, with y increasing
/// upward (GSLIB `pixelplt` orientation).
#[derive(Debug, Clone, PartialEq)]
pub struct PixelMap {
    values: Vec<f64>,
    nx: usize,
    range: Range<f64>,
}

impl PixelMap {
    /// # Errors
    /// [`Error::InvalidInput`] if `nx` is zero, `values.len()` is not a
    /// multiple of `nx`, or there are no finite values.
    pub fn new(values: Vec<f64>, nx: usize) -> Result<Self> {
        if nx == 0 || values.is_empty() || !values.len().is_multiple_of(nx) {
            return Err(Error::invalid(
                "pixel map",
                "values must fill whole rows of nx > 0 cells",
            ));
        }
        let range = extent(values.iter().copied())
            .ok_or_else(|| Error::invalid("pixel map", "no finite values"))?;
        Ok(Self { values, nx, range })
    }

    /// Blue → green → red ramp for `v` within the map's range.
    fn colour(&self, v: f64) -> RGBColor {
        let t = ((v - self.range.start) / (self.range.end - self.range.start)).clamp(0.0, 1.0);
        let channel =
            |x: f64| u8::try_from(floor_index(x * 255.0).unwrap_or(0).min(255)).unwrap_or(u8::MAX);
        RGBColor(
            channel(t),
            channel(1.0 - 2.0_f64.mul_add(t, -1.0).abs()),
            channel(1.0 - t),
        )
    }
}

impl Figure for PixelMap {
    fn size(&self) -> (u32, u32) {
        (800, 800)
    }

    fn draw<DB: DrawingBackend>(&self, root: &DrawingArea<DB, Shift>) -> DrawResult<DB> {
        let ny = self.values.len().checked_div(self.nx).unwrap_or(0);
        root.fill(&WHITE)?;
        let mut chart = ChartBuilder::on(root)
            .margin(20)
            .x_label_area_size(30)
            .y_label_area_size(30)
            .build_cartesian_2d(0.0..to_f64(self.nx), 0.0..to_f64(ny))?;
        chart.configure_mesh().disable_mesh().draw()?;
        chart.draw_series(self.values.chunks_exact(self.nx).enumerate().flat_map(
            |(iy, row)| {
                row.iter().enumerate().map(move |(ix, &v)| {
                    let (x, y) = (to_f64(ix), to_f64(iy));
                    Rectangle::new([(x, y), (x + 1.0, y + 1.0)], self.colour(v).filled())
                })
            },
        ))?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn histogram_counts_and_edges() {
        let h = HistogramPlot::new(&[0.0, 0.4, 0.6, 1.0, 9.0, -3.0], 2, Some(0.0..1.0)).unwrap();
        assert_eq!(h.counts, vec![3, 3]);
        assert!(HistogramPlot::new(&[1.0], 0, None).is_err(), "zero bins");
        assert!(
            HistogramPlot::new(&[f64::NAN], 3, None).is_err(),
            "no finite values"
        );
    }

    #[test]
    fn regression_line_fits_exact_line() {
        let s = ScatterPlot {
            points: vec![(0.0, 1.0), (1.0, 3.0), (2.0, 5.0)],
            regression: true,
        };
        let (m, b) = s.regression_line().unwrap();
        assert_abs_diff_eq!(m, 2.0, epsilon = 1e-12);
        assert_abs_diff_eq!(b, 1.0, epsilon = 1e-12);
        let vertical = ScatterPlot {
            points: vec![(1.0, 0.0), (1.0, 2.0)],
            regression: true,
        };
        assert_eq!(vertical.regression_line(), None);
    }

    #[test]
    fn extent_widens_degenerate_ranges() {
        assert_eq!(extent([2.0, 2.0]), Some(1.0..3.0));
        assert_eq!(extent([f64::NAN]), None);
    }

    #[test]
    fn renders_png_and_svg() {
        let dir = tempfile::tempdir().unwrap();
        let fig = VariogramPlot {
            experimental: vec![
                VariogramSeries {
                    label: "direction 1".into(),
                    points: vec![(1.0, 0.5), (2.0, 0.9)],
                },
                VariogramSeries {
                    label: "direction 2".into(),
                    points: vec![(1.0, 0.7)],
                },
            ],
            model: Some(vec![(0.0, 0.0), (2.0, 1.0)]),
        };
        for name in ["v.png", "v.svg"] {
            let path = dir.path().join(name);
            save(&fig, &path).unwrap();
            assert!(
                std::fs::metadata(&path).unwrap().len() > 0,
                "{name} is empty"
            );
        }
        let map = PixelMap::new(vec![0.0, 1.0, 2.0, 3.0], 2).unwrap();
        save(&map, &dir.path().join("m.png")).unwrap();
    }
}
