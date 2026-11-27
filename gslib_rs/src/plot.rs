//! Plotting utilities for GSLIB‐RS using the `plotters` crate.

use plotters::prelude::*;
use std::ffi::OsStr;
use std::path::Path;

/// Plot an empirical variogram (lags → semivariances) and optional model curve.
pub fn plot_variogram(
    lags: &[f64],
    semivs: &[f64],
    model_lags: Option<&[f64]>,
    model_semivs: Option<&[f64]>,
    out_path: &Path,
) -> Result<(), Box<dyn std::error::Error>> {
    // 1) Check inputs
    assert_eq!(
        lags.len(),
        semivs.len(),
        "lags and semivariances lengths must match"
    );

    // 2) Compute x_max, y_max from both empirical and model (if any)
    let mut all_x: Vec<f64> = lags.to_vec();
    if let Some(m) = model_lags {
        all_x.extend_from_slice(m);
    }
    let x_max = all_x.into_iter().fold(f64::NEG_INFINITY, f64::max);

    let mut all_y: Vec<f64> = semivs.to_vec();
    if let Some(m) = model_semivs {
        all_y.extend_from_slice(m);
    }
    let y_max = all_y.into_iter().fold(f64::NEG_INFINITY, f64::max);

    // 3) Choose backend
    let size = (800, 600);
    let ext = out_path
        .extension()
        .and_then(OsStr::to_str)
        .unwrap_or("")
        .to_lowercase();

    // 4a) SVG
    if ext == "svg" {
        let root = SVGBackend::new(out_path, size).into_drawing_area();
        root.fill(&WHITE)
            .map_err(Box::<dyn std::error::Error>::from)?;
        let mut chart = ChartBuilder::on(&root)
            .caption("Empirical Variogram", ("sans-serif", 24))
            .margin(20)
            .x_label_area_size(40)
            .y_label_area_size(40)
            .build_cartesian_2d(0.0..x_max, 0.0..y_max)?;
        chart
            .configure_mesh()
            .draw()
            .map_err(Box::<dyn std::error::Error>::from)?;

        chart
            .draw_series(
                lags.iter()
                    .zip(semivs.iter())
                    .map(|(&h, &g)| Circle::new((h, g), 3, RED.filled())),
            )
            .map_err(Box::<dyn std::error::Error>::from)?;

        if let (Some(ml), Some(ms)) = (model_lags, model_semivs) {
            chart
                .draw_series(LineSeries::new(
                    ml.iter().cloned().zip(ms.iter().cloned()),
                    BLUE.stroke_width(2),
                ))
                .map_err(Box::<dyn std::error::Error>::from)?;
        }

        root.present().map_err(Box::<dyn std::error::Error>::from)?;
        return Ok(());
    }

    // 4b) Bitmap
    {
        let root = BitMapBackend::new(out_path, size).into_drawing_area();
        root.fill(&WHITE)
            .map_err(Box::<dyn std::error::Error>::from)?;
        let mut chart = ChartBuilder::on(&root)
            .caption("Empirical Variogram", ("sans-serif", 24))
            .margin(20)
            .x_label_area_size(40)
            .y_label_area_size(40)
            .build_cartesian_2d(0.0..x_max, 0.0..y_max)?;
        chart
            .configure_mesh()
            .draw()
            .map_err(Box::<dyn std::error::Error>::from)?;

        chart
            .draw_series(
                lags.iter()
                    .zip(semivs.iter())
                    .map(|(&h, &g)| Circle::new((h, g), 3, RED.filled())),
            )
            .map_err(Box::<dyn std::error::Error>::from)?;

        if let (Some(ml), Some(ms)) = (model_lags, model_semivs) {
            chart
                .draw_series(LineSeries::new(
                    ml.iter().cloned().zip(ms.iter().cloned()),
                    BLUE.stroke_width(2),
                ))
                .map_err(Box::<dyn std::error::Error>::from)?;
        }

        root.present().map_err(Box::<dyn std::error::Error>::from)?;
    }

    Ok(())
}

/// Scatter plot of y vs. x. Optionally overlay regression.
pub fn plot_scatter(
    x: &[f64],
    y: &[f64],
    out_path: &Path,
    with_regression: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    assert_eq!(x.len(), y.len(), "x and y lengths must match");

    // compute data ranges
    let x_min = x.iter().cloned().fold(f64::INFINITY, f64::min);
    let x_max = x.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let y_min = y.iter().cloned().fold(f64::INFINITY, f64::min);
    let y_max = y.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

    let (x_lo, x_hi) = if (x_max - x_min).abs() < f64::EPSILON {
        (x_min - 1.0, x_max + 1.0)
    } else {
        (x_min, x_max)
    };
    let (y_lo, y_hi) = if (y_max - y_min).abs() < f64::EPSILON {
        (y_min - 1.0, y_max + 1.0)
    } else {
        (y_min, y_max)
    };

    let size = (600, 600);
    let ext = out_path
        .extension()
        .and_then(OsStr::to_str)
        .unwrap_or("")
        .to_lowercase();

    // SVG
    if ext == "svg" {
        let root = SVGBackend::new(out_path, size).into_drawing_area();
        root.fill(&WHITE)
            .map_err(Box::<dyn std::error::Error>::from)?;
        let mut chart = ChartBuilder::on(&root)
            .caption("Scatter Plot", ("sans-serif", 24))
            .margin(20)
            .x_label_area_size(40)
            .y_label_area_size(40)
            .build_cartesian_2d(x_lo..x_hi, y_lo..y_hi)?;
        chart
            .configure_mesh()
            .draw()
            .map_err(Box::<dyn std::error::Error>::from)?;
        chart
            .draw_series(
                x.iter()
                    .zip(y.iter())
                    .map(|(&xi, &yi)| Circle::new((xi, yi), 3, RED.filled())),
            )
            .map_err(Box::<dyn std::error::Error>::from)?;

        if with_regression {
            let n = x.len() as f64;
            let mean_x = x.iter().sum::<f64>() / n;
            let mean_y = y.iter().sum::<f64>() / n;
            let cov_xy = x
                .iter()
                .zip(y.iter())
                .map(|(&xi, &yi)| (xi - mean_x) * (yi - mean_y))
                .sum::<f64>();
            let var_x = x.iter().map(|&xi| (xi - mean_x).powi(2)).sum::<f64>();
            let slope = cov_xy / var_x;
            let intercept = mean_y - slope * mean_x;
            let line_pts = vec![
                (x_lo, intercept + slope * x_lo),
                (x_hi, intercept + slope * x_hi),
            ];
            chart
                .draw_series(LineSeries::new(line_pts, BLUE.stroke_width(2)))
                .map_err(Box::<dyn std::error::Error>::from)?;
        }

        root.present().map_err(Box::<dyn std::error::Error>::from)?;
        return Ok(());
    }

    // Bitmap
    {
        let root = BitMapBackend::new(out_path, size).into_drawing_area();
        root.fill(&WHITE)
            .map_err(Box::<dyn std::error::Error>::from)?;
        let mut chart = ChartBuilder::on(&root)
            .caption("Scatter Plot", ("sans-serif", 24))
            .margin(20)
            .x_label_area_size(40)
            .y_label_area_size(40)
            .build_cartesian_2d(x_lo..x_hi, y_lo..y_hi)?;
        chart
            .configure_mesh()
            .draw()
            .map_err(Box::<dyn std::error::Error>::from)?;
        chart
            .draw_series(
                x.iter()
                    .zip(y.iter())
                    .map(|(&xi, &yi)| Circle::new((xi, yi), 3, RED.filled())),
            )
            .map_err(Box::<dyn std::error::Error>::from)?;

        if with_regression {
            let n = x.len() as f64;
            let mean_x = x.iter().sum::<f64>() / n;
            let mean_y = y.iter().sum::<f64>() / n;
            let cov_xy = x
                .iter()
                .zip(y.iter())
                .map(|(&xi, &yi)| (xi - mean_x) * (yi - mean_y))
                .sum::<f64>();
            let var_x = x.iter().map(|&xi| (xi - mean_x).powi(2)).sum::<f64>();
            let slope = cov_xy / var_x;
            let intercept = mean_y - slope * mean_x;
            let line_pts = vec![
                (x_lo, intercept + slope * x_lo),
                (x_hi, intercept + slope * x_hi),
            ];
            chart
                .draw_series(LineSeries::new(line_pts, BLUE.stroke_width(2)))
                .map_err(Box::<dyn std::error::Error>::from)?;
        }

        root.present().map_err(Box::<dyn std::error::Error>::from)?;
    }

    Ok(())
}

/// Histogram plot
pub fn plot_histogram(
    data: &[f64],
    bins: usize,
    min: f64,
    max: f64,
    out_path: &Path,
) -> Result<(), Box<dyn std::error::Error>> {
    // build counts
    let mut counts = vec![0; bins];
    let span = max - min;
    for &v in data {
        let idx = if v <= min {
            0
        } else if v >= max {
            bins - 1
        } else {
            (((v - min) / span) * bins as f64).floor() as usize
        };
        counts[idx] += 1;
    }
    let y_max = *counts.iter().max().unwrap_or(&0) as f64;

    let size = (600, 400);
    let ext = out_path
        .extension()
        .and_then(OsStr::to_str)
        .unwrap_or("")
        .to_lowercase();

    // SVG
    if ext == "svg" {
        let root = SVGBackend::new(out_path, size).into_drawing_area();
        root.fill(&WHITE)
            .map_err(Box::<dyn std::error::Error>::from)?;
        let mut chart = ChartBuilder::on(&root)
            .caption("Histogram", ("sans-serif", 24))
            .margin(20)
            .x_label_area_size(40)
            .y_label_area_size(40)
            .build_cartesian_2d(min..max, 0.0..y_max)?;
        chart
            .configure_mesh()
            .draw()
            .map_err(Box::<dyn std::error::Error>::from)?;
        let bw = (max - min) / bins as f64;
        chart
            .draw_series(counts.iter().enumerate().map(|(i, &c)| {
                let x0 = min + i as f64 * bw;
                Rectangle::new([(x0, 0.0), (x0 + bw, c as f64)], BLUE.filled())
            }))
            .map_err(Box::<dyn std::error::Error>::from)?;
        root.present().map_err(Box::<dyn std::error::Error>::from)?;
        return Ok(());
    }

    // Bitmap
    {
        let root = BitMapBackend::new(out_path, size).into_drawing_area();
        root.fill(&WHITE)
            .map_err(Box::<dyn std::error::Error>::from)?;
        let mut chart = ChartBuilder::on(&root)
            .caption("Histogram", ("sans-serif", 24))
            .margin(20)
            .x_label_area_size(40)
            .y_label_area_size(40)
            .build_cartesian_2d(min..max, 0.0..y_max)?;
        chart
            .configure_mesh()
            .draw()
            .map_err(Box::<dyn std::error::Error>::from)?;
        let bw = (max - min) / bins as f64;
        chart
            .draw_series(counts.iter().enumerate().map(|(i, &c)| {
                let x0 = min + i as f64 * bw;
                Rectangle::new([(x0, 0.0), (x0 + bw, c as f64)], BLUE.filled())
            }))
            .map_err(Box::<dyn std::error::Error>::from)?;
        root.present().map_err(Box::<dyn std::error::Error>::from)?;
    }

    Ok(())
}

/// Pixel map
pub fn plot_pixelmap<Fc>(
    data: &[f64],
    nx: usize,
    ny: usize,
    out_path: &Path,
    palette: Fc,
) -> Result<(), Box<dyn std::error::Error>>
where
    Fc: Fn(f64) -> RGBColor,
{
    assert_eq!(data.len(), nx * ny, "data.len must be nx*ny");

    let size = (800, 800);
    let ext = out_path
        .extension()
        .and_then(OsStr::to_str)
        .unwrap_or("")
        .to_lowercase();

    // SVG
    if ext == "svg" {
        let root = SVGBackend::new(out_path, size).into_drawing_area();
        root.fill(&WHITE)
            .map_err(Box::<dyn std::error::Error>::from)?;
        let x_scale = size.0 as f64 / nx as f64;
        let y_scale = size.1 as f64 / ny as f64;
        for ix in 0..nx {
            for iy in 0..ny {
                let v = data[ix + iy * nx];
                let c = palette(v);
                let x0 = (ix as f64 * x_scale) as i32;
                let y0 = (iy as f64 * y_scale) as i32;
                let x1 = ((ix + 1) as f64 * x_scale) as i32;
                let y1 = ((iy + 1) as f64 * y_scale) as i32;
                root.draw(&Rectangle::new([(x0, y0), (x1, y1)], c.filled()))
                    .map_err(Box::<dyn std::error::Error>::from)?;
            }
        }
        root.present().map_err(Box::<dyn std::error::Error>::from)?;
        return Ok(());
    }

    // Bitmap
    {
        let root = BitMapBackend::new(out_path, size).into_drawing_area();
        root.fill(&WHITE)
            .map_err(Box::<dyn std::error::Error>::from)?;
        let x_scale = size.0 as f64 / nx as f64;
        let y_scale = size.1 as f64 / ny as f64;
        for ix in 0..nx {
            for iy in 0..ny {
                let v = data[ix + iy * nx];
                let c = palette(v);
                let x0 = (ix as f64 * x_scale) as i32;
                let y0 = (iy as f64 * y_scale) as i32;
                let x1 = ((ix + 1) as f64 * x_scale) as i32;
                let y1 = ((iy + 1) as f64 * y_scale) as i32;
                root.draw(&Rectangle::new([(x0, y0), (x1, y1)], c.filled()))
                    .map_err(Box::<dyn std::error::Error>::from)?;
            }
        }
        root.present().map_err(Box::<dyn std::error::Error>::from)?;
    }

    Ok(())
}
