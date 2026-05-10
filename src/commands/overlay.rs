use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use crate::sieve::load_numbers;
use crate::depth::compute_m;
use crate::stats::{loglog_slope, invert_curve_at_y};

pub fn cmd_overlay(n: usize, seed: Option<&PathBuf>, outdir: &PathBuf) {
    eprintln!("Loading {} numbers...", n);
    let numbers = load_numbers(n, seed, false);
    eprintln!("Computing m-values...");
    let m_values = compute_m(&numbers);

    let max_m = *m_values.iter().max().unwrap_or(&0);

    let mut classes: Vec<Vec<usize>> = vec![Vec::new(); max_m as usize + 1];
    for (i, &m) in m_values.iter().enumerate() {
        classes[m as usize].push(i + 1);
    }

    // Per class, build a downsampled log-log curve at geometric k = 1, 2, 4, ..., len.
    let log_curves: Vec<Vec<(f64, f64)>> = classes
        .iter()
        .map(|cls| {
            if cls.len() < 8 {
                return Vec::new();
            }
            let mut pts = Vec::new();
            let mut k = 1usize;
            while k <= cls.len() {
                pts.push(((k as f64).ln(), (cls[k - 1] as f64).ln()));
                k = (k * 2).max(k + 1);
            }
            let last = cls.len();
            let last_ln = (last as f64).ln();
            if pts.last().map(|p| p.0).unwrap_or(0.0) < last_ln {
                pts.push((last_ln, (cls[last - 1] as f64).ln()));
            }
            pts
        })
        .collect();

    let mut fits: Vec<Option<(f64, f64)>> = Vec::with_capacity(log_curves.len());
    for curve in &log_curves {
        if curve.len() < 3 {
            fits.push(None);
            continue;
        }
        let slope = loglog_slope(curve);
        let mean_x = curve.iter().map(|p| p.0).sum::<f64>() / curve.len() as f64;
        let mean_y = curve.iter().map(|p| p.1).sum::<f64>() / curve.len() as f64;
        let intercept = mean_y - slope * mean_x;
        fits.push(Some((slope, intercept)));
    }

    println!("\nClass log-log fits  log(pi) = slope * log(k) + intercept\n");
    println!("{:<6} {:>10} {:>10} {:>12}", "m", "size", "slope", "intercept");
    println!("{}", "-".repeat(42));
    for (level, fit) in fits.iter().enumerate() {
        let size = classes[level].len();
        match fit {
            Some((s, b)) => println!("{:<6} {:>10} {:>10.4} {:>12.4}", level, size, s, b),
            None => println!("{:<6} {:>10} {:>10} {:>12}", level, size, "—", "—"),
        }
    }

    println!("\nHorizontal-shift test between consecutive classes\n");
    println!(
        "If 'parallel lines' holds, dx_y should be ~constant across the overlap range,\n\
         and the residual after subtracting the mean shift should be small.\n"
    );
    println!(
        "{:<10} {:>14} {:>14} {:>14} {:>14}",
        "pair", "mean dx", "stddev dx", "y-overlap", "n samples"
    );
    println!("{}", "-".repeat(70));

    fs::create_dir_all(outdir).unwrap();
    let path = outdir.join("overlay_shifts.tsv");
    let mut w = BufWriter::new(File::create(&path).unwrap());
    writeln!(w, "m_lo\tm_hi\ty\tx_lo\tx_hi\tdx").unwrap();

    for level in 0..(log_curves.len().saturating_sub(1)) {
        let lo = &log_curves[level];
        let hi = &log_curves[level + 1];
        if lo.len() < 3 || hi.len() < 3 {
            continue;
        }

        let y_lo_min = lo.first().unwrap().1;
        let y_lo_max = lo.last().unwrap().1;
        let y_hi_min = hi.first().unwrap().1;
        let y_hi_max = hi.last().unwrap().1;
        let y_start = y_lo_min.max(y_hi_min);
        let y_end = y_lo_max.min(y_hi_max);
        if y_end <= y_start {
            continue;
        }

        let n_samples = 50usize;
        let mut shifts = Vec::with_capacity(n_samples);
        for i in 0..n_samples {
            let y = y_start + (y_end - y_start) * (i as f64) / ((n_samples - 1) as f64);
            let x_lo = invert_curve_at_y(lo, y);
            let x_hi = invert_curve_at_y(hi, y);
            if let (Some(xl), Some(xh)) = (x_lo, x_hi) {
                let dx = xh - xl;
                shifts.push((y, xl, xh, dx));
                writeln!(w, "{}\t{}\t{}\t{}\t{}\t{}", level, level + 1, y, xl, xh, dx).unwrap();
            }
        }

        if shifts.is_empty() {
            continue;
        }
        let nf = shifts.len() as f64;
        let mean: f64 = shifts.iter().map(|(_, _, _, d)| d).sum::<f64>() / nf;
        let var: f64 = shifts.iter().map(|(_, _, _, d)| (d - mean).powi(2)).sum::<f64>() / nf;
        let std = var.sqrt();
        let pair = format!("m{}->m{}", level, level + 1);
        let yspan = format!("{:.2}-{:.2}", y_start, y_end);
        println!(
            "{:<10} {:>14.4} {:>14.4} {:>14} {:>14}",
            pair, mean, std, yspan, shifts.len()
        );
    }

    eprintln!("Wrote {}", path.display());

    println!(
        "\nInterpretation:\n  - If stddev/mean is small (<<1), classes are well-described as horizontal shifts\n    of each other on log-log axes — the 'parallel lines' / pure-thinning hypothesis.\n  - If stddev/mean is large, the shapes differ and the construction is doing more than thinning."
    );
}
