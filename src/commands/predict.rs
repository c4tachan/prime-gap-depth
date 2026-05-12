use std::path::PathBuf;

use crate::sieve::load_numbers;
use crate::depth::compute_m;
use crate::stats::{quadratic_fit, linear_fit};

pub fn cmd_predict(n: usize, seed: Option<&PathBuf>, from_generator: bool, m_min: u32, m_max: u32) {
    eprintln!("Loading {} numbers...", n);
    let numbers = load_numbers(n, seed, from_generator, false);
    eprintln!("Computing m-values...");
    let m_values = compute_m(&numbers);

    let observed_max_m = *m_values.iter().max().unwrap_or(&0);

    // first_idx[m] = 1-based index of the smallest element with that m-value.
    // numbers is sorted ascending, so the first occurrence is the smallest.
    let mut first_idx: Vec<Option<usize>> = vec![None; (observed_max_m as usize) + 1];
    for (i, &m) in m_values.iter().enumerate() {
        let slot = &mut first_idx[m as usize];
        if slot.is_none() {
            *slot = Some(i + 1);
        }
    }

    let mut data: Vec<(f64, f64)> = Vec::new(); // (m, log_idx) for the fit
    println!("\nFirst-occurrence indices (observed):\n");
    println!("{:<6} {:>20} {:>14}", "m", "first index", "ln(index)");
    println!("{}", "-".repeat(44));
    for m in 0..=observed_max_m {
        if let Some(idx) = first_idx[m as usize] {
            let lnv = (idx as f64).ln();
            println!("{:<6} {:>20} {:>14.4}", m, idx, lnv);
            if m >= m_min {
                data.push((m as f64, lnv));
            }
        }
    }

    if data.len() < 3 {
        println!(
            "\nNeed at least 3 data points with m >= {} to fit a quadratic. Got {}.",
            m_min, data.len()
        );
        println!("Hint: rerun with a larger -n so more m-classes appear.");
        return;
    }

    let (a, b, c) = quadratic_fit(&data);
    let (b1, c1) = linear_fit(&data);

    println!("\nQuadratic fit  ln(first_idx) = a*m^2 + b*m + c");
    println!("  a = {:>12.6}", a);
    println!("  b = {:>12.6}", b);
    println!("  c = {:>12.6}", c);
    let mut sse_q = 0.0;
    let mut sse_l = 0.0;
    for &(m, y) in &data {
        sse_q += (y - (a * m * m + b * m + c)).powi(2);
        sse_l += (y - (b1 * m + c1)).powi(2);
    }
    println!("  SSE (quadratic) = {:.6}", sse_q);
    println!("\nLinear fit (for reference)  ln(first_idx) = b*m + c");
    println!("  b = {:>12.6}", b1);
    println!("  c = {:>12.6}", c1);
    println!("  SSE (linear)    = {:.6}", sse_l);

    println!("\nResiduals at observed m (m >= {}):", m_min);
    println!("{:<6} {:>14} {:>14} {:>14}", "m", "ln(observed)", "quad pred", "residual");
    println!("{}", "-".repeat(50));
    for &(m, y) in &data {
        let pred = a * m * m + b * m + c;
        println!("{:<6} {:>14.4} {:>14.4} {:>14.4}", m as u32, y, pred, y - pred);
    }

    println!("\nProjections for higher m:");
    println!(
        "{:<6} {:>14} {:>26} {:>26}",
        "m", "ln(idx) pred", "predicted element index", "approx element value"
    );
    println!("{}", "-".repeat(76));
    for m in m_min..=m_max {
        let mf = m as f64;
        let ln_pred = a * mf * mf + b * mf + c;
        let idx_pred = ln_pred.exp();
        let p_approx = idx_pred * idx_pred.ln().max(1.0);
        let observed_marker = if m <= observed_max_m { " (obs)" } else { "" };
        let idx_str = if idx_pred.is_finite() && idx_pred < 1e30 {
            format!("{:.3e}", idx_pred)
        } else {
            "overflow".to_string()
        };
        let p_str = if p_approx.is_finite() && p_approx < 1e35 {
            format!("{:.3e}", p_approx)
        } else {
            "overflow".to_string()
        };
        println!(
            "{:<6} {:>14.4} {:>26} {:>26}{}",
            m, ln_pred, idx_str, p_str, observed_marker
        );
    }

    println!(
        "\nNotes:\n  - 'predicted element index' is what the quadratic extrapolation gives for the smallest\n    element in that m-class. Compare against observed values where available.\n  - 'approx element value' uses the PNT-style estimate p ~ n * ln(n) as a rough scale guide.\n  - Confidence is highest near the fit range (m={}..{}). Each step beyond that compounds error."
        , m_min, observed_max_m
    );
}
