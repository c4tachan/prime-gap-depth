use std::fs;
use std::io::Write;
use std::path::PathBuf;

use crate::sieve::sieve_first_n;
use crate::depth::compute_m;

/// Find the minimum window size W such that, when computing m-values on the
/// suffix primes[n-W..n], every prime in that window matches its global m-value.
///
/// Algorithm: start with W=1 anchored at the end. Step backwards one prime at
/// a time. If the new prime doesn't match with the current window, expand W
/// until it does — then re-verify all primes already in the window with the
/// new W. Keep expanding until the whole window is consistent. Stop if W
/// reaches the index of the current prime (full prefix needed).
pub fn cmd_locality(n: usize, outdir: &PathBuf) {
    eprintln!("Loading {} primes...", n);
    let primes = sieve_first_n(n);
    eprintln!("Computing global m-values...");
    let global_m = compute_m(&primes);

    let out_path = outdir.join("locality.csv");
    fs::create_dir_all(outdir).expect("cannot create output directory");
    let mut file = fs::File::create(&out_path).expect("cannot create locality.csv");
    writeln!(file, "index,prime,global_m,window_at_acceptance").unwrap();

    println!();
    println!("{:>10}  {:>14}  {:>8}  {:>10}", "index", "prime", "global_m", "window");
    println!("{}", "-".repeat(48));

    // `w` is the current window size — the suffix primes[n-w..n].
    // Invariant: all primes in the current window agree with global_m
    // when computed together on primes[n-w..n].
    let mut w = 0usize;
    // Track window size at the point each prime was accepted, in reverse order.
    let mut accepted: Vec<(usize, usize)> = Vec::with_capacity(n); // (global_index, window)
    let mut stopped_at: Option<usize> = None;

    'outer: for target in (0..n).rev() {
        // Try expanding the window until the entire suffix primes[n-w..n]
        // (which now includes `target`) is consistent.
        loop {
            w += 1;
            let lo = n - w;
            // w must not exceed target+1 (i.e. lo must not go past target).
            // Since we're adding target, lo == target when w == n - target.
            if lo < target {
                // Shouldn't happen in normal flow, but guard anyway.
                w = n - target;
            }

            let window = &primes[n - w..n];
            let local_m = compute_m(window);

            // Check every prime in the window against global.
            let all_ok = (0..w).all(|i| local_m[i] == global_m[n - w + i]);

            if all_ok {
                // Record acceptance for each prime newly added (just `target` here,
                // but re-verifications of earlier primes don't need re-recording).
                println!("{:>10}  {:>14}  {:>8}  {:>10}", target, primes[target], global_m[target], w);
                writeln!(file, "{},{},{},{}", target, primes[target], global_m[target], w).unwrap();
                accepted.push((target, w));
                break;
            }

            // Window not yet consistent — but check stopping condition before expanding.
            if w >= target + 1 {
                println!("{:>10}  {:>14}  {:>8}  {:>10}  (full set needed, stopping)",
                    target, primes[target], global_m[target], w);
                writeln!(file, "{},{},{},{}", target, primes[target], global_m[target], w).unwrap();
                accepted.push((target, w));
                stopped_at = Some(target);
                break 'outer;
            }
        }
    }

    // Summary
    let windows: Vec<usize> = accepted.iter().map(|&(_, w)| w).collect();
    let count = windows.len();
    if count == 0 {
        println!("\nNo data.");
        return;
    }
    let mut sorted = windows.clone();
    sorted.sort_unstable();
    let mean = sorted.iter().sum::<usize>() as f64 / count as f64;
    let median = sorted[count / 2];
    let p90 = sorted[(count as f64 * 0.90) as usize];
    let p99 = sorted[(count as f64 * 0.99) as usize];
    let max = *sorted.last().unwrap();

    println!();
    println!("=== Summary ({} primes examined) ===", count);
    println!("  mean window : {:.1}", mean);
    println!("  median      : {}", median);
    println!("  p90         : {}", p90);
    println!("  p99         : {}", p99);
    println!("  max         : {}", max);
    if let Some(idx) = stopped_at {
        println!("  stopped at index {} (full prefix required)", idx);
    } else {
        println!("  completed all {} primes without hitting full-set condition", n);
    }
    println!();
    println!("Output written to {}", out_path.display());
}
