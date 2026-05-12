use std::fs;
use std::io::Write;
use std::path::PathBuf;

use crate::sieve::load_numbers;
use crate::depth::compute_m;

/// Scan from the last prime backwards, maintaining a growing window size `w`.
///
/// For each prime at index `i`, compute m on the local window
/// `primes[i+1-w ..= i]` and check whether the last element's local m matches
/// `global_m[i]`. If it does, record the acceptance window and move to `i-1`.
/// If it doesn't, increment `w` and restart the scan from `i = n-1`.
///
/// Each prime's recorded window is `w` at the time it was first accepted.
pub fn cmd_locality(n: usize, seed: Option<&PathBuf>, use_primes: bool, outdir: &PathBuf) {
    eprintln!("Loading {} numbers...", n);
    let primes = load_numbers(n, seed, use_primes, false);
    eprintln!("Computing global m-values...");
    let global_m = compute_m(&primes);

    let out_path = outdir.join("locality.csv");
    fs::create_dir_all(outdir).expect("cannot create output directory");
    let mut file = fs::File::create(&out_path).expect("cannot create locality.csv");
    writeln!(file, "index,prime,global_m,window_at_acceptance").unwrap();

    // window_at[i] = w when primes[i] was first accepted; 0 = not yet accepted.
    let mut window_at = vec![0usize; n];

    let mut w = 1usize;
    let mut i = n - 1;

    loop {
        let lo = i.saturating_sub(w - 1);
        let window = &primes[lo..=i];
        let local_m = compute_m(window);

        if local_m[window.len() - 1] == global_m[i] {
            if window_at[i] == 0 {
                window_at[i] = w;
            }
            if i == 0 {
                break;
            }
            i -= 1;
        } else {
            w += 1;
            if w > n {
                // Give up — record remaining unaccepted primes as needing full set.
                for j in 0..=i {
                    if window_at[j] == 0 {
                        window_at[j] = n;
                    }
                }
                break;
            }
            i = n - 1;
        }
    }

    println!();
    println!("{:>10}  {:>14}  {:>8}  {:>10}", "index", "prime", "global_m", "window");
    println!("{}", "-".repeat(48));

    for idx in 0..n {
        println!("{:>10}  {:>14}  {:>8}  {:>10}", idx, primes[idx], global_m[idx], window_at[idx]);
        writeln!(file, "{},{},{},{}", idx, primes[idx], global_m[idx], window_at[idx]).unwrap();
    }

    let mut sorted = window_at.clone();
    sorted.sort_unstable();
    let mean = sorted.iter().sum::<usize>() as f64 / n as f64;
    let median = sorted[n / 2];
    let p90 = sorted[(n as f64 * 0.90) as usize];
    let p99 = sorted[(n as f64 * 0.99) as usize];
    let max = *sorted.last().unwrap();

    println!();
    println!("=== Summary ({} primes examined) ===", n);
    println!("  mean window : {:.1}", mean);
    println!("  median      : {}", median);
    println!("  p90         : {}", p90);
    println!("  p99         : {}", p99);
    println!("  max         : {}", max);
    println!();
    println!("Output written to {}", out_path.display());
}

