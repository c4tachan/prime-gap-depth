use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use crate::sieve::load_numbers;
use crate::depth::compute_pi_chain;
use crate::stats::{build_histogram, print_histogram};

pub fn cmd_pi_chain(n: usize, seed: Option<&PathBuf>, use_primes: bool, outdir: &PathBuf) {
    eprintln!("Loading {} numbers...", n);
    let primes = load_numbers(n, seed, use_primes, false);
    eprintln!("Computing pi-chain depths...");
    let depths = compute_pi_chain(&primes);

    let hist = build_histogram(&depths);
    let max_d = *depths.iter().max().unwrap_or(&0);

    println!("\n=== Pi-chain depth: family counts ===");
    print_histogram(&hist);

    println!("\n=== First appearances (smallest prime achieving each pi-chain depth) ===");
    println!("{:<6} {:>16} {:>16}", "m", "pi (index)", "p");
    println!("{}", "-".repeat(42));
    for level in 0..=max_d {
        if let Some(pos) = depths.iter().position(|&d| d == level) {
            println!("{:<6} {:>16} {:>16}", level, pos + 1, primes[pos]);
        }
    }

    // Per-class ordered list of (1-based prime index, prime value).
    let mut classes: Vec<Vec<(usize, u64)>> = vec![Vec::new(); max_d as usize + 1];
    for (i, (&p, &d)) in primes.iter().zip(depths.iter()).enumerate() {
        classes[d as usize].push((i + 1, p));
    }

    fs::create_dir_all(outdir).unwrap();
    let path = outdir.join("pichain_C.tsv");
    let mut w = BufWriter::new(File::create(&path).unwrap());
    writeln!(w, "m\tk\tp_k\tC").unwrap();

    let ordinals: [usize; 5] = [10, 100, 1_000, 10_000, 100_000];

    println!("\n=== C(m, k) = p_k / k^2 across each pi-chain depth class ===");
    println!("(Tracks how p_k scales with k within a depth-m family.)\n");
    println!("{:<6} {:>10} {:>16} {:>22}", "m", "k", "p_k", "C = p_k / k^2");
    println!("{}", "-".repeat(58));

    for level in 1..=max_d {
        let cls = &classes[level as usize];
        if cls.is_empty() {
            continue;
        }
        let mut emitted_last = false;
        for &k in &ordinals {
            if k <= cls.len() {
                let (_idx, pv) = cls[k - 1];
                let c = pv as f64 / (k as f64).powi(2);
                writeln!(w, "{}\t{}\t{}\t{}", level, k, pv, c).unwrap();
                println!("{:<6} {:>10} {:>16} {:>22}", level, k, pv, format!("{}", c));
                if k == cls.len() {
                    emitted_last = true;
                }
            }
        }
        // Always include the last element of the class, so callers see
        // the asymptotic-most C value even if the class is small or large.
        if !emitted_last {
            let k_last = cls.len();
            let (_idx, pv) = cls[k_last - 1];
            let c = pv as f64 / (k_last as f64).powi(2);
            writeln!(w, "{}\t{}\t{}\t{}", level, k_last, pv, c).unwrap();
            println!("{:<6} {:>10} {:>16} {:>22}", level, k_last, pv, format!("{}", c));
        }
        println!();
    }
    eprintln!("Wrote {}", path.display());

    // Ratios: between consecutive pi-chain depth classes, compare last-C values.
    println!("=== C-ratios between consecutive depth classes (last entry of each class) ===");
    println!("{:<10} {:>20} {:>20} {:>14}", "pair", "C(m, last)", "C(m+1, last)", "ratio");
    println!("{}", "-".repeat(66));
    for level in 1..max_d {
        let cls_a = &classes[level as usize];
        let cls_b = &classes[(level + 1) as usize];
        if cls_a.is_empty() || cls_b.is_empty() {
            continue;
        }
        let ka = cls_a.len() as f64;
        let kb = cls_b.len() as f64;
        let pa = cls_a.last().unwrap().1 as f64;
        let pb = cls_b.last().unwrap().1 as f64;
        let ca = pa / ka.powi(2);
        let cb = pb / kb.powi(2);
        let ratio = if cb > 0.0 { cb / ca } else { f64::NAN };
        let pair = format!("m{}->m{}", level, level + 1);
        println!("{:<10} {:>20.6e} {:>20.6e} {:>14.4}", pair, ca, cb, ratio);
    }
}
