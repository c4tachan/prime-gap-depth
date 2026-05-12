use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use crate::sieve::load_numbers;
use crate::depth::compute_m;
use crate::stats::loglog_slope;

pub fn cmd_class_quantiles(n: usize, seed: Option<&PathBuf>, from_generator: bool, outdir: &PathBuf) {
    eprintln!("Loading {} numbers...", n);
    let numbers = load_numbers(n, seed, from_generator, false);
    eprintln!("Computing m-values...");
    let m_values = compute_m(&numbers);

    let max_m = *m_values.iter().max().unwrap_or(&0);
    let mut classes: Vec<Vec<(usize, u64)>> = vec![Vec::new(); max_m as usize + 1];
    for (i, (&p, &m)) in numbers.iter().zip(m_values.iter()).enumerate() {
        classes[m as usize].push((i + 1, p));
    }

    let ordinals: [usize; 10] = [1, 2, 5, 10, 100, 1_000, 10_000, 100_000, 1_000_000, 10_000_000];

    println!("\nFor each m-level, the k-th element to reach that level.");
    println!("Columns are 'index (value)' — index is the 1-based rank in the input sequence.\n");

    print!("{:<6} {:>10}", "m", "size");
    for &k in &ordinals {
        print!("  {:>22}", format!("k={}", k));
    }
    println!();
    println!("{}", "-".repeat(6 + 10 + ordinals.len() * 24));

    for level in 0..=max_m {
        let cls = &classes[level as usize];
        print!("{:<6} {:>10}", level, cls.len());
        for &k in &ordinals {
            if k <= cls.len() {
                let (pi, p) = cls[k - 1];
                print!("  {:>22}", format!("{} ({})", pi, p));
            } else {
                print!("  {:>22}", "—");
            }
        }
        println!();
    }

    fs::create_dir_all(outdir).unwrap();
    let path = outdir.join("class_quantiles.tsv");
    let mut w = BufWriter::new(File::create(&path).unwrap());
    writeln!(w, "m\tsize\tk\telement_index\telement_value").unwrap();
    for level in 0..=max_m {
        let cls = &classes[level as usize];
        for &k in &ordinals {
            if k <= cls.len() {
                let (pi, p) = cls[k - 1];
                writeln!(w, "{}\t{}\t{}\t{}\t{}", level, cls.len(), k, pi, p).unwrap();
            }
        }
    }
    eprintln!("Wrote {}", path.display());

    println!("\nLog-log fit  log(element_index) vs log(k)  within each class:");
    println!("{:<6} {:>8} {:>12}", "m", "slope", "intercept");
    println!("{}", "-".repeat(30));
    for level in 0..=max_m {
        let cls = &classes[level as usize];
        if cls.len() < 5 {
            continue;
        }
        let pts: Vec<(f64, f64)> = (1..=cls.len())
            .filter(|&k| k.is_power_of_two() || k == 1 || k % 100 == 0 || k == cls.len())
            .take(2000)
            .map(|k| ((k as f64).ln(), (cls[k - 1].0 as f64).ln()))
            .collect();
        if pts.len() < 2 {
            continue;
        }
        let slope = loglog_slope(&pts);
        let mean_x = pts.iter().map(|p| p.0).sum::<f64>() / pts.len() as f64;
        let mean_y = pts.iter().map(|p| p.1).sum::<f64>() / pts.len() as f64;
        let intercept = mean_y - slope * mean_x;
        println!("{:<6} {:>8.4} {:>12.4}", level, slope, intercept);
    }
}
