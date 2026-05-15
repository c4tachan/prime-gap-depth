use std::collections::{BTreeMap, HashSet};
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use crate::sieve::load_numbers;
use crate::depth::compute_m;
use crate::stats::{build_histogram, loglog_slope};

pub fn cmd_growth(seed: Option<&PathBuf>, from_generator: bool, outdir: &PathBuf) {
    let ns: Vec<usize> = vec![100, 1_000, 10_000, 100_000, 1_000_000, 10_000_000];
    println!("Growth analysis across N values\n");

    let mut all_hists: Vec<(usize, BTreeMap<u32, usize>)> = Vec::new();
    for &n in &ns {
        let numbers = load_numbers(n, seed, from_generator, false);
        let m_values = compute_m::<u32>(&numbers);
        let hist = build_histogram(&m_values);
        let max_m = hist.keys().max().copied().unwrap_or(0);
        eprint!("N={:>10}  max_m={}  ", n, max_m);
        for (&m, &cnt) in &hist {
            eprint!("m{}={} ", m, cnt);
        }
        eprintln!();
        all_hists.push((n, hist));
    }

    let mut all_levels: Vec<u32> = all_hists
        .iter()
        .flat_map(|(_, h)| h.keys().copied())
        .collect::<HashSet<_>>()
        .into_iter()
        .collect();
    all_levels.sort_unstable();

    let col_w = 12;
    print!("{:<12}", "N");
    for &m in &all_levels {
        print!("{:>width$}", format!("m={}", m), width = col_w);
    }
    println!();
    println!("{}", "-".repeat(12 + col_w * all_levels.len()));
    for (n, hist) in &all_hists {
        print!("{:<12}", n);
        for &m in &all_levels {
            let cnt = hist.get(&m).copied().unwrap_or(0);
            print!("{:>width$}", cnt, width = col_w);
        }
        println!();
    }

    fs::create_dir_all(outdir).unwrap();
    let path = outdir.join("growth.tsv");
    let file = File::create(&path).unwrap();
    let mut w = BufWriter::new(file);
    let hdrs: Vec<String> = all_levels.iter().map(|m| format!("m={}", m)).collect();
    writeln!(w, "N\t{}", hdrs.join("\t")).unwrap();
    for (n, hist) in &all_hists {
        let vals: Vec<String> = all_levels
            .iter()
            .map(|m| hist.get(m).copied().unwrap_or(0).to_string())
            .collect();
        writeln!(w, "{}\t{}", n, vals.join("\t")).unwrap();
    }
    eprintln!("Wrote {}", path.display());

    println!("\nEstimated growth exponents (log-log linear fit count ~ N^alpha):");
    for &m in &all_levels {
        let points: Vec<(f64, f64)> = all_hists
            .iter()
            .filter_map(|(n, hist)| {
                let cnt = *hist.get(&m)?;
                if cnt == 0 { None } else { Some(((*n as f64).ln(), (cnt as f64).ln())) }
            })
            .collect();
        if points.len() >= 2 {
            let alpha = loglog_slope(&points);
            println!("  m={}: alpha ~ {:.4}", m, alpha);
        }
    }
}
