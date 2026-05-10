use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use crate::sieve::load_numbers;
use crate::depth::compute_m;
use crate::stats::chi2_p_value_approx;

pub fn cmd_mod_residue(n: usize, seed: Option<&PathBuf>, outdir: &PathBuf, modulus: u64) {
    let numbers = load_numbers(n, seed, false);
    let m_values = compute_m(&numbers);

    // Baseline: all numbers mod modulus
    let mut baseline_counts: HashMap<u64, usize> = HashMap::new();
    for &p in &numbers {
        *baseline_counts.entry(p % modulus).or_insert(0) += 1;
    }
    let total = numbers.len() as f64;

    let max_m = *m_values.iter().max().unwrap_or(&0);

    // Residues that actually appear in the seed set
    let mut residues: Vec<u64> = baseline_counts.keys().copied().collect();
    residues.sort_unstable();

    println!("Mod-residue analysis (mod {})", modulus);
    println!("m-level counts and residue distributions vs baseline\n");

    fs::create_dir_all(outdir).unwrap();
    let path = outdir.join(format!("mod_{}.tsv", modulus));
    let file = File::create(&path).unwrap();
    let mut w = BufWriter::new(file);

    let header_parts: Vec<String> = residues.iter().map(|r| format!("r={}", r)).collect();
    writeln!(w, "m\tcount\tchi2_p\t{}", header_parts.join("\t")).unwrap();

    for level in 0..=max_m {
        let level_nums: Vec<u64> = numbers
            .iter()
            .zip(m_values.iter())
            .filter(|(_, &m)| m == level)
            .map(|(&p, _)| p)
            .collect();
        let level_n = level_nums.len() as f64;
        if level_n == 0.0 {
            continue;
        }

        let mut level_counts: HashMap<u64, usize> = HashMap::new();
        for p in &level_nums {
            *level_counts.entry(p % modulus).or_insert(0) += 1;
        }

        // Chi-squared: sum((obs - exp)^2 / exp) where exp = level_n * baseline_frac
        let mut chi2 = 0.0f64;
        let mut df = 0usize;
        for &r in &residues {
            let obs = *level_counts.get(&r).unwrap_or(&0) as f64;
            let base_frac = *baseline_counts.get(&r).unwrap_or(&0) as f64 / total;
            let exp = level_n * base_frac;
            if exp > 0.0 {
                chi2 += (obs - exp).powi(2) / exp;
                df += 1;
            }
        }
        let p_val = chi2_p_value_approx(chi2, df.saturating_sub(1));

        if level <= 5 {
            println!("m={}: n={:.0}  chi2={:.2}  p~{:.4}", level, level_n, chi2, p_val);
        }

        let pcts: Vec<String> = residues
            .iter()
            .map(|r| {
                let cnt = *level_counts.get(r).unwrap_or(&0) as f64;
                format!("{:.4}", cnt / level_n * 100.0)
            })
            .collect();
        writeln!(w, "{}\t{}\t{:.6}\t{}", level, level_n as usize, p_val, pcts.join("\t")).unwrap();
    }

    eprintln!("Wrote {}", path.display());
}
