use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use crate::sieve::load_numbers;
use crate::depth::compute_m;
use crate::stats::build_histogram;

pub fn cmd_oeis_export(n: usize, seed: Option<&PathBuf>, use_primes: bool, outdir: &PathBuf) {
    let numbers = load_numbers(n, seed, use_primes, false);
    let m_values = compute_m(&numbers);
    let hist = build_histogram(&m_values);

    fs::create_dir_all(outdir).unwrap();

    for (&level, _count) in &hist {
        let path = outdir.join(format!("oeis_m{}.txt", level));
        let file = File::create(&path).unwrap();
        let mut w = BufWriter::new(file);
        writeln!(w, "# m-class {} of prime gap depth construction", level).unwrap();
        writeln!(w, "# N = {} numbers used", n).unwrap();
        writeln!(w, "# Format: index value (1-based index within the m-class)").unwrap();
        let mut idx = 1usize;
        for (&p, &m) in numbers.iter().zip(m_values.iter()) {
            if m == level {
                writeln!(w, "{} {}", idx, p).unwrap();
                idx += 1;
            }
        }
        eprintln!("Wrote {}", path.display());
    }
}
