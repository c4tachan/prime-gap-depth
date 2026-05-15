use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

use crate::sieve::load_numbers;
use crate::depth::compute_m;
use crate::stats::build_histogram;

// ---------------------------------------------------------------------------
// Helper: count data rows (lines minus header) without allocating
// ---------------------------------------------------------------------------

fn count_data_lines(path: &Path) -> u64 {
    let mut file = File::open(path).unwrap_or_else(|e| {
        eprintln!("Cannot open {}: {}", path.display(), e);
        std::process::exit(1);
    });
    let mut buf = vec![0u8; 65536];
    let mut newlines: u64 = 0;
    loop {
        let n = file.read(&mut buf).unwrap();
        if n == 0 { break; }
        newlines += buf[..n].iter().filter(|&&b| b == b'\n').count() as u64;
    }
    newlines.saturating_sub(1) // subtract header line
}

// ---------------------------------------------------------------------------
// Export from an existing results.csv (streams — O(1) memory)
// ---------------------------------------------------------------------------

pub fn cmd_oeis_export_from_csv(results_path: &Path, outdir: &PathBuf) {
    eprintln!("Scanning {} for row count ...", results_path.display());
    let total_n = count_data_lines(results_path);
    eprintln!("  {} data rows.", total_n);

    eprintln!("Streaming {} ...", results_path.display());
    let file = File::open(results_path).unwrap_or_else(|e| {
        eprintln!("Cannot open {}: {}", results_path.display(), e);
        std::process::exit(1);
    });
    let reader = BufReader::with_capacity(8 * 1024 * 1024, file);

    fs::create_dir_all(outdir).unwrap();

    const MAX_M: usize = 10;
    let mut writers: Vec<Option<BufWriter<File>>> = (0..MAX_M).map(|_| None).collect();
    let mut counters: Vec<u64> = vec![0u64; MAX_M];
    let mut row_count: u64 = 0;

    let mut lines = reader.lines();
    lines.next(); // skip header

    for line in lines {
        let line = match line {
            Ok(l)  => l,
            Err(e) => { eprintln!("Read error: {}", e); break; }
        };

        // Format: index,value,m_gap[,m_pichain]
        let mut parts = line.splitn(4, ',');
        let _index: u64 = match parts.next().and_then(|s| s.parse().ok()) {
            Some(v) => v, None => continue,
        };
        let value: u64 = match parts.next().and_then(|s| s.parse().ok()) {
            Some(v) => v, None => continue,
        };
        let m: usize = match parts.next().and_then(|s| s.trim().parse().ok()) {
            Some(v) => v, None => continue,
        };

        if m < MAX_M {
            counters[m] += 1;
            let idx = counters[m];
            let path_clone = outdir.join(format!("oeis_m{}.txt", m));
            let w = writers[m].get_or_insert_with(|| {
                let f = File::create(&path_clone).unwrap();
                let mut bw = BufWriter::new(f);
                writeln!(bw, "# m-class {} of gap-depth construction", m).unwrap();
                writeln!(bw, "# N = {} numbers used", total_n).unwrap();
                writeln!(bw, "# Format: index value (1-based index within the m-class)").unwrap();
                bw
            });
            writeln!(w, "{} {}", idx, value).unwrap();
        }

        row_count += 1;
        if row_count % 25_000_000 == 0 {
            eprint!("\r  {:.0}M rows processed ...", row_count as f64 / 1_000_000.0);
        }
    }
    eprintln!("\r  {} rows processed.              ", row_count);

    for (m, w_opt) in writers.iter_mut().enumerate() {
        if let Some(w) = w_opt {
            w.flush().unwrap();
            eprintln!("Wrote oeis_m{}.txt  ({} entries)", m, counters[m]);
        }
    }
}

// ---------------------------------------------------------------------------
// Export by recomputing from scratch
// ---------------------------------------------------------------------------

pub fn cmd_oeis_export(n: usize, seed: Option<&PathBuf>, from_generator: bool, outdir: &PathBuf) {
    let numbers = load_numbers(n, seed, from_generator, false);
    let m_values = compute_m::<u32>(&numbers);
    let hist = build_histogram(&m_values);

    fs::create_dir_all(outdir).unwrap();

    for (&level, _count) in &hist {
        let path = outdir.join(format!("oeis_m{}.txt", level));
        let file = File::create(&path).unwrap();
        let mut w = BufWriter::new(file);
        writeln!(w, "# m-class {} of gap-depth construction", level).unwrap();
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
