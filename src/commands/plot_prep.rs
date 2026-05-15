use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

// ---------------------------------------------------------------------------
// Minimal xorshift64 — no external dependency needed
// ---------------------------------------------------------------------------

fn xorshift64(state: &mut u64) -> u64 {
    *state ^= *state << 13;
    *state ^= *state >> 7;
    *state ^= *state << 17;
    *state
}

// ---------------------------------------------------------------------------
// Reservoir (Algorithm R)
// ---------------------------------------------------------------------------

struct Reservoir {
    samples: Vec<(u64, u64, u64)>, // (exact_k, prime_index, prime_value)
    total:   u64,
    cap:     usize,
    rng:     u64,
}

impl Reservoir {
    fn new(cap: usize, seed: u64) -> Self {
        Self { samples: Vec::with_capacity(cap), total: 0, cap, rng: seed }
    }

    fn push(&mut self, k: u64, prime_index: u64, prime_value: u64) {
        self.total += 1;
        if self.samples.len() < self.cap {
            self.samples.push((k, prime_index, prime_value));
        } else {
            // Replace a random existing slot with probability cap/total
            let j = (xorshift64(&mut self.rng) % self.total) as usize;
            if j < self.cap {
                self.samples[j] = (k, prime_index, prime_value);
            }
        }
    }

    fn finalise(&mut self) {
        // Restore ascending prime_index order (reservoir sampling shuffles this)
        self.samples.sort_unstable_by_key(|&(_, idx, _)| idx);
    }
}

// ---------------------------------------------------------------------------
// Public command
// ---------------------------------------------------------------------------

pub fn cmd_plot_prep(results_path: &Path, samples: usize, split: bool, outdir: &Path) {
    eprintln!("Reading {} ...", results_path.display());

    let file = File::open(results_path).unwrap_or_else(|e| {
        eprintln!("Cannot open {}: {}", results_path.display(), e);
        std::process::exit(1);
    });

    // 8 MB read buffer helps with the large file
    let reader = BufReader::with_capacity(8 * 1024 * 1024, file);

    const MAX_M: usize = 10;
    // Use distinct seeds per class so the random decisions are independent
    let mut res: Vec<Reservoir> = (0..MAX_M)
        .map(|m| Reservoir::new(samples, 0xcafe_babe_dead_beef_u64.wrapping_add(m as u64 * 6_364_136_223_846_793_005)))
        .collect();

    let mut row_count: u64 = 0;
    let mut lines = reader.lines();
    lines.next(); // skip header line

    for line in lines {
        let line = match line {
            Ok(l)  => l,
            Err(e) => { eprintln!("Read error: {}", e); break; }
        };

        // Format: index,value,m_gap  (no spaces, comma-separated)
        let mut parts = line.splitn(3, ',');
        let prime_index: u64 = match parts.next().and_then(|s| s.parse().ok()) {
            Some(v) => v,
            None    => continue,
        };
        let prime_value: u64 = match parts.next().and_then(|s| s.parse().ok()) {
            Some(v) => v,
            None    => continue,
        };
        let m: usize = match parts.next().and_then(|s| s.trim().parse().ok()) {
            Some(v) => v,
            None    => continue,
        };

        if m < MAX_M {
            let k = res[m].total + 1; // exact 1-based rank within this m-class
            res[m].push(k, prime_index, prime_value);
        }

        row_count += 1;
        if row_count % 25_000_000 == 0 {
            eprint!("\r  {:.0}M rows processed ...", row_count as f64 / 1_000_000.0);
        }
    }
    eprintln!("\r  {} rows processed.              ", row_count);

    // Write output TSV(s)
    fs::create_dir_all(outdir).unwrap();

    // Helper to write one class into a BufWriter and return the path used
    fn write_class(
        w: &mut BufWriter<File>,
        m: usize,
        r: &mut Reservoir,
    ) {
        r.finalise();
        if r.samples.is_empty() { return; }

        let class_size = r.total;
        let samples_owned = std::mem::take(&mut r.samples);

        for (k, pi, pv) in samples_owned.iter() {
            writeln!(w, "{}\t{}\t{}\t{}\t{}", m, k, class_size, pi, pv).unwrap();
        }
    }

    if split {
        // One file per m-class
        for (m, r) in res.iter_mut().enumerate() {
            if r.total == 0 { continue; }
            let out_path = outdir.join(format!("plot_data_m{}.tsv", m));
            let f = File::create(&out_path).unwrap();
            let mut w = BufWriter::with_capacity(4 * 1024 * 1024, f);
            writeln!(w, "m\trank\tclass_size\tprime_index\tprime_value").unwrap();
            write_class(&mut w, m, r);
            eprintln!("Wrote {}", out_path.display());
        }
        eprintln!();
        eprintln!("Plot with (e.g.):");
        eprintln!("  python3 scripts/plot.py \"{}\"", outdir.join("plot_data_m2.tsv").display());
    } else {
        // Single combined file
        let out_path = outdir.join("plot_data.tsv");
        let f = File::create(&out_path).unwrap();
        let mut w = BufWriter::with_capacity(4 * 1024 * 1024, f);
        writeln!(w, "m\trank\tclass_size\tprime_index\tprime_value").unwrap();
        for (m, r) in res.iter_mut().enumerate() {
            write_class(&mut w, m, r);
        }
        eprintln!("Wrote {}", out_path.display());
        eprintln!();
        eprintln!("Plot with:");
        eprintln!("  python3 scripts/plot.py \"{}\"", out_path.display());
    }
}
