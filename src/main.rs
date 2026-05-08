use clap::{Parser, Subcommand};
use std::collections::{BTreeMap, HashMap};
use std::fs::{self, File};
use std::io::{self, BufRead, BufWriter, Write};
use std::path::PathBuf;

// ---------------------------------------------------------------------------
// CLI
// ---------------------------------------------------------------------------

#[derive(Parser)]
#[command(name = "pgd", about = "Prime Gap Depth — iterated-regrouping construction on primes")]
struct Cli {
    #[command(subcommand)]
    command: Option<Command>,

    /// Number of primes to use (default 1_000_000)
    #[arg(short = 'n', long, default_value_t = 1_000_000)]
    count: usize,

    /// Supply your own ascending integers (one per line) instead of primes
    #[arg(long, value_name = "FILE")]
    seed_set: Option<PathBuf>,

    /// Output directory for CSV/TSV files
    #[arg(short, long, default_value = "out")]
    outdir: PathBuf,
}

#[derive(Subcommand)]
enum Command {
    /// Verify m(p) is dataset-independent across cutoffs
    Stability,
    /// Distribution of primes mod MOD per m-class with chi-squared p-values
    ModResidue {
        #[arg(default_value_t = 30)]
        modulus: u64,
    },
    /// Show how m-class counts grow with N (100..100M)
    Growth,
    /// Export each small m-class as an OEIS b-file
    OeisExport,
    /// Find the first prime to achieve each m-value up to MAX_M
    FirstAt {
        /// Search up to this m-value (default: 6)
        #[arg(default_value_t = 6)]
        max_m: u32,
    },
    /// Within each m-class, show the 1st, 10th, 100th, ... prime that hits that level
    ClassQuantiles,
    /// Overlay log-log CDFs of m-classes; estimate horizontal shift between consecutive classes
    Overlay,
    /// Fit intercept(m) for converged classes and project forward to higher m
    Predict {
        /// Lowest m to include in the fit (default: 3 — earlier classes are transient)
        #[arg(long, default_value_t = 3)]
        m_min: u32,
        /// Project up to this m
        #[arg(long, default_value_t = 10)]
        m_max: u32,
    },
}

// ---------------------------------------------------------------------------
// Sieve
// ---------------------------------------------------------------------------

/// Simple segmented sieve returning the first `n` primes.
fn sieve_first_n(n: usize) -> Vec<u64> {
    if n == 0 {
        return vec![];
    }
    // Upper bound via prime counting function approximation
    let limit: u64 = if n < 6 {
        15
    } else {
        let nl = n as f64;
        let bound = nl * (nl.ln() + nl.ln().ln() + 2.0);
        bound.ceil() as u64 + 100
    };
    sieve_up_to(limit, Some(n))
}

/// Sieve of Eratosthenes up to `limit`, returning at most `max_count` primes.
/// Uses a segmented approach with 512 KiB segments (fits in L2 cache).
fn sieve_up_to(limit: u64, max_count: Option<usize>) -> Vec<u64> {
    const SEG: usize = 1 << 19; // 512 Ki bits = 512 KiB
    let limit = limit as usize;
    let sqrt_limit = (limit as f64).sqrt() as usize + 1;

    // Small sieve for primes up to sqrt(limit)
    let mut small = vec![true; sqrt_limit + 1];
    small[0] = false;
    if sqrt_limit >= 1 {
        small[1] = false;
    }
    for i in 2..=sqrt_limit {
        if small[i] {
            let mut j = i * i;
            while j <= sqrt_limit {
                small[j] = false;
                j += i;
            }
        }
    }
    let small_primes: Vec<usize> = (2..=sqrt_limit).filter(|&i| small[i]).collect();

    let mut primes: Vec<u64> = Vec::with_capacity(max_count.unwrap_or(50_000_000));
    let cap = max_count.unwrap_or(usize::MAX);

    let mut low = 0usize;
    while low <= limit && primes.len() < cap {
        let high = (low + SEG).min(limit + 1);
        let seg_len = high - low;
        let mut sieve = vec![true; seg_len];

        if low == 0 {
            if seg_len > 0 {
                sieve[0] = false;
            }
            if seg_len > 1 {
                sieve[1] = false;
            }
        }

        for &p in &small_primes {
            let start = if p * p >= low {
                p * p
            } else {
                let rem = low % p;
                if rem == 0 { low } else { low + p - rem }
            };
            let mut j = start;
            while j < high {
                sieve[j - low] = false;
                j += p;
            }
        }

        for i in 0..seg_len {
            if sieve[i] {
                primes.push((low + i) as u64);
                if primes.len() >= cap {
                    break;
                }
            }
        }
        low += SEG;
    }
    primes
}

// ---------------------------------------------------------------------------
// Core construction: compute m(p) for every prime
// ---------------------------------------------------------------------------

/// Returns a Vec<u32> of the same length as `numbers`, where each entry is
/// the m-value (depth level) for that number.
fn compute_m(numbers: &[u64]) -> Vec<u32> {
    let n = numbers.len();
    let mut m_values = vec![u32::MAX; n];

    // Map value -> index for O(1) lookups
    // Work queue: each entry is (level, indices_of_row_members)
    // We use index-based rows to avoid cloning values
    let mut queue: Vec<(u32, Vec<usize>)> = vec![(0, (0..n).collect())];

    while let Some((level, row)) = queue.pop() {
        if row.is_empty() {
            continue;
        }
        // First element of row gets this level
        let first_idx = row[0];
        m_values[first_idx] = level;

        if row.len() == 1 {
            continue;
        }

        // Compute gaps and group remaining by gap value
        // Use a BTreeMap so buckets are processed in consistent order
        let mut buckets: BTreeMap<u64, Vec<usize>> = BTreeMap::new();
        for i in 1..row.len() {
            let gap = numbers[row[i]] - numbers[row[i - 1]];
            buckets.entry(gap).or_default().push(row[i]);
        }

        for (_, bucket) in buckets {
            queue.push((level + 1, bucket));
        }
    }

    m_values
}

// ---------------------------------------------------------------------------
// Analysis helpers
// ---------------------------------------------------------------------------

fn build_histogram(m_values: &[u32]) -> BTreeMap<u32, usize> {
    let mut hist = BTreeMap::new();
    for &m in m_values {
        *hist.entry(m).or_insert(0) += 1;
    }
    hist
}

fn print_histogram(hist: &BTreeMap<u32, usize>) {
    println!("\nm-value histogram:");
    println!("{:<8} {:>12}", "m", "count");
    println!("{}", "-".repeat(22));
    for (&m, &cnt) in hist {
        println!("{:<8} {:>12}", m, cnt);
    }
    let max_m = hist.keys().max().copied().unwrap_or(0);
    println!("\nmax m = {}", max_m);
}

fn print_per_level(numbers: &[u64], m_values: &[u32], hist: &BTreeMap<u32, usize>) {
    println!("\nPer m-level (first 10 primes):");
    println!("{:<6} {:>10}  first 10", "m", "count");
    println!("{}", "-".repeat(60));

    let max_m = hist.keys().max().copied().unwrap_or(0);
    for level in 0..=max_m {
        let count = hist.get(&level).copied().unwrap_or(0);
        let mut first10: Vec<u64> = numbers
            .iter()
            .zip(m_values.iter())
            .filter(|(_, &m)| m == level)
            .map(|(&p, _)| p)
            .take(10)
            .collect();
        first10.sort_unstable();
        let s: Vec<String> = first10.iter().map(|p| p.to_string()).collect();
        println!("{:<6} {:>10}  [{}]", level, count, s.join(", "));
    }
}

fn write_csv(
    outdir: &PathBuf,
    numbers: &[u64],
    m_values: &[u32],
) -> io::Result<()> {
    fs::create_dir_all(outdir)?;
    let path = outdir.join("results.csv");
    let file = File::create(&path)?;
    let mut w = BufWriter::new(file);
    writeln!(w, "pi,p,m")?;
    for (i, (&p, &m)) in numbers.iter().zip(m_values.iter()).enumerate() {
        writeln!(w, "{},{},{}", i + 1, p, m)?;
    }
    eprintln!("Wrote {}", path.display());
    Ok(())
}

// ---------------------------------------------------------------------------
// Subcommands
// ---------------------------------------------------------------------------

fn cmd_stability() {
    let cutoffs = [1_000usize, 10_000, 100_000, 1_000_000];
    println!("Stability check: computing m for first 1000 primes at each cutoff");
    println!("{:<12} {:>12} {:>16}", "cutoff", "unstable", "first 1000 m-vals hash");
    println!("{}", "-".repeat(45));

    let mut baseline: Option<Vec<u32>> = None;
    let mut unstable_total = 0usize;

    for &cutoff in &cutoffs {
        let primes = sieve_first_n(cutoff);
        let m_full = compute_m(&primes);
        let m1000: Vec<u32> = m_full.into_iter().take(1000).collect();
        let unstable = match &baseline {
            None => 0,
            Some(base) => base.iter().zip(m1000.iter()).filter(|(a, b)| a != b).count(),
        };
        unstable_total += unstable;
        // Simple checksum
        let chk: u64 = m1000.iter().enumerate().map(|(i, &m)| (i as u64 + 1) * m as u64).sum();
        println!("{:<12} {:>12} {:>16}", cutoff, unstable, chk);
        baseline = Some(m1000);
    }
    println!("\nTotal unstable primes across cutoffs: {}", unstable_total);
    if unstable_total == 0 {
        println!("PASS: m(p) is dataset-independent for the first 1000 primes.");
    } else {
        println!("FAIL: some m-values changed across cutoffs!");
    }
}

fn cmd_mod_residue(n: usize, seed: Option<&PathBuf>, outdir: &PathBuf, modulus: u64) {
    let numbers = load_numbers(n, seed);
    let m_values = compute_m(&numbers);

    // Baseline: all numbers mod modulus
    let mut baseline_counts: HashMap<u64, usize> = HashMap::new();
    for &p in &numbers {
        *baseline_counts.entry(p % modulus).or_insert(0) += 1;
    }
    let total = numbers.len() as f64;

    let max_m = *m_values.iter().max().unwrap_or(&0);

    // Residues that actually appear among primes
    let mut residues: Vec<u64> = baseline_counts.keys().copied().collect();
    residues.sort_unstable();

    println!("Mod-residue analysis (mod {})", modulus);
    println!("m-level counts and residue distributions vs baseline\n");

    fs::create_dir_all(outdir).unwrap();
    let path = outdir.join(format!("mod_{}.tsv", modulus));
    let file = File::create(&path).unwrap();
    let mut w = BufWriter::new(file);

    // Header
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

        // Print summary to stdout for first few levels
        if level <= 5 {
            println!("m={}: n={:.0}  chi2={:.2}  p≈{:.4}", level, level_n, chi2, p_val);
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

/// Very rough chi-squared p-value via Wilson-Hilferty normal approximation.
fn chi2_p_value_approx(chi2: f64, df: usize) -> f64 {
    if df == 0 || chi2 <= 0.0 {
        return 1.0;
    }
    let k = df as f64;
    // Normal approx: z = ((chi2/k)^(1/3) - (1 - 2/(9k))) / sqrt(2/(9k))
    let cbrt_ratio = (chi2 / k).cbrt();
    let mu = 1.0 - 2.0 / (9.0 * k);
    let sigma = (2.0 / (9.0 * k)).sqrt();
    let z = (cbrt_ratio - mu) / sigma;
    // Upper-tail normal CDF approximation (Abramowitz & Stegun 26.2.17)
    1.0 - normal_cdf(z)
}

fn normal_cdf(x: f64) -> f64 {
    // Hart approximation for Phi(x)
    let t = 1.0 / (1.0 + 0.2316419 * x.abs());
    let poly = t * (0.319381530
        + t * (-0.356563782
            + t * (1.781477937 + t * (-1.821255978 + t * 1.330274429))));
    let pdf = (-x * x / 2.0).exp() / (2.0 * std::f64::consts::PI).sqrt();
    let tail = pdf * poly;
    if x >= 0.0 { 1.0 - tail } else { tail }
}

fn cmd_growth(seed: Option<&PathBuf>, outdir: &PathBuf) {
    let ns: Vec<usize> = vec![100, 1_000, 10_000, 100_000, 1_000_000, 10_000_000];
    // 100M takes too long interactively; include if seed is provided or flag set
    println!("Growth analysis across N values\n");

    // Collect all histograms
    let mut all_hists: Vec<(usize, BTreeMap<u32, usize>)> = Vec::new();
    for &n in &ns {
        let numbers = load_numbers(n, seed);
        let m_values = compute_m(&numbers);
        let hist = build_histogram(&m_values);
        let max_m = hist.keys().max().copied().unwrap_or(0);
        eprint!("N={:>10}  max_m={}  ", n, max_m);
        for (&m, &cnt) in &hist {
            eprint!("m{}={} ", m, cnt);
        }
        eprintln!();
        all_hists.push((n, hist));
    }

    // Find all m-levels observed
    let mut all_levels: Vec<u32> = all_hists
        .iter()
        .flat_map(|(_, h)| h.keys().copied())
        .collect::<std::collections::HashSet<_>>()
        .into_iter()
        .collect::<Vec<_>>();
    all_levels.sort_unstable();

    // Print table
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

    // Write TSV
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

    // Log-log exponents
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
            println!("  m={}: alpha ≈ {:.4}", m, alpha);
        }
    }
}

fn loglog_slope(points: &[(f64, f64)]) -> f64 {
    let n = points.len() as f64;
    let sx: f64 = points.iter().map(|(x, _)| x).sum::<f64>();
    let sy: f64 = points.iter().map(|(_, y)| y).sum::<f64>();
    let sxy: f64 = points.iter().map(|(x, y)| x * y).sum::<f64>();
    let sxx: f64 = points.iter().map(|(x, _)| x * x).sum::<f64>();
    (n * sxy - sx * sy) / (n * sxx - sx * sx)
}

fn cmd_oeis_export(n: usize, seed: Option<&PathBuf>, outdir: &PathBuf) {
    let numbers = load_numbers(n, seed);
    let m_values = compute_m(&numbers);
    let hist = build_histogram(&m_values);

    fs::create_dir_all(outdir).unwrap();

    for (&level, &count) in &hist {
        if count >= 10_000 {
            continue;
        }
        let path = outdir.join(format!("oeis_m{}.txt", level));
        let file = File::create(&path).unwrap();
        let mut w = BufWriter::new(file);
        writeln!(w, "# m-class {} of prime gap depth construction", level).unwrap();
        writeln!(w, "# N = {} primes used", n).unwrap();
        writeln!(w, "# Format: index value (1-based index of the prime)").unwrap();
        let mut idx = 1usize;
        for (i, (&p, &m)) in numbers.iter().zip(m_values.iter()).enumerate() {
            if m == level {
                writeln!(w, "{} {}", idx, p).unwrap();
                let _ = i;
                idx += 1;
            }
        }
        eprintln!("Wrote {}", path.display());
    }
}

// ---------------------------------------------------------------------------
// first-at: find the smallest prime (and its index) for each m-value
// ---------------------------------------------------------------------------

fn cmd_first_at(max_m: u32) {
    // For each target m-level, we need enough primes that the level appears.
    // Strategy: double N until all levels 0..=max_m have been seen, capping at 100M.
    // Since m(p) is dataset-independent, the first prime to achieve level k is stable.

    println!("Searching for first prime at each m-level 0..={}\n", max_m);

    // results[k] = Some((pi, p)) — 1-based prime index and prime value
    let mut results: Vec<Option<(usize, u64)>> = vec![None; max_m as usize + 1];
    let mut n = 1_000usize;
    const N_MAX: usize = 100_000_000;

    loop {
        eprint!("  sieving {} primes... ", n);
        let primes = sieve_first_n(n);
        let m_values = compute_m(&primes);

        // For each level, find the minimum prime that has that m-value.
        // Because the list is already sorted ascending, the first match is the minimum.
        for level in 0..=max_m {
            if results[level as usize].is_some() {
                continue;
            }
            if let Some(pos) = m_values.iter().position(|&m| m == level) {
                // pos is a 0-based index into the sorted primes array.
                // pi = pos + 1 (1-based prime index, i.e. π(p) = pi means p is the pi-th prime).
                results[level as usize] = Some((pos + 1, primes[pos]));
            }
        }

        let found = results.iter().filter(|r| r.is_some()).count();
        eprintln!("found {}/{} levels", found, max_m + 1);

        if results.iter().all(|r| r.is_some()) || n >= N_MAX {
            break;
        }
        n = (n * 10).min(N_MAX);
    }

    println!("{:<6} {:>16} {:>16}", "m", "pi (index)", "p (prime)");
    println!("{}", "-".repeat(42));
    for (level, result) in results.iter().enumerate() {
        match result {
            Some((pi, p)) => println!("{:<6} {:>16} {:>16}", level, pi, p),
            None => println!("{:<6} {:>16} {:>16}", level, "not found", "—"),
        }
    }
}

// ---------------------------------------------------------------------------
// class-quantiles: within each m-class, show ordinal occurrences
// ---------------------------------------------------------------------------

fn cmd_class_quantiles(n: usize, seed: Option<&PathBuf>, outdir: &PathBuf) {
    eprintln!("Loading {} numbers...", n);
    let numbers = load_numbers(n, seed);
    eprintln!("Computing m-values...");
    let m_values = compute_m(&numbers);

    // Group: for each level, the ordered list of (prime_index, prime_value).
    // numbers is sorted ascending, so we just walk in order.
    let max_m = *m_values.iter().max().unwrap_or(&0);
    let mut classes: Vec<Vec<(usize, u64)>> = vec![Vec::new(); max_m as usize + 1];
    for (i, (&p, &m)) in numbers.iter().zip(m_values.iter()).enumerate() {
        classes[m as usize].push((i + 1, p));
    }

    // Ordinals to report: 1, 2, 5, 10, 100, 1k, 10k, 100k, 1M, 10M
    let ordinals: [usize; 10] = [1, 2, 5, 10, 100, 1_000, 10_000, 100_000, 1_000_000, 10_000_000];

    println!("\nFor each m-level, the k-th prime to reach that level.");
    println!("Columns are 'prime-index (prime-value)' — index is π(p), the prime's 1-based rank.\n");

    // Console table
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

    // TSV output: one row per (level, k) with both index and value, easier for plotting.
    fs::create_dir_all(outdir).unwrap();
    let path = outdir.join("class_quantiles.tsv");
    let mut w = BufWriter::new(File::create(&path).unwrap());
    writeln!(w, "m\tsize\tk\tprime_index\tprime_value").unwrap();
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

    // Log-log slope of prime_index vs k within each class.
    // If pi(k) ~ A * k^beta, then a slope of 1 means roughly linear, < 1 means thinning.
    println!("\nLog-log fit  log(prime_index) vs log(k)  within each class:");
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

// ---------------------------------------------------------------------------
// overlay: test whether m-class log-log CDFs are horizontal shifts of each other
// ---------------------------------------------------------------------------

fn cmd_overlay(n: usize, seed: Option<&PathBuf>, outdir: &PathBuf) {
    eprintln!("Loading {} numbers...", n);
    let numbers = load_numbers(n, seed);
    eprintln!("Computing m-values...");
    let m_values = compute_m(&numbers);

    let max_m = *m_values.iter().max().unwrap_or(&0);

    // For each level, build the sorted sequence of prime indices (1-based) in that class.
    // numbers is already sorted ascending, so iterating in order yields ascending indices.
    let mut classes: Vec<Vec<usize>> = vec![Vec::new(); max_m as usize + 1];
    for (i, &m) in m_values.iter().enumerate() {
        classes[m as usize].push(i + 1);
    }

    // For each class with enough samples, build a downsampled log-log curve.
    // Sample at geometric k = 1, 2, 4, 8, ..., len.
    let log_curves: Vec<Vec<(f64, f64)>> = classes
        .iter()
        .map(|cls| {
            if cls.len() < 8 {
                return Vec::new();
            }
            let mut pts = Vec::new();
            let mut k = 1usize;
            while k <= cls.len() {
                pts.push(((k as f64).ln(), (cls[k - 1] as f64).ln()));
                k = (k * 2).max(k + 1);
            }
            // Always include the last point.
            let last = cls.len();
            let last_ln = (last as f64).ln();
            if pts.last().map(|p| p.0).unwrap_or(0.0) < last_ln {
                pts.push((last_ln, (cls[last - 1] as f64).ln()));
            }
            pts
        })
        .collect();

    // Per-class slope/intercept (single linear fit log_pi = slope * log_k + intercept).
    let mut fits: Vec<Option<(f64, f64)>> = Vec::with_capacity(log_curves.len());
    for curve in &log_curves {
        if curve.len() < 3 {
            fits.push(None);
            continue;
        }
        let slope = loglog_slope(curve);
        let mean_x = curve.iter().map(|p| p.0).sum::<f64>() / curve.len() as f64;
        let mean_y = curve.iter().map(|p| p.1).sum::<f64>() / curve.len() as f64;
        let intercept = mean_y - slope * mean_x;
        fits.push(Some((slope, intercept)));
    }

    println!("\nClass log-log fits  log(pi) = slope * log(k) + intercept\n");
    println!("{:<6} {:>10} {:>10} {:>12}", "m", "size", "slope", "intercept");
    println!("{}", "-".repeat(42));
    for (level, fit) in fits.iter().enumerate() {
        let size = classes[level].len();
        match fit {
            Some((s, b)) => println!("{:<6} {:>10} {:>10.4} {:>12.4}", level, size, s, b),
            None => println!("{:<6} {:>10} {:>10} {:>12}", level, size, "—", "—"),
        }
    }

    // Test the "parallel lines" / horizontal-shift hypothesis between consecutive classes.
    // For consecutive (m, m+1), find the horizontal shift dx that best aligns curve(m+1)
    // onto curve(m), i.e. minimizing sum (y_{m+1}(x) - y_m(x - dx))^2 over the overlapping y range.
    //
    // Equivalent and simpler: at each y, find x_m(y) and x_{m+1}(y) and look at the difference.
    // We invert each curve via linear interpolation in y.
    println!("\nHorizontal-shift test between consecutive classes\n");
    println!(
        "If 'parallel lines' holds, dx_y should be ~constant across the overlap range,\n\
         and the residual after subtracting the mean shift should be small.\n"
    );
    println!(
        "{:<10} {:>14} {:>14} {:>14} {:>14}",
        "pair", "mean dx", "stddev dx", "y-overlap", "n samples"
    );
    println!("{}", "-".repeat(70));

    fs::create_dir_all(outdir).unwrap();
    let path = outdir.join("overlay_shifts.tsv");
    let mut w = BufWriter::new(File::create(&path).unwrap());
    writeln!(w, "m_lo\tm_hi\ty\tx_lo\tx_hi\tdx").unwrap();

    for level in 0..(log_curves.len().saturating_sub(1)) {
        let lo = &log_curves[level];
        let hi = &log_curves[level + 1];
        if lo.len() < 3 || hi.len() < 3 {
            continue;
        }

        // Determine y overlap: y is the second coord (log of prime index).
        let y_lo_min = lo.first().unwrap().1;
        let y_lo_max = lo.last().unwrap().1;
        let y_hi_min = hi.first().unwrap().1;
        let y_hi_max = hi.last().unwrap().1;
        let y_start = y_lo_min.max(y_hi_min);
        let y_end = y_lo_max.min(y_hi_max);
        if y_end <= y_start {
            continue;
        }

        // Sample y uniformly across overlap, find x in each curve via linear interpolation
        // on (y, x) — i.e. invert each curve.
        let n_samples = 50usize;
        let mut shifts = Vec::with_capacity(n_samples);
        for i in 0..n_samples {
            let y = y_start + (y_end - y_start) * (i as f64) / ((n_samples - 1) as f64);
            let x_lo = invert_curve_at_y(lo, y);
            let x_hi = invert_curve_at_y(hi, y);
            match (x_lo, x_hi) {
                (Some(xl), Some(xh)) => {
                    let dx = xh - xl;
                    shifts.push((y, xl, xh, dx));
                    writeln!(w, "{}\t{}\t{}\t{}\t{}\t{}", level, level + 1, y, xl, xh, dx)
                        .unwrap();
                }
                _ => {}
            }
        }

        if shifts.is_empty() {
            continue;
        }
        let n = shifts.len() as f64;
        let mean: f64 = shifts.iter().map(|(_, _, _, d)| d).sum::<f64>() / n;
        let var: f64 = shifts.iter().map(|(_, _, _, d)| (d - mean).powi(2)).sum::<f64>() / n;
        let std = var.sqrt();
        let pair = format!("m{}->m{}", level, level + 1);
        let yspan = format!("{:.2}-{:.2}", y_start, y_end);
        println!(
            "{:<10} {:>14.4} {:>14.4} {:>14} {:>14}",
            pair, mean, std, yspan, shifts.len()
        );
    }

    eprintln!("Wrote {}", path.display());

    // Summary interpretation
    println!(
        "\nInterpretation:\n  - If stddev/mean is small (<<1), classes are well-described as horizontal shifts\n    of each other on log-log axes — the 'parallel lines' / pure-thinning hypothesis.\n  - If stddev/mean is large, the shapes differ and the construction is doing more than thinning."
    );
}

/// Given a strictly-monotone-in-y curve as (x, y) points, return the x at which the
/// piecewise-linear interpolation crosses the given y value, if y is in range.
fn invert_curve_at_y(curve: &[(f64, f64)], y: f64) -> Option<f64> {
    if curve.len() < 2 {
        return None;
    }
    // Curve is (x, y) sorted by x ascending. Because pi is monotone in k for ordered classes,
    // y is also non-decreasing. Find the segment containing y.
    if y < curve[0].1 || y > curve[curve.len() - 1].1 {
        return None;
    }
    let mut lo = 0usize;
    let mut hi = curve.len() - 1;
    while hi - lo > 1 {
        let mid = (lo + hi) / 2;
        if curve[mid].1 <= y {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let (x0, y0) = curve[lo];
    let (x1, y1) = curve[hi];
    if (y1 - y0).abs() < 1e-12 {
        return Some(x0);
    }
    let t = (y - y0) / (y1 - y0);
    Some(x0 + t * (x1 - x0))
}

// ---------------------------------------------------------------------------
// predict: fit intercept(m) on converged classes and extrapolate to higher m
// ---------------------------------------------------------------------------

fn cmd_predict(n: usize, seed: Option<&PathBuf>, m_min: u32, m_max: u32) {
    eprintln!("Loading {} numbers...", n);
    let numbers = load_numbers(n, seed);
    eprintln!("Computing m-values...");
    let m_values = compute_m(&numbers);

    let observed_max_m = *m_values.iter().max().unwrap_or(&0);

    // Reuse: per-class first prime index (1-based).
    // Because numbers is sorted ascending, the first occurrence of m == level is
    // the smallest prime in that class — its 0-based position + 1 is its prime index.
    let mut first_idx: Vec<Option<usize>> = vec![None; (observed_max_m as usize) + 1];
    for (i, &m) in m_values.iter().enumerate() {
        let slot = &mut first_idx[m as usize];
        if slot.is_none() {
            *slot = Some(i + 1);
        }
    }

    // log_idx[m] = ln(first prime index at level m), for m in m_min..=observed_max_m.
    let mut data: Vec<(f64, f64)> = Vec::new(); // (m, log_idx)
    println!("\nFirst-occurrence prime indices (observed):\n");
    println!("{:<6} {:>20} {:>14}", "m", "first prime index", "ln(index)");
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

    // Quadratic least-squares fit: y = a*m^2 + b*m + c.
    let (a, b, c) = quadratic_fit(&data);

    // Linear fit for comparison: y = b1*m + c1.
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

    // Show residuals on the fit data
    println!("\nResiduals at observed m (m >= {}):", m_min);
    println!("{:<6} {:>14} {:>14} {:>14}", "m", "ln(observed)", "quad pred", "residual");
    println!("{}", "-".repeat(50));
    for &(m, y) in &data {
        let pred = a * m * m + b * m + c;
        println!("{:<6} {:>14.4} {:>14.4} {:>14.4}", m as u32, y, pred, y - pred);
    }

    // Project forward
    println!("\nProjections for higher m:");
    println!(
        "{:<6} {:>14} {:>26} {:>26}",
        "m", "ln(idx) pred", "predicted prime index", "approx prime value"
    );
    println!("{}", "-".repeat(76));
    for m in (m_min)..=m_max {
        let mf = m as f64;
        let ln_pred = a * mf * mf + b * mf + c;
        // Convert log prime-index back to an integer; skip if it overflows.
        let idx_pred = ln_pred.exp();
        // Approx p ~ idx * ln(idx) (PNT). Beyond ~10^18 the formatting loses meaning.
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
        "\nNotes:\n  - 'predicted prime index' is what the quadratic extrapolation gives for the smallest\n    prime in that m-class. Compare against observed values where available.\n  - 'approx prime value' uses the PNT estimate p ~ n * ln(n) — only an order of magnitude.\n  - Confidence is highest near the fit range (m={}..{}). Each step beyond that compounds error."
        , m_min, observed_max_m
    );
}

/// Least-squares quadratic fit y = a*x^2 + b*x + c using the normal equations.
fn quadratic_fit(data: &[(f64, f64)]) -> (f64, f64, f64) {
    // Build sums for the 3x3 normal-equation matrix.
    let n = data.len() as f64;
    let mut s = [0.0f64; 5]; // sum x^k for k=1..=4
    let mut sy = 0.0f64;
    let mut sxy = 0.0f64;
    let mut sxxy = 0.0f64;
    for &(x, y) in data {
        let x2 = x * x;
        let x3 = x2 * x;
        let x4 = x3 * x;
        s[0] += x;
        s[1] += x2;
        s[2] += x3;
        s[3] += x4;
        sy += y;
        sxy += x * y;
        sxxy += x2 * y;
    }
    // Normal equations:
    // [ sxx4 sxx3 sxx2 ] [a]   [sxxy]
    // [ sxx3 sxx2 sxx  ] [b] = [sxy ]
    // [ sxx2 sxx  n    ] [c]   [sy  ]
    let m = [
        [s[3], s[2], s[1]],
        [s[2], s[1], s[0]],
        [s[1], s[0], n],
    ];
    let v = [sxxy, sxy, sy];
    solve3(&m, &v).unwrap_or((0.0, 0.0, 0.0))
}

/// Solve a 3x3 linear system by Cramer's rule. Returns None if matrix is singular.
fn solve3(m: &[[f64; 3]; 3], v: &[f64; 3]) -> Option<(f64, f64, f64)> {
    let det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    if det.abs() < 1e-18 {
        return None;
    }
    let detx = v[0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (v[1] * m[2][2] - m[1][2] * v[2])
        + m[0][2] * (v[1] * m[2][1] - m[1][1] * v[2]);
    let dety = m[0][0] * (v[1] * m[2][2] - m[1][2] * v[2])
        - v[0] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * v[2] - v[1] * m[2][0]);
    let detz = m[0][0] * (m[1][1] * v[2] - v[1] * m[2][1])
        - m[0][1] * (m[1][0] * v[2] - v[1] * m[2][0])
        + v[0] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    Some((detx / det, dety / det, detz / det))
}

fn linear_fit(data: &[(f64, f64)]) -> (f64, f64) {
    let n = data.len() as f64;
    let sx: f64 = data.iter().map(|p| p.0).sum();
    let sy: f64 = data.iter().map(|p| p.1).sum();
    let sxy: f64 = data.iter().map(|p| p.0 * p.1).sum();
    let sxx: f64 = data.iter().map(|p| p.0 * p.0).sum();
    let slope = (n * sxy - sx * sy) / (n * sxx - sx * sx);
    let intercept = (sy - slope * sx) / n;
    (slope, intercept)
}

// ---------------------------------------------------------------------------
// Loader: either sieve primes or read seed file
// ---------------------------------------------------------------------------

fn load_numbers(n: usize, seed: Option<&PathBuf>) -> Vec<u64> {
    match seed {
        Some(path) => {
            let file = File::open(path).expect("cannot open seed file");
            let reader = io::BufReader::new(file);
            let mut nums: Vec<u64> = reader
                .lines()
                .filter_map(|l| l.ok())
                .filter_map(|l| l.trim().parse::<u64>().ok())
                .take(n)
                .collect();
            nums.sort_unstable();
            nums.dedup();
            nums
        }
        None => sieve_first_n(n),
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_n100_histogram() {
        let primes = sieve_first_n(100);
        let m_values = compute_m(&primes);
        let hist = build_histogram(&m_values);

        assert_eq!(hist.get(&0).copied().unwrap_or(0), 1,  "m=0 count");
        assert_eq!(hist.get(&1).copied().unwrap_or(0), 9,  "m=1 count");
        assert_eq!(hist.get(&2).copied().unwrap_or(0), 42, "m=2 count");
        assert_eq!(hist.get(&3).copied().unwrap_or(0), 42, "m=3 count");
        assert_eq!(hist.get(&4).copied().unwrap_or(0), 6,  "m=4 count");
        assert_eq!(hist.get(&5), None, "no m=5");
        assert_eq!(*hist.keys().max().unwrap(), 4, "max m");
    }

    #[test]
    fn test_n100_m0() {
        let primes = sieve_first_n(100);
        let m_values = compute_m(&primes);
        let m0: Vec<u64> = primes.iter().zip(m_values.iter())
            .filter(|(_, &m)| m == 0).map(|(&p, _)| p).collect();
        assert_eq!(m0, vec![2]);
    }

    #[test]
    fn test_n100_m1() {
        let primes = sieve_first_n(100);
        let m_values = compute_m(&primes);
        let mut m1: Vec<u64> = primes.iter().zip(m_values.iter())
            .filter(|(_, &m)| m == 1).map(|(&p, _)| p).collect();
        m1.sort_unstable();
        assert_eq!(m1, vec![3, 5, 11, 29, 97, 127, 149, 211, 541]);
    }

    #[test]
    fn test_n100_m4() {
        let primes = sieve_first_n(100);
        let m_values = compute_m(&primes);
        let mut m4: Vec<u64> = primes.iter().zip(m_values.iter())
            .filter(|(_, &m)| m == 4).map(|(&p, _)| p).collect();
        m4.sort_unstable();
        assert_eq!(m4, vec![113, 199, 271, 283, 313, 461]);
    }

    #[test]
    fn test_sieve_correctness() {
        let p = sieve_first_n(10);
        assert_eq!(p, vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29]);
    }
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

fn main() {
    let cli = Cli::parse();
    let n = cli.count;
    let seed = cli.seed_set.as_ref();
    let outdir = &cli.outdir;

    match &cli.command {
        None => {
            // Default: run full analysis
            eprintln!("Loading {} numbers...", n);
            let numbers = load_numbers(n, seed);
            eprintln!("Computing m-values...");
            let m_values = compute_m(&numbers);

            let hist = build_histogram(&m_values);
            print_histogram(&hist);
            print_per_level(&numbers, &m_values, &hist);

            write_csv(outdir, &numbers, &m_values).expect("failed writing CSV");
        }
        Some(Command::Stability) => {
            cmd_stability();
        }
        Some(Command::ModResidue { modulus }) => {
            cmd_mod_residue(n, seed, outdir, *modulus);
        }
        Some(Command::Growth) => {
            cmd_growth(seed, outdir);
        }
        Some(Command::OeisExport) => {
            cmd_oeis_export(n, seed, outdir);
        }
        Some(Command::FirstAt { max_m }) => {
            cmd_first_at(*max_m);
        }
        Some(Command::ClassQuantiles) => {
            cmd_class_quantiles(n, seed, outdir);
        }
        Some(Command::Overlay) => {
            cmd_overlay(n, seed, outdir);
        }
        Some(Command::Predict { m_min, m_max }) => {
            cmd_predict(n, seed, *m_min, *m_max);
        }
    }
}
