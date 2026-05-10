use crate::sieve::sieve_first_n;
use crate::depth::compute_m;

pub fn cmd_stability() {
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
