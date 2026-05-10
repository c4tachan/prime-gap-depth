use crate::sieve::sieve_first_n;
use crate::depth::compute_m;

pub fn cmd_first_at(max_m: u32) {
    println!("Searching for first prime at each m-level 0..={}\n", max_m);

    let mut results: Vec<Option<(usize, u64)>> = vec![None; max_m as usize + 1];
    let mut n = 1_000usize;
    const N_MAX: usize = 100_000_000;

    loop {
        eprint!("  sieving {} primes... ", n);
        let primes = sieve_first_n(n);
        let m_values = compute_m(&primes);

        for level in 0..=max_m {
            if results[level as usize].is_some() {
                continue;
            }
            if let Some(pos) = m_values.iter().position(|&m| m == level) {
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
