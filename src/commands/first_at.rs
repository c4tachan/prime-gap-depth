use std::path::PathBuf;

use crate::sieve::load_numbers;
use crate::depth::compute_m;

pub fn cmd_first_at(max_m: u32, n: usize, seed: Option<&PathBuf>, use_primes: bool) {
    println!("Searching for first prime at each m-level 0..={}\n", max_m);

    let mut results: Vec<Option<(usize, u64)>> = vec![None; max_m as usize + 1];
    let mut batch_n = 1_000usize;
    const N_MAX: usize = 100_000_000;
    let fixed_input = seed.is_some() || !use_primes;

    loop {
        let current_n = if fixed_input { n } else { batch_n };
        eprint!("  loading {} numbers... ", current_n);
        let numbers = load_numbers(current_n, seed, use_primes, false);
        let m_values = compute_m(&numbers);

        for level in 0..=max_m {
            if results[level as usize].is_some() {
                continue;
            }
            if let Some(pos) = m_values.iter().position(|&m| m == level) {
                results[level as usize] = Some((pos + 1, numbers[pos]));
            }
        }

        let found = results.iter().filter(|r| r.is_some()).count();
        eprintln!("found {}/{} levels", found, max_m + 1);

        if results.iter().all(|r| r.is_some()) || fixed_input || batch_n >= N_MAX {
            break;
        }
        batch_n = (batch_n * 10).min(N_MAX);
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
