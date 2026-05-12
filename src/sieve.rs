use std::fs::File;
use std::io::{self, BufRead, Read};
use std::path::PathBuf;

/// Returns the first `n` primes using Sieve of Eratosthenes.
pub fn sieve_first_n(n: usize) -> Vec<u64> {
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
pub fn sieve_up_to(limit: u64, max_count: Option<usize>) -> Vec<u64> {
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

/// Load numbers from a seed file or from the prime generator.
///
/// When `preserve_order` is false (the default for callers), seed-file inputs
/// are sorted ascending and de-duplicated — this is what the empirical/
/// statistical commands expect. When true, the file is read as-is, allowing
/// non-monotone or duplicate-bearing sequences to flow through unchanged.
/// Generated primes are always monotone regardless.
/// Load a binary `.gaps` file produced by `gen_gap_file`.
///
/// Format: 8-byte u64 LE first prime, then (N-1) u16 LE gaps.
/// Reconstructs the full prime sequence up to `n` entries.
fn load_gap_file(path: &PathBuf, n: usize) -> Vec<u64> {
    let mut file = File::open(path).expect("cannot open .gaps file");
    let mut bytes = Vec::new();
    file.read_to_end(&mut bytes).expect("cannot read .gaps file");

    if bytes.len() < 8 {
        panic!("gap file too short — missing u64 header");
    }

    let first = u64::from_le_bytes(bytes[0..8].try_into().unwrap());
    let gap_bytes = &bytes[8..];

    if gap_bytes.len() % 2 != 0 {
        panic!("gap file has odd number of gap bytes");
    }

    let num_gaps = gap_bytes.len() / 2;
    let capacity = (num_gaps + 1).min(n);
    let mut primes = Vec::with_capacity(capacity);
    primes.push(first);

    let mut current = first;
    for chunk in gap_bytes.chunks_exact(2) {
        if primes.len() >= n {
            break;
        }
        let gap = u16::from_le_bytes([chunk[0], chunk[1]]) as u64;
        current += gap;
        primes.push(current);
    }

    primes
}

pub fn load_numbers(
    n: usize,
    seed: Option<&PathBuf>,
    from_generator: bool,
    preserve_order: bool,
) -> Vec<u64> {
    match seed {
        Some(path) => {
            if path.extension().and_then(|e| e.to_str()) == Some("gaps") {
                return load_gap_file(path, n);
            }
            let file = File::open(path).expect("cannot open seed file");
            let reader = io::BufReader::new(file);
            let mut nums: Vec<u64> = reader
                .lines()
                .filter_map(|l| l.ok())
                .filter_map(|l| l.trim().parse::<u64>().ok())
                .take(n)
                .collect();
            if !preserve_order {
                nums.sort_unstable();
                nums.dedup();
            }
            nums
        }
        None if from_generator => sieve_first_n(n),
        None => panic!("no input source selected: pass --seed-file FILE or --generator primes"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sieve_correctness() {
        let p = sieve_first_n(10);
        assert_eq!(p, vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29]);
    }
}

