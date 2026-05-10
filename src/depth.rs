use std::collections::BTreeMap;

/// Compute m(p) for each number in the input array.
///
/// m(p) is the "depth" — the number of recursive gap-regrouping iterations
/// before p becomes the leader (smallest element) of its group.
///
/// Returns a Vec<u32> of the same length as `numbers`, where each entry is
/// the m-value (depth level) for that number.
pub fn compute_m(numbers: &[u64]) -> Vec<u32> {
    let n = numbers.len();
    let mut m_values = vec![u32::MAX; n];

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

/// Pi-chain depth: m(p) = number of recursive pi-applications before chain ends.
///
/// For each prime p, compute how many times we can apply the prime-counting
/// function pi(p) before reaching a non-prime. E.g.:
/// - p=2: pi(2)=1 not prime -> 1 step
/// - p=11: pi(11)=5, pi(5)=3, pi(3)=2, pi(2)=1 -> 4 steps
struct PiChainContext {
    is_prime: Vec<bool>, // indexed by value, size = max_v + 1
    pi_table: Vec<u64>,  // pi_table[v] = #primes <= v
    max_v: u64,
}

impl PiChainContext {
    fn new(primes: &[u64]) -> Self {
        // pi values can never exceed primes.len() after the first hop; on the first
        // hop they equal the 1-based index of the prime, which is <= primes.len().
        // So max_v = primes.len() is sufficient for all chain values after p itself.
        let max_v = primes.len() as u64;
        let mut is_prime = vec![false; (max_v + 1) as usize];
        // primes is sorted ascending; mark each one that fits in [0, max_v]
        for &p in primes {
            if p <= max_v {
                is_prime[p as usize] = true;
            } else {
                break;
            }
        }
        // pi_table: prefix count
        let mut pi_table = vec![0u64; (max_v + 1) as usize];
        let mut count = 0u64;
        for v in 0..=max_v as usize {
            if is_prime[v] {
                count += 1;
            }
            pi_table[v] = count;
        }
        PiChainContext { is_prime, pi_table, max_v }
    }

    /// Compute pi-chain depth for the prime at 0-based index `idx` in the sorted primes.
    fn depth_for_index(&self, idx: usize) -> u32 {
        // First hop: pi(p) = idx + 1 (1-based index in the primes array).
        let mut current = (idx as u64) + 1;
        let mut steps = 1u32;
        // After the first hop, current is in [1, primes.len()]. Subsequent hops use pi_table.
        loop {
            // Is `current` prime?
            if current > self.max_v || !self.is_prime[current as usize] {
                return steps;
            }
            // Recurse: pi(current) = pi_table[current]
            current = self.pi_table[current as usize];
            steps += 1;
        }
    }
}

pub fn compute_pi_chain(primes: &[u64]) -> Vec<u32> {
    let ctx = PiChainContext::new(primes);
    (0..primes.len()).map(|i| ctx.depth_for_index(i)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sieve::sieve_first_n;

    #[test]
    fn test_n100_histogram() {
        let primes = sieve_first_n(100);
        let m_values = compute_m(&primes);
        let mut hist: BTreeMap<u32, usize> = BTreeMap::new();
        for &m in &m_values {
            *hist.entry(m).or_insert(0) += 1;
        }

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
    fn test_pi_chain_small() {
        let primes = sieve_first_n(10);
        let depths = compute_pi_chain(&primes);
        assert_eq!(depths, vec![1, 2, 3, 1, 4, 1, 2, 1, 1, 1]);
    }
}

