/// Storage type for depth levels.
///
/// - `u8`: compact (1 byte/prime), max level 255. Good for prime inputs where
///   max observed m ≤ 6 at N=10^8. Use with `--compact-m` / `compute_m::<u8>`.
/// - `u32`: wide (4 bytes/prime), max level ~4 billion. Required for arbitrary
///   sequences (e.g. digits of π) where m can reach thousands.
pub trait MLevel: Copy + Default + Into<u32> + Ord + 'static {
    /// Convert an internal u32 level counter to this type.
    /// Panics on overflow with a clear message.
    fn from_level(v: u32) -> Self;
}

impl MLevel for u8 {
    fn from_level(v: u32) -> Self {
        u8::try_from(v).unwrap_or_else(|_| panic!(
            "m-level {v} exceeded u8::MAX (255). \
             Use --wide-m or compute_m::<u32>() for non-prime sequences."
        ))
    }
}

impl MLevel for u32 {
    fn from_level(v: u32) -> Self { v }
}

/// Compute m(p) for each number in the input array.
///
/// m(p) is the "depth" — the number of recursive gap-regrouping iterations
/// before p becomes the leader (smallest-position element) of its group. See
/// `docs/algorithm.md` for the formal construction.
///
/// Index convention: the formal definition uses 1-indexed rows with
/// destination-indexed gaps `g_R(j) := s_{i_j} − s_{i_{j-1}}` for `2 ≤ j ≤ k`,
/// and `i_j` (the destination) is what gets bucketed. The loop below uses a
/// 0-indexed Rust range `1..row.len()`, which is the same iteration: index
/// `i` here corresponds to math index `j = i+1`, and `row[i]` is `i_j`.
///
/// Gaps are computed as signed `i64` so the algorithm handles non-monotone
/// inputs (rows ordered by *position*, not value). For monotone inputs all
/// gaps are positive and behavior is identical to the original.
///
/// The type parameter `M` controls storage width. Use `u8` for prime inputs
/// (compact, max m ≤ 255) or `u32` for general sequences.
///
/// # Memory layout
///
/// Uses BFS level-by-level with a Compressed Sparse Row (CSR) representation.
/// Only two level buffers exist at any moment; the previous level is freed
/// before the next is built. Member indices are stored as `u32` (half the
/// footprint of `usize` on 64-bit) and are valid for N < 2^32.
pub fn compute_m<M: MLevel>(numbers: &[u64]) -> Vec<M> {
    let n = numbers.len();
    if n == 0 {
        return vec![];
    }
    assert!(
        n <= u32::MAX as usize,
        "N={} exceeds u32 index capacity; reduce the input size",
        n
    );

    let mut m_values: Vec<M> = vec![M::default(); n];

    // CSR layout: cur_members[cur_offsets[r]..cur_offsets[r+1]] holds the
    // member indices of row r at the current BFS level.
    let mut cur_members: Vec<u32> = (0..n as u32).collect();
    let mut cur_offsets: Vec<u32> = vec![0, n as u32];
    let mut level: u32 = 0;  // internal counter; converted to M when stored

    // Per-row scratch buffer reused across rows to avoid repeated allocation.
    let mut scratch: Vec<(i64, usize, u32)> = Vec::new();

    loop {
        let num_rows = cur_offsets.len() - 1;
        // Reserve capacity: next level has at most (active - rows) members
        // because each row loses exactly its leader.
        let next_cap = cur_members.len().saturating_sub(num_rows);
        let mut next_members: Vec<u32> = Vec::with_capacity(next_cap);
        let mut next_offsets: Vec<u32> = vec![0u32];

        for row_idx in 0..num_rows {
            let start = cur_offsets[row_idx] as usize;
            let end   = cur_offsets[row_idx + 1] as usize;
            let row   = &cur_members[start..end];

            // First member is the leader; assign current level.
            m_values[row[0] as usize] = M::from_level(level);

            if row.len() == 1 {
                continue;
            }

            // Build (gap, parent-position, member-idx) triples.
            // The parent-position tiebreaker preserves original row order
            // within each gap bucket, matching the former BTreeMap behaviour.
            scratch.clear();
            for i in 1..row.len() {
                let gap = numbers[row[i] as usize] as i64
                    - numbers[row[i - 1] as usize] as i64;
                scratch.push((gap, i, row[i]));
            }
            scratch.sort_unstable_by_key(|&(g, pos, _)| (g, pos));

            // Append child rows to the next-level CSR buffers.
            let mut i = 0;
            while i < scratch.len() {
                let gap = scratch[i].0;
                while i < scratch.len() && scratch[i].0 == gap {
                    next_members.push(scratch[i].2);
                    i += 1;
                }
                next_offsets.push(next_members.len() as u32);
            }
        }

        // Free the current level entirely before advancing.
        drop(cur_members);
        drop(cur_offsets);

        if next_members.is_empty() {
            break;
        }

        cur_members = next_members;
        cur_offsets = next_offsets;
        level += 1;
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
    fn depth_for_index<M: MLevel>(&self, idx: usize) -> M {
        // First hop: pi(p) = idx + 1 (1-based index in the primes array).
        let mut current = (idx as u64) + 1;
        let mut steps = 1u32;
        // After the first hop, current is in [1, primes.len()]. Subsequent hops use pi_table.
        loop {
            // Is `current` prime?
            if current > self.max_v || !self.is_prime[current as usize] {
                return M::from_level(steps);
            }
            // Recurse: pi(current) = pi_table[current]
            current = self.pi_table[current as usize];
            steps += 1;
        }
    }
}

pub fn compute_pi_chain<M: MLevel>(primes: &[u64]) -> Vec<M> {
    let ctx = PiChainContext::new(primes);
    (0..primes.len()).map(|i| ctx.depth_for_index::<M>(i)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeMap;
    use crate::sieve::sieve_first_n;

    #[test]
    fn test_n100_histogram() {
        let primes = sieve_first_n(100);
        let m_values = compute_m::<u8>(&primes);
        let mut hist: BTreeMap<u8, usize> = BTreeMap::new();
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
        let m_values = compute_m::<u8>(&primes);
        let m0: Vec<u64> = primes.iter().zip(m_values.iter())
            .filter(|(_, &m)| m == 0).map(|(&p, _)| p).collect();
        assert_eq!(m0, vec![2]);
    }

    #[test]
    fn test_n100_m1() {
        let primes = sieve_first_n(100);
        let m_values = compute_m::<u8>(&primes);
        let mut m1: Vec<u64> = primes.iter().zip(m_values.iter())
            .filter(|(_, &m)| m == 1).map(|(&p, _)| p).collect();
        m1.sort_unstable();
        assert_eq!(m1, vec![3, 5, 11, 29, 97, 127, 149, 211, 541]);
    }

    #[test]
    fn test_n100_m4() {
        let primes = sieve_first_n(100);
        let m_values = compute_m::<u8>(&primes);
        let mut m4: Vec<u64> = primes.iter().zip(m_values.iter())
            .filter(|(_, &m)| m == 4).map(|(&p, _)| p).collect();
        m4.sort_unstable();
        assert_eq!(m4, vec![113, 199, 271, 283, 313, 461]);
    }

    #[test]
    fn test_pi_chain_small() {
        let primes = sieve_first_n(10);
        let depths: Vec<u8> = compute_pi_chain(&primes);
        assert_eq!(depths, vec![1u8, 2, 3, 1, 4, 1, 2, 1, 1, 1]);
    }
}

