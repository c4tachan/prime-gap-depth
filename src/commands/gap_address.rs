use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;

/// Write a CSV recording the gap value selecting each number's bucket at every
/// iteration — i.e. the gap-path address `addr(s_i) = (g_1, …, g_m)` from
/// `docs/algorithm.md` §4.
///
/// One row per input number. Columns: `index, value, m, g_1, g_2, …, g_M`,
/// where `M = max m`. Cells beyond a number's leader-level are empty.
///
/// Numbers are processed in the order given. For monotone inputs all gaps are
/// positive; for non-monotone inputs gaps are signed.
pub fn cmd_gap_address(numbers: &[u64], outdir: &Path) {
    if numbers.is_empty() {
        eprintln!("No numbers to process.");
        return;
    }

    let n = numbers.len();
    eprintln!("Computing gap-path addresses for {} numbers...", n);

    let mut address: Vec<Vec<i64>> = vec![Vec::new(); n];
    let mut m_values: Vec<u32> = vec![0u32; n];

    let mut current: Vec<Vec<usize>> = vec![(0..n).collect()];
    let mut level: u32 = 0;

    while !current.is_empty() {
        let mut next: Vec<Vec<usize>> = Vec::new();
        for row in &current {
            m_values[row[0]] = level;
            if row.len() <= 1 {
                continue;
            }
            let mut buckets: BTreeMap<i64, Vec<usize>> = BTreeMap::new();
            for i in 1..row.len() {
                let gap = numbers[row[i]] as i64 - numbers[row[i - 1]] as i64;
                buckets.entry(gap).or_default().push(row[i]);
            }
            for (g, bucket) in buckets {
                for &p in &bucket {
                    address[p].push(g);
                }
                next.push(bucket);
            }
        }
        current = next;
        level += 1;
    }

    let max_m = *m_values.iter().max().unwrap_or(&0);

    fs::create_dir_all(outdir).expect("failed to create output directory");
    let csv_path = outdir.join("gap_address.csv");
    let file = File::create(&csv_path).expect("failed to create gap_address.csv");
    let mut csv = BufWriter::new(file);

    write!(csv, "index,value,m").unwrap();
    for j in 1..=max_m {
        write!(csv, ",g_{}", j).unwrap();
    }
    writeln!(csv).unwrap();

    for i in 0..n {
        write!(csv, "{},{},{}", i + 1, numbers[i], m_values[i]).unwrap();
        for j in 0..max_m as usize {
            if j < address[i].len() {
                write!(csv, ",{}", address[i][j]).unwrap();
            } else {
                write!(csv, ",").unwrap();
            }
        }
        writeln!(csv).unwrap();
    }

    csv.flush().unwrap();
    eprintln!("Wrote {} (max m = {})", csv_path.display(), max_m);
}
