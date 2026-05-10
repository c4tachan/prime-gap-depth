use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;

/// For each iteration level ℓ, print the rows of `L_ℓ` (and write a CSV).
///
/// Numbers are processed in the order given. Rows are ordered by position
/// (index into `numbers`), not by value, so non-monotone inputs are accepted;
/// in-row gaps are signed integers in that case.
pub fn cmd_iterations(numbers: &[u64], outdir: &Path) {
    if numbers.is_empty() {
        eprintln!("No numbers to process.");
        return;
    }

    let n = numbers.len();
    eprintln!("Computing iteration rows for {} numbers...", n);

    fs::create_dir_all(outdir).expect("failed to create output directory");
    let csv_path = outdir.join("iterations.csv");
    let file = File::create(&csv_path).expect("failed to create iterations.csv");
    let mut csv = BufWriter::new(file);
    writeln!(
        csv,
        "level,row_id,row_size,producing_gap,pos_in_row,index,value,in_row_gap"
    )
    .unwrap();

    // Each entry: (gap that produced this row from its parent — None at level 0, row positions).
    let mut current: Vec<(Option<i64>, Vec<usize>)> = vec![(None, (0..n).collect())];
    let mut level: u32 = 0;

    println!("\n=== Iteration rows ({} numbers) ===", n);

    while !current.is_empty() {
        println!(
            "\n--- level {} ({} row{}) ---",
            level,
            current.len(),
            if current.len() == 1 { "" } else { "s" }
        );

        for (row_id, (producing_gap, row)) in current.iter().enumerate() {
            print_row(numbers, level, row_id, *producing_gap, row);
            write_row_csv(&mut csv, numbers, level, row_id, *producing_gap, row);
        }

        let mut next: Vec<(Option<i64>, Vec<usize>)> = Vec::new();
        for (_, row) in &current {
            if row.len() <= 1 {
                continue;
            }
            let mut buckets: BTreeMap<i64, Vec<usize>> = BTreeMap::new();
            for i in 1..row.len() {
                let gap = numbers[row[i]] as i64 - numbers[row[i - 1]] as i64;
                buckets.entry(gap).or_default().push(row[i]);
            }
            for (g, bucket) in buckets {
                next.push((Some(g), bucket));
            }
        }
        current = next;
        level += 1;
    }

    csv.flush().unwrap();
    eprintln!("Wrote {} ({} levels)", csv_path.display(), level);
}

fn print_row(numbers: &[u64], level: u32, row_id: usize, producing_gap: Option<i64>, row: &[usize]) {
    let header = match producing_gap {
        Some(g) => format!("level {} row {}  gap={}  k={}", level, row_id, g, row.len()),
        None => format!("level {} row {}  k={}", level, row_id, row.len()),
    };
    println!("  {}", header);

    let values: Vec<u64> = row.iter().map(|&i| numbers[i]).collect();
    println!("    values: {}", format_seq_u64(&values));

    if row.len() >= 2 {
        let gaps: Vec<i64> = (1..row.len())
            .map(|i| numbers[row[i]] as i64 - numbers[row[i - 1]] as i64)
            .collect();
        println!("    gaps:   {}", format_seq_i64(&gaps));
    }
}

fn format_seq_u64(xs: &[u64]) -> String {
    format_seq_with(xs, |x| x.to_string())
}

fn format_seq_i64(xs: &[i64]) -> String {
    format_seq_with(xs, |x| x.to_string())
}

fn format_seq_with<T, F: Fn(&T) -> String>(xs: &[T], to_s: F) -> String {
    const MAX_INLINE: usize = 32;
    if xs.len() <= MAX_INLINE {
        let parts: Vec<String> = xs.iter().map(&to_s).collect();
        format!("[{}]", parts.join(", "))
    } else {
        let head: Vec<String> = xs.iter().take(8).map(&to_s).collect();
        let tail_rev: Vec<String> = xs.iter().rev().take(4).map(&to_s).collect();
        let tail: Vec<String> = tail_rev.into_iter().rev().collect();
        format!(
            "[{}, … ({} more) …, {}]",
            head.join(", "),
            xs.len() - 12,
            tail.join(", ")
        )
    }
}

fn write_row_csv<W: Write>(
    csv: &mut W,
    numbers: &[u64],
    level: u32,
    row_id: usize,
    producing_gap: Option<i64>,
    row: &[usize],
) {
    let producing = producing_gap
        .map(|g| g.to_string())
        .unwrap_or_default();
    for (j, &idx) in row.iter().enumerate() {
        let in_row_gap = if j == 0 {
            String::new()
        } else {
            (numbers[idx] as i64 - numbers[row[j - 1]] as i64).to_string()
        };
        writeln!(
            csv,
            "{},{},{},{},{},{},{},{}",
            level,
            row_id,
            row.len(),
            producing,
            j + 1,
            idx + 1,
            numbers[idx],
            in_row_gap
        )
        .unwrap();
    }
}
