use std::collections::HashSet;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// Scan an existing gap_address.csv and report the cumulative number of unique
/// values in each g_* column at every 10^k row-prefix checkpoint.
///
/// The scan is a single O(n) pass: unique sets are accumulated live and their
/// sizes are snapshotted at each milestone. A direct-index column map and
/// `HashSet<i64>` avoid per-row string allocations.
pub fn cmd_gap_address_scan(
    csv_path: &Path,
    pow: u32,
    output: Option<&Path>,
    show_progress: bool,
) {
    let file = File::open(csv_path).unwrap_or_else(|e| {
        eprintln!("error: cannot open {}: {}", csv_path.display(), e);
        std::process::exit(2);
    });
    // 1 MiB read buffer cuts syscall overhead on large files
    let mut reader = BufReader::with_capacity(1 << 20, file);

    // ── Header ───────────────────────────────────────────────────────────────
    let mut header_line = String::new();
    reader
        .read_line(&mut header_line)
        .expect("failed to read CSV header");
    let header = header_line.trim_end_matches(|c| c == '\n' || c == '\r');
    let headers: Vec<&str> = header.split(',').collect();

    // Collect (field_position, column_name) for every g_* column, in order
    let gap_cols: Vec<(usize, &str)> = headers
        .iter()
        .enumerate()
        .filter(|(_, name)| name.starts_with("g_"))
        .map(|(i, name)| (i, *name))
        .collect();

    if gap_cols.is_empty() {
        eprintln!("error: no gap-address columns found (expected headers named g_1, g_2, …)");
        std::process::exit(2);
    }

    let n_gap = gap_cols.len();

    // Build a direct-index lookup so the inner loop does O(1) dispatch:
    //   col_to_seen[field_index] = Some(index into seen[])
    let max_col = gap_cols.iter().map(|(i, _)| *i).max().unwrap();
    let mut col_to_seen: Vec<Option<usize>> = vec![None; max_col + 1];
    for (seen_idx, &(col_pos, _)) in gap_cols.iter().enumerate() {
        col_to_seen[col_pos] = Some(seen_idx);
    }

    // ── Milestones and unique-value accumulators ──────────────────────────────
    let milestones: Vec<u64> = (1..=pow).map(|e| 10u64.pow(e)).collect();
    let mut seen: Vec<HashSet<i64>> = (0..n_gap).map(|_| HashSet::new()).collect();
    let mut row_count: u64 = 0;
    let mut milestone_index: usize = 0;

    // Report buffer (small — one line per milestone)
    let col_header: String = gap_cols
        .iter()
        .map(|(_, name)| format!("{}_unique", name))
        .collect::<Vec<_>>()
        .join(",");
    let mut report_lines: Vec<String> = Vec::with_capacity(pow as usize + 1);
    report_lines.push(format!("power,S,total_unique,{}", col_header));

    if show_progress {
        eprint!("\rrows=0  next=10^1 ");
        let _ = io::stderr().flush();
    }

    // ── Single O(n) pass ─────────────────────────────────────────────────────
    let mut line = String::new();
    loop {
        line.clear();
        match reader.read_line(&mut line) {
            Ok(0) => break,
            Ok(_) => {}
            Err(e) => {
                eprintln!("\nerror reading CSV: {}", e);
                std::process::exit(2);
            }
        }

        row_count += 1;
        let row = line.trim_end_matches(|c: char| c == '\n' || c == '\r');

        // Update unique-value sets; stop walking fields once past the last g_* column
        for (field_idx, val) in row.split(',').enumerate() {
            if field_idx >= col_to_seen.len() {
                break;
            }
            if let Some(seen_idx) = col_to_seen[field_idx] {
                if !val.is_empty() {
                    if let Ok(v) = val.parse::<i64>() {
                        seen[seen_idx].insert(v);
                    }
                }
            }
        }

        // Snapshot unique counts at each 10^k milestone
        while milestone_index < milestones.len() && row_count >= milestones[milestone_index] {
            let exp = milestone_index + 1;
            let s = milestones[milestone_index];
            let col_sizes: Vec<usize> = seen.iter().map(|h| h.len()).collect();
            let total_unique: usize = col_sizes.iter().sum();
            let counts: String = col_sizes.iter().map(|c| c.to_string()).collect::<Vec<_>>().join(",");
            report_lines.push(format!("{},{},{},{}", exp, s, total_unique, counts));
            if show_progress {
                let next_label = if milestone_index + 1 < milestones.len() {
                    format!("10^{}", milestone_index + 2)
                } else {
                    "(done)".to_string()
                };
                eprint!(
                    "\rrows={}  reached 10^{}  next={} ",
                    row_count, exp, next_label
                );
                let _ = io::stderr().flush();
            }
            milestone_index += 1;
        }

        // Periodic progress tick — every 1 M rows
        if show_progress && row_count % 1_000_000 == 0 && milestone_index < milestones.len() {
            eprint!("\rrows={}  next=10^{} ", row_count, milestone_index + 1);
            let _ = io::stderr().flush();
        }
    }

    if show_progress {
        eprintln!("\rrows={}  done                              ", row_count);
    }

    if milestone_index < milestones.len() {
        eprintln!(
            "warning: file ended after {} rows; only reported through 10^{}",
            row_count, milestone_index
        );
    }

    // ── Write report ─────────────────────────────────────────────────────────
    let report = report_lines.join("\n") + "\n";

    match output {
        Some(out_path) => {
            if let Some(parent) = out_path.parent() {
                if !parent.as_os_str().is_empty() {
                    fs::create_dir_all(parent).expect("failed to create output directory");
                }
            }
            let out_file = File::create(out_path).unwrap_or_else(|e| {
                eprintln!("error: cannot create {}: {}", out_path.display(), e);
                std::process::exit(2);
            });
            let mut w = BufWriter::new(out_file);
            w.write_all(report.as_bytes()).expect("failed to write report");
            w.flush().expect("failed to flush report");
            eprintln!("Wrote {}", out_path.display());
        }
        None => print!("{}", report),
    }
}
