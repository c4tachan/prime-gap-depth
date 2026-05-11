/// gen_gap_file: converts a one-prime-per-line ASCII file into a compact binary
/// gap file for fast loading by pgd.
///
/// Binary format:
///   [8 bytes] first prime as u64 little-endian
///   [2 bytes each] predecessor gaps as u16 little-endian (N-1 values for N primes)
///
/// Usage: gen_gap_file <input.txt> <output.gaps>
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() != 3 {
        eprintln!("Usage: gen_gap_file <input.txt> <output.gaps>");
        std::process::exit(1);
    }

    let input_path = &args[1];
    let output_path = &args[2];

    let reader = BufReader::with_capacity(8 * 1024 * 1024, File::open(input_path)
        .unwrap_or_else(|e| { eprintln!("Cannot open {}: {}", input_path, e); std::process::exit(1); }));

    let out_file = File::create(output_path)
        .unwrap_or_else(|e| { eprintln!("Cannot create {}: {}", output_path, e); std::process::exit(1); });
    let mut writer = BufWriter::with_capacity(8 * 1024 * 1024, out_file);

    let mut prev: Option<u64> = None;
    let mut count: u64 = 0;
    let mut max_gap: u64 = 0;

    for (line_no, line) in reader.lines().enumerate() {
        let line = line.unwrap_or_else(|e| { eprintln!("Read error at line {}: {}", line_no + 1, e); std::process::exit(1); });
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        let p: u64 = trimmed.parse().unwrap_or_else(|_| {
            eprintln!("Parse error at line {}: {:?}", line_no + 1, trimmed);
            std::process::exit(1);
        });

        match prev {
            None => {
                // Write header: first prime as u64 LE
                writer.write_all(&p.to_le_bytes())
                    .expect("write error");
            }
            Some(prev_p) => {
                let gap = p.checked_sub(prev_p).unwrap_or_else(|| {
                    eprintln!("Non-monotone input at line {}: {} after {}", line_no + 1, p, prev_p);
                    std::process::exit(1);
                });
                if gap > u16::MAX as u64 {
                    eprintln!("Gap {} at line {} exceeds u16::MAX — cannot encode", gap, line_no + 1);
                    std::process::exit(1);
                }
                if gap > max_gap { max_gap = gap; }
                writer.write_all(&(gap as u16).to_le_bytes())
                    .expect("write error");
            }
        }

        prev = Some(p);
        count += 1;
        if count % 50_000_000 == 0 {
            eprint!("\r  processed {}M primes...", count / 1_000_000);
        }
    }

    writer.flush().expect("flush error");
    eprintln!("\rDone. {} primes written. Max gap: {}. Output: {}", count, max_gap, output_path);
}
