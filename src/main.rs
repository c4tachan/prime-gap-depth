use clap::{Parser, Subcommand};
use std::path::PathBuf;

mod sieve;
mod depth;
mod stats;
mod commands;

use sieve::load_numbers;
use depth::{compute_m, compute_pi_chain};
use stats::{build_histogram, print_histogram, print_per_level, write_csv};
use commands::*;

// ---------------------------------------------------------------------------
// CLI
// ---------------------------------------------------------------------------

#[derive(Parser)]
#[command(name = "pgd", about = "Prime Gap Depth — iterated-regrouping construction on primes")]
struct Cli {
    #[command(subcommand)]
    command: Option<Command>,

    /// Number of primes to use (default 1_000_000)
    #[arg(short = 'n', long, default_value_t = 1_000_000, global = true)]
    count: usize,

    /// Supply your own ascending integers (one per line) instead of primes
    #[arg(long, value_name = "FILE", global = true)]
    seed_set: Option<PathBuf>,

    /// Output directory for CSV/TSV files
    #[arg(short, long, default_value = "out", global = true)]
    outdir: PathBuf,
}

#[derive(Subcommand)]
enum Command {
    /// Verify m(p) is dataset-independent across cutoffs
    Stability,
    /// Distribution of primes mod MOD per m-class with chi-squared p-values
    ModResidue {
        #[arg(default_value_t = 30)]
        modulus: u64,
    },
    /// Show how m-class counts grow with N (100..100M)
    Growth,
    /// Export each small m-class as an OEIS b-file
    OeisExport,
    /// Find the first prime to achieve each m-value up to MAX_M
    FirstAt {
        /// Search up to this m-value (default: 6)
        #[arg(default_value_t = 6)]
        max_m: u32,
    },
    /// Within each m-class, show the 1st, 10th, 100th, ... prime that hits that level
    ClassQuantiles,
    /// Overlay log-log CDFs of m-classes; estimate horizontal shift between consecutive classes
    Overlay,
    /// Fit intercept(m) for converged classes and project forward to higher m
    Predict {
        /// Lowest m to include in the fit (default: 3 — earlier classes are transient)
        #[arg(long, default_value_t = 3)]
        m_min: u32,
        /// Project up to this m
        #[arg(long, default_value_t = 10)]
        m_max: u32,
    },
    /// Pi-chain depth measurement layer: family counts, first appearances, C(m,k), ratios
    PiChain,
    /// Print the rows at each iteration level for a given set of numbers
    Iterations {
        /// Optional explicit list of strictly-ascending positive integers.
        /// If omitted, falls back to -n / --seed-set.
        #[arg(value_name = "NUMBER")]
        numbers: Vec<u64>,
    },
    /// Write a CSV of each number's gap-path address (gap selected at each iteration)
    GapAddress {
        /// Optional explicit list of strictly-ascending positive integers.
        /// If omitted, falls back to -n / --seed-set.
        #[arg(value_name = "NUMBER")]
        numbers: Vec<u64>,
    },
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

fn main() {
    let cli = Cli::parse();
    let n = cli.count;
    let seed = cli.seed_set.as_ref();
    let outdir = &cli.outdir;

    match &cli.command {
        None => {
            eprintln!("Loading {} numbers...", n);
            let numbers = load_numbers(n, seed);
            eprintln!("Computing gap-depth m-values...");
            let m_values = compute_m(&numbers);

            let pichain = if seed.is_none() {
                eprintln!("Computing pi-chain depth...");
                Some(compute_pi_chain(&numbers))
            } else {
                None
            };

            println!("\n=== Gap-depth ===");
            let hist = build_histogram(&m_values);
            print_histogram(&hist);
            print_per_level(&numbers, &m_values, &hist);

            if let Some(pc) = &pichain {
                println!("\n=== Pi-chain depth ===");
                let pc_hist = build_histogram(pc);
                print_histogram(&pc_hist);
                print_per_level(&numbers, pc, &pc_hist);
            }

            write_csv(outdir, &numbers, &m_values, pichain.as_deref())
                .expect("failed writing CSV");
        }
        Some(Command::Stability) => {
            cmd_stability();
        }
        Some(Command::ModResidue { modulus }) => {
            cmd_mod_residue(n, seed, outdir, *modulus);
        }
        Some(Command::Growth) => {
            cmd_growth(seed, outdir);
        }
        Some(Command::OeisExport) => {
            cmd_oeis_export(n, seed, outdir);
        }
        Some(Command::FirstAt { max_m }) => {
            cmd_first_at(*max_m);
        }
        Some(Command::ClassQuantiles) => {
            cmd_class_quantiles(n, seed, outdir);
        }
        Some(Command::Overlay) => {
            cmd_overlay(n, seed, outdir);
        }
        Some(Command::Predict { m_min, m_max }) => {
            cmd_predict(n, seed, *m_min, *m_max);
        }
        Some(Command::PiChain) => {
            cmd_pi_chain(n, outdir);
        }
        Some(Command::Iterations { numbers }) => {
            let nums = if numbers.is_empty() {
                load_numbers(n, seed)
            } else {
                let mut v = numbers.clone();
                v.sort_unstable();
                v.dedup();
                v
            };
            cmd_iterations(&nums, outdir);
        }
        Some(Command::GapAddress { numbers }) => {
            let nums = if numbers.is_empty() {
                load_numbers(n, seed)
            } else {
                let mut v = numbers.clone();
                v.sort_unstable();
                v.dedup();
                v
            };
            cmd_gap_address(&nums, outdir);
        }
    }
}
