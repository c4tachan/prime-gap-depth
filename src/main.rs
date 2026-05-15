use clap::{Parser, Subcommand, ValueEnum};
use std::path::PathBuf;

mod sieve;
mod depth;
mod stats;
mod commands;

use sieve::load_numbers;
use depth::{compute_m, compute_pi_chain, MLevel};
use stats::{build_histogram, print_histogram, print_per_level, write_csv};
use commands::*;

// ---------------------------------------------------------------------------
// CLI
// ---------------------------------------------------------------------------

#[derive(Copy, Clone, Debug, ValueEnum)]
enum Generator {
    Primes,
}

#[derive(Parser)]
#[command(name = "pgd", about = "Prime Gap Depth — iterated-regrouping construction on an input set")]
struct Cli {
    #[command(subcommand)]
    command: Option<Command>,

    /// Number of input elements to use (default 1_000_000)
    #[arg(short = 'n', long, default_value_t = 1_000_000, global = true)]
    count: usize,

    /// Supply your own integer sequence (one per line).
    /// By default the file is sorted+deduped on load (required by the
    /// empirical/statistical commands). Pass --preserve-order to use the
    /// file as-is for non-monotone inputs.
    #[arg(long = "seed-file", alias = "seed-set", value_name = "FILE", global = true, conflicts_with = "generator")]
    seed_file: Option<PathBuf>,

    /// Choose a built-in generator when no seed file is provided.
    #[arg(long, value_enum, global = true, conflicts_with = "seed_file")]
    generator: Option<Generator>,

    /// Output directory for CSV/TSV files
    #[arg(short, long, default_value = "out", global = true)]
    outdir: PathBuf,

    /// Preserve input ordering (skip sort+dedup). Only meaningful with
    /// --seed-file; honored by `iterations`, `gap-address`, and the default
    /// run. The empirical/statistical commands always sort+dedup regardless,
    /// since their analyses assume a monotone sequence.
    #[arg(long, global = true)]
    preserve_order: bool,

    /// Use u32 instead of u8 for m-values. Required for non-prime inputs
    /// where gap-depth m can exceed 255 (e.g. digits of \u03c0 at N\u22652 600).
    /// Default (u8) suffices for prime inputs (max m \u2248 6 at N=10^8) and
    /// uses 4\u00d7 less memory at large N.
    #[arg(long, global = true)]
    wide_m: bool,
}

#[derive(Subcommand)]
enum Command {
    /// Verify m(s) is dataset-independent across cutoffs
    Stability,
    /// Distribution of input values mod MOD per m-class with chi-squared p-values
    ModResidue {
        #[arg(default_value_t = 30)]
        modulus: u64,
    },
    /// Show how m-class counts grow with N (100..100M)
    Growth,
    /// Export each small m-class as an OEIS b-file
    OeisExport {
        /// Read from an existing results.csv instead of recomputing from scratch
        #[arg(long, value_name = "CSV")]
        from_csv: Option<PathBuf>,
    },
    /// Find the first input value to achieve each m-value up to MAX_M
    FirstAt {
        /// Search up to this m-value (default: 6)
        #[arg(default_value_t = 6)]
        max_m: u32,
    },
    /// Within each m-class, show the 1st, 10th, 100th, ... element that hits that level
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
    /// Scan from the last element backwards with a growing window size W. For each
    /// element, compute m on its local predecessor window and check it matches the
    /// global m-value. On mismatch, increment W and restart from the end.
    /// Reports the window size at first acceptance for every element.
    Locality,
    /// Print the rows at each iteration level for a given set of numbers.
    /// Positional NUMBERs are used in the order given (no sorting). Non-monotone
    /// sequences are accepted; in-row gaps become signed integers in that case.
    Iterations {
        /// Optional explicit list of integers, used in the order given.
        /// If omitted, falls back to -n with --seed-file/--generator.
        #[arg(value_name = "NUMBER")]
        numbers: Vec<u64>,
    },
    /// Write a CSV of each number's gap-path address (gap selected at each iteration).
    /// Positional NUMBERs are used in the order given (no sorting). Non-monotone
    /// sequences are accepted; gaps become signed integers in that case.
    GapAddress {
        /// Optional explicit list of integers, used in the order given.
        /// If omitted, falls back to -n with --seed-file/--generator.
        #[arg(value_name = "NUMBER")]
        numbers: Vec<u64>,
    },
    /// Scan an existing gap_address.csv and report the number of unique values per
    /// gap column at each 10^k row-prefix checkpoint (10^1 … 10^pow).
    /// Single O(n) pass — no re-reading the file.
    GapAddressScan {
        /// Path to the gap_address.csv file to scan
        #[arg(value_name = "CSV")]
        csv_path: PathBuf,
        /// Report at prefixes 10^1 through 10^pow
        #[arg(value_name = "POW")]
        pow: u32,
        /// Write the report to this file instead of stdout
        #[arg(short = 'O', long, value_name = "FILE")]
        output: Option<PathBuf>,
        /// Suppress the progress display on stderr
        #[arg(long)]
        no_progress: bool,
    },
    /// Stream an existing results.csv and reservoir-sample each m-class,
    /// writing a compact plot_data.tsv for use with scripts/plot.py.
    /// Does not require --seed-file or --generator.
    PlotPrep {
        /// Path to the results.csv to read
        #[arg(value_name = "CSV")]
        results_path: PathBuf,
        /// Number of reservoir samples per m-class (default: 2000)
        #[arg(short = 's', long, default_value_t = 2000)]
        samples: usize,
        /// Write one plot_data_m{N}.tsv per m-class instead of a single combined file
        #[arg(long)]
        split: bool,
    },
}

// ---------------------------------------------------------------------------
// Default-run helper (generic over m-level storage width)
// ---------------------------------------------------------------------------

fn run_compute<M: MLevel>(
    numbers: &[u64],
    from_generator: bool,
    seed: Option<&PathBuf>,
    outdir: &PathBuf,
) {
    eprintln!("Computing gap-depth m-values...");
    let m_values: Vec<M> = compute_m(numbers);

    let pichain: Option<Vec<M>> = if from_generator && seed.is_none() {
        eprintln!("Computing pi-chain depth...");
        Some(compute_pi_chain(numbers))
    } else {
        None
    };

    println!("\n=== Gap-depth ===");
    let hist = build_histogram(&m_values);
    print_histogram(&hist);
    print_per_level(numbers, &m_values, &hist);

    if let Some(pc) = &pichain {
        println!("\n=== Pi-chain depth ===");
        let pc_hist = build_histogram(pc);
        print_histogram(&pc_hist);
        print_per_level(numbers, pc, &pc_hist);
    }

    write_csv(outdir, numbers, &m_values, pichain.as_deref())
        .expect("failed writing CSV");
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

fn main() {
    let cli = Cli::parse();
    let n = cli.count;
    let seed = cli.seed_file.as_ref();
    let from_generator = cli.generator.is_some();
    let outdir = &cli.outdir;
    let preserve_order = cli.preserve_order;
    let wide_m = cli.wide_m;

    // Commands that read their own files and need no seed/generator
    if let Some(Command::GapAddressScan { csv_path, pow, output, no_progress }) = &cli.command {
        cmd_gap_address_scan(csv_path, *pow, output.as_deref(), !no_progress);
        return;
    }
    if let Some(Command::PlotPrep { results_path, samples, split }) = &cli.command {
        cmd_plot_prep(results_path, *samples, *split, outdir);
        return;
    }
    if let Some(Command::OeisExport { from_csv: Some(csv_path) }) = &cli.command {
        cmd_oeis_export_from_csv(csv_path, outdir);
        return;
    }

    if seed.is_none() && cli.generator.is_none() {
        eprintln!("error: choose an input source: pass --seed-file FILE (or --seed-set FILE) or --generator primes");
        std::process::exit(2);
    }

    match &cli.command {
        None => {
            eprintln!("Loading {} numbers...", n);
            let numbers = load_numbers(n, seed, from_generator, preserve_order);
            if wide_m {
                run_compute::<u32>(&numbers, from_generator, seed, outdir);
            } else {
                run_compute::<u8>(&numbers, from_generator, seed, outdir);
            }
        }
        Some(Command::Stability) => {
            cmd_stability(seed, from_generator);
        }
        Some(Command::ModResidue { modulus }) => {
            cmd_mod_residue(n, seed, from_generator, outdir, *modulus);
        }
        Some(Command::Growth) => {
            cmd_growth(seed, from_generator, outdir);
        }
        Some(Command::OeisExport { from_csv: None }) => {
            cmd_oeis_export(n, seed, from_generator, outdir);
        }
        Some(Command::OeisExport { from_csv: Some(_) }) => unreachable!(),
        Some(Command::FirstAt { max_m }) => {
            cmd_first_at(*max_m, n, seed, from_generator);
        }
        Some(Command::ClassQuantiles) => {
            cmd_class_quantiles(n, seed, from_generator, outdir);
        }
        Some(Command::Overlay) => {
            cmd_overlay(n, seed, from_generator, outdir);
        }
        Some(Command::Predict { m_min, m_max }) => {
            cmd_predict(n, seed, from_generator, *m_min, *m_max);
        }
        Some(Command::PiChain) => {
            cmd_pi_chain(n, seed, from_generator, outdir);
        }
        Some(Command::Iterations { numbers }) => {
            let nums = if numbers.is_empty() {
                load_numbers(n, seed, from_generator, preserve_order)
            } else {
                numbers.clone()
            };
            cmd_iterations(&nums, outdir);
        }
        Some(Command::Locality) => {
            cmd_locality(n, seed, from_generator, outdir);
        }
        Some(Command::GapAddress { numbers }) => {
            let nums = if numbers.is_empty() {
                load_numbers(n, seed, from_generator, preserve_order)
            } else {
                numbers.clone()
            };
            cmd_gap_address(&nums, outdir);
        }
        Some(Command::GapAddressScan { .. }) | Some(Command::PlotPrep { .. }) => {
            // handled above before the seed guard
            unreachable!()
        }
    }
}
