pub mod sieve;
pub mod depth;
pub mod stats;
pub mod commands;

pub use sieve::{sieve_first_n, sieve_up_to, load_numbers};
pub use depth::{compute_m, compute_pi_chain, MLevel};
pub use stats::{
    build_histogram, print_histogram, print_per_level, write_csv,
    chi2_p_value_approx, normal_cdf, loglog_slope,
    quadratic_fit, linear_fit, solve3, invert_curve_at_y,
};
