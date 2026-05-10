use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::{self, BufWriter, Write};
use std::path::PathBuf;

pub fn build_histogram(m_values: &[u32]) -> BTreeMap<u32, usize> {
    let mut hist = BTreeMap::new();
    for &m in m_values {
        *hist.entry(m).or_insert(0) += 1;
    }
    hist
}

pub fn print_histogram(hist: &BTreeMap<u32, usize>) {
    println!("\nm-value histogram:");
    println!("{:<8} {:>12}", "m", "count");
    println!("{}", "-".repeat(22));
    for (&m, &cnt) in hist {
        println!("{:<8} {:>12}", m, cnt);
    }
    let max_m = hist.keys().max().copied().unwrap_or(0);
    println!("\nmax m = {}", max_m);
}

pub fn print_per_level(numbers: &[u64], m_values: &[u32], hist: &BTreeMap<u32, usize>) {
    println!("\nPer m-level (first 10 primes):");
    println!("{:<6} {:>10}  first 10", "m", "count");
    println!("{}", "-".repeat(60));

    let max_m = hist.keys().max().copied().unwrap_or(0);
    for level in 0..=max_m {
        let count = hist.get(&level).copied().unwrap_or(0);
        let mut first10: Vec<u64> = numbers
            .iter()
            .zip(m_values.iter())
            .filter(|(_, &m)| m == level)
            .map(|(&p, _)| p)
            .take(10)
            .collect();
        first10.sort_unstable();
        let s: Vec<String> = first10.iter().map(|p| p.to_string()).collect();
        println!("{:<6} {:>10}  [{}]", level, count, s.join(", "));
    }
}

pub fn write_csv(
    outdir: &PathBuf,
    numbers: &[u64],
    m_values: &[u32],
    m_pichain: Option<&[u32]>,
) -> io::Result<()> {
    fs::create_dir_all(outdir)?;
    let path = outdir.join("results.csv");
    let file = File::create(&path)?;
    let mut w = BufWriter::new(file);
    match m_pichain {
        Some(_) => writeln!(w, "pi,p,m_gap,m_pichain")?,
        None => writeln!(w, "pi,p,m_gap")?,
    }
    for (i, (&p, &m)) in numbers.iter().zip(m_values.iter()).enumerate() {
        match m_pichain {
            Some(pc) => writeln!(w, "{},{},{},{}", i + 1, p, m, pc[i])?,
            None => writeln!(w, "{},{},{}", i + 1, p, m)?,
        }
    }
    eprintln!("Wrote {}", path.display());
    Ok(())
}

pub fn chi2_p_value_approx(chi2: f64, df: usize) -> f64 {
    if df == 0 || chi2 <= 0.0 {
        return 1.0;
    }
    let k = df as f64;
    // Wilson-Hilferty normal approximation: z = ((chi2/k)^(1/3) - (1 - 2/(9k))) / sqrt(2/(9k))
    let cbrt_ratio = (chi2 / k).cbrt();
    let mu = 1.0 - 2.0 / (9.0 * k);
    let sigma = (2.0 / (9.0 * k)).sqrt();
    let z = (cbrt_ratio - mu) / sigma;
    1.0 - normal_cdf(z)
}

pub fn normal_cdf(x: f64) -> f64 {
    // Hart approximation for Phi(x)
    let t = 1.0 / (1.0 + 0.2316419 * x.abs());
    let poly = t * (0.319381530
        + t * (-0.356563782
            + t * (1.781477937 + t * (-1.821255978 + t * 1.330274429))));
    let pdf = (-x * x / 2.0).exp() / (2.0 * std::f64::consts::PI).sqrt();
    let tail = pdf * poly;
    if x >= 0.0 { 1.0 - tail } else { tail }
}

pub fn loglog_slope(points: &[(f64, f64)]) -> f64 {
    let n = points.len() as f64;
    let sx: f64 = points.iter().map(|(x, _)| x).sum::<f64>();
    let sy: f64 = points.iter().map(|(_, y)| y).sum::<f64>();
    let sxy: f64 = points.iter().map(|(x, y)| x * y).sum::<f64>();
    let sxx: f64 = points.iter().map(|(x, _)| x * x).sum::<f64>();
    (n * sxy - sx * sy) / (n * sxx - sx * sx)
}

pub fn quadratic_fit(data: &[(f64, f64)]) -> (f64, f64, f64) {
    let n = data.len() as f64;
    let mut s = [0.0f64; 4]; // sum x^k for k=1..=4
    let mut sy = 0.0f64;
    let mut sxy = 0.0f64;
    let mut sxxy = 0.0f64;
    for &(x, y) in data {
        let x2 = x * x;
        let x3 = x2 * x;
        let x4 = x3 * x;
        s[0] += x;
        s[1] += x2;
        s[2] += x3;
        s[3] += x4;
        sy += y;
        sxy += x * y;
        sxxy += x2 * y;
    }
    let m = [
        [s[3], s[2], s[1]],
        [s[2], s[1], s[0]],
        [s[1], s[0], n],
    ];
    let v = [sxxy, sxy, sy];
    solve3(&m, &v).unwrap_or((0.0, 0.0, 0.0))
}

pub fn solve3(m: &[[f64; 3]; 3], v: &[f64; 3]) -> Option<(f64, f64, f64)> {
    let det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    if det.abs() < 1e-18 {
        return None;
    }
    let detx = v[0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (v[1] * m[2][2] - m[1][2] * v[2])
        + m[0][2] * (v[1] * m[2][1] - m[1][1] * v[2]);
    let dety = m[0][0] * (v[1] * m[2][2] - m[1][2] * v[2])
        - v[0] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * v[2] - v[1] * m[2][0]);
    let detz = m[0][0] * (m[1][1] * v[2] - v[1] * m[2][1])
        - m[0][1] * (m[1][0] * v[2] - v[1] * m[2][0])
        + v[0] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    Some((detx / det, dety / det, detz / det))
}

pub fn linear_fit(data: &[(f64, f64)]) -> (f64, f64) {
    let n = data.len() as f64;
    let sx: f64 = data.iter().map(|p| p.0).sum();
    let sy: f64 = data.iter().map(|p| p.1).sum();
    let sxy: f64 = data.iter().map(|p| p.0 * p.1).sum();
    let sxx: f64 = data.iter().map(|p| p.0 * p.0).sum();
    let slope = (n * sxy - sx * sy) / (n * sxx - sx * sx);
    let intercept = (sy - slope * sx) / n;
    (slope, intercept)
}

pub fn invert_curve_at_y(curve: &[(f64, f64)], y: f64) -> Option<f64> {
    for i in 0..curve.len() - 1 {
        let (x0, y0) = curve[i];
        let (x1, y1) = curve[i + 1];
        if (y0 - y).abs() < 1e-10 {
            return Some(x0);
        }
        if (y0 <= y && y <= y1) || (y1 <= y && y <= y0) {
            // Linear interpolation
            if (y1 - y0).abs() < 1e-10 {
                return Some((x0 + x1) / 2.0);
            }
            let t = (y - y0) / (y1 - y0);
            return Some(x0 + t * (x1 - x0));
        }
    }
    None
}
