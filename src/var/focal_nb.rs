//! Negative binomial model of per-gene target read counts, for focal CN.
//!
//! Reads on gene g follow NB(mean mu_g, variance mu_g + alpha * mu_g^2) with
//! mu_g = r * (L_g + l): r is reads per base, L_g the gene's target length and
//! l a length offset (a read is counted when any aligned base falls in the
//! target, so reads starting up to about one read length upstream count too).
//! r, l and alpha are fitted by maximum likelihood on the diploid baseline
//! genes; a focal gene's count is then compared with its expected count at 2n.

use std::collections::HashMap;
use std::io::{BufRead, BufReader};

use crate::karyotype::peakfit::ln_gamma;

/// Fewer baseline genes than this and the NB is not fitted.
pub const MIN_NB_BASELINE: usize = 20;
/// Upper bound of the length offset l (bases).
const MAX_OFFSET: f64 = 30_000.0;
const MIN_ALPHA: f64 = 1e-6;
const MAX_ALPHA: f64 = 1.0;

/// Primary reads with >= 1 aligned base in a gene's targets, and the gene's
/// total target length (as in `{prefix}.target_counts.tsv`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GeneReadCount {
    pub reads: u64,
    pub length: u64,
}

/// Fitted NB parameters.
#[derive(Debug, Clone, Copy)]
pub struct NbFit {
    /// Reads per base
    pub r: f64,
    /// Length offset (bases)
    pub l: f64,
    /// Dispersion (variance = mu + alpha mu^2)
    pub alpha: f64,
}

impl NbFit {
    /// Expected reads at 2n for a gene with this target length.
    pub fn mu(&self, length: u64) -> f64 {
        self.r * (length as f64 + self.l)
    }
}

/// NB test of one gene against its expected count at 2n.
#[derive(Debug, Clone, Copy)]
pub struct NbGeneTest {
    pub reads: u64,
    /// Expected reads at 2n
    pub mu: f64,
    /// Two-sided p-value of the count under 2n
    pub p_2n: f64,
}

pub fn test_gene(fit: &NbFit, count: GeneReadCount) -> NbGeneTest {
    let mu = fit.mu(count.length);
    NbGeneTest {
        reads: count.reads,
        mu,
        p_2n: two_sided_p(count.reads, mu, fit.alpha),
    }
}

/// log P(X = n), X ~ NB(mean mu, variance mu + alpha mu^2).
pub fn nb_logpmf(n: u64, mu: f64, alpha: f64) -> f64 {
    let size = 1.0 / alpha;
    let nf = n as f64;
    ln_gamma(nf + size) - ln_gamma(size) - ln_gamma(nf + 1.0)
        - size * (mu / size).ln_1p()
        + nf * (mu / (size + mu)).ln()
}

/// P(X <= n) and P(X >= n). Each tail is summed from n away from the mode,
/// where terms shrink, so the smaller tail keeps its precision; the other is
/// one minus the neighbouring small tail.
fn tails(n: u64, mu: f64, alpha: f64) -> (f64, f64) {
    let size = 1.0 / alpha;
    let q = mu / (size + mu);
    let mode = if size > 1.0 { ((size - 1.0) * q / (1.0 - q)).floor() as u64 } else { 0 };
    // P(X <= k) for k <= mode, summing pmf(k), pmf(k-1), ... via the pmf ratio
    let lower = |k: u64| -> f64 {
        let lp = nb_logpmf(k, mu, alpha);
        let (mut term, mut sum) = (1.0, 1.0);
        let mut j = k;
        while j > 0 {
            // pmf(j-1) / pmf(j) = j / ((j - 1 + size) q)
            term *= j as f64 / ((j as f64 - 1.0 + size) * q);
            sum += term;
            j -= 1;
            if term < 1e-17 * sum {
                break;
            }
        }
        (lp + sum.ln()).exp().min(1.0)
    };
    // P(X >= k) for k >= mode
    let upper = |k: u64| -> f64 {
        let lp = nb_logpmf(k, mu, alpha);
        let (mut term, mut sum) = (1.0, 1.0);
        let mut j = k;
        for _ in 0..50_000_000u64 {
            // pmf(j+1) / pmf(j) = (j + size) q / (j + 1)
            term *= (j as f64 + size) * q / (j as f64 + 1.0);
            sum += term;
            j += 1;
            if term < 1e-17 * sum {
                break;
            }
        }
        (lp + sum.ln()).exp().min(1.0)
    };
    if n <= mode {
        let lo = lower(n);
        let hi = if n == 0 { 1.0 } else { 1.0 - lower(n - 1) };
        (lo, hi)
    } else {
        let hi = upper(n);
        let lo = 1.0 - upper(n + 1);
        (lo, hi)
    }
}

/// Two-sided p-value: min(1, 2 min(P(X <= n), P(X >= n))).
pub fn two_sided_p(n: u64, mu: f64, alpha: f64) -> f64 {
    let (lo, hi) = tails(n, mu, alpha);
    (2.0 * lo.min(hi)).clamp(0.0, 1.0)
}

/// Maps the unconstrained optimizer vector to (r, l, alpha).
fn unpack(x: &[f64; 3]) -> NbFit {
    NbFit {
        r: x[0].exp(),
        l: MAX_OFFSET / (1.0 + (-x[1]).exp()),
        alpha: x[2].clamp(MIN_ALPHA.ln(), MAX_ALPHA.ln()).exp(),
    }
}

/// Maximum-likelihood fit of (r, l, alpha) to the baseline genes' counts.
/// None if there are no reads or the fit is not finite.
pub fn fit_nb(counts: &[GeneReadCount]) -> Option<NbFit> {
    let total_reads: f64 = counts.iter().map(|c| c.reads as f64).sum();
    let total_len: f64 = counts.iter().map(|c| c.length as f64 + 5000.0).sum();
    if counts.is_empty() || total_reads <= 0.0 || total_len <= 0.0 {
        return None;
    }
    let nll = |x: &[f64; 3]| -> f64 {
        let f = unpack(x);
        -counts.iter().map(|c| nb_logpmf(c.reads, f.mu(c.length), f.alpha)).sum::<f64>()
    };
    // l starts at 5 kb, alpha at 0.01
    let x0 = [(total_reads / total_len).ln(), (5000.0f64 / 25000.0).ln(), 0.01f64.ln()];
    let x = nelder_mead(nll, x0, 4000, 1e-6, 1e-6);
    let f = unpack(&x);
    (f.r.is_finite() && f.l.is_finite() && f.alpha.is_finite() && nll(&x).is_finite()).then_some(f)
}

/// Nelder–Mead as in scipy.optimize.minimize(method="Nelder-Mead") with
/// default coefficients, so the fit matches the Python analysis.
fn nelder_mead<F: Fn(&[f64; 3]) -> f64>(f: F, x0: [f64; 3], maxiter: usize, xatol: f64, fatol: f64) -> [f64; 3] {
    const N: usize = 3;
    let (rho, chi, psi, sigma) = (1.0, 2.0, 0.5, 0.5);
    let mut sim = [x0; N + 1];
    for k in 0..N {
        sim[k + 1][k] = if x0[k] != 0.0 { 1.05 * x0[k] } else { 0.00025 };
    }
    let mut fsim: Vec<f64> = sim.iter().map(&f).collect();
    let sort = |sim: &mut [[f64; N]; N + 1], fsim: &mut Vec<f64>| {
        let mut idx: Vec<usize> = (0..=N).collect();
        idx.sort_by(|&a, &b| fsim[a].partial_cmp(&fsim[b]).unwrap_or(std::cmp::Ordering::Equal));
        let (s, fs) = (*sim, fsim.clone());
        for (i, &j) in idx.iter().enumerate() {
            sim[i] = s[j];
            fsim[i] = fs[j];
        }
    };
    sort(&mut sim, &mut fsim);
    let comb = |a: f64, x: &[f64; N], b: f64, y: &[f64; N]| -> [f64; N] {
        std::array::from_fn(|i| a * x[i] + b * y[i])
    };

    for _ in 0..maxiter {
        let xspread = (1..=N).flat_map(|j| (0..N).map(move |i| (j, i)))
            .map(|(j, i)| (sim[j][i] - sim[0][i]).abs())
            .fold(0.0, f64::max);
        let fspread = (1..=N).map(|j| (fsim[0] - fsim[j]).abs()).fold(0.0, f64::max);
        if xspread <= xatol && fspread <= fatol {
            break;
        }
        let xbar: [f64; N] = std::array::from_fn(|i| sim[..N].iter().map(|s| s[i]).sum::<f64>() / N as f64);
        let worst = sim[N];
        let xr = comb(1.0 + rho, &xbar, -rho, &worst);
        let fxr = f(&xr);
        let mut shrink = false;
        if fxr < fsim[0] {
            let xe = comb(1.0 + rho * chi, &xbar, -rho * chi, &worst);
            let fxe = f(&xe);
            if fxe < fxr {
                sim[N] = xe;
                fsim[N] = fxe;
            } else {
                sim[N] = xr;
                fsim[N] = fxr;
            }
        } else if fxr < fsim[N - 1] {
            sim[N] = xr;
            fsim[N] = fxr;
        } else if fxr < fsim[N] {
            let xc = comb(1.0 + psi * rho, &xbar, -psi * rho, &worst);
            let fxc = f(&xc);
            if fxc <= fxr {
                sim[N] = xc;
                fsim[N] = fxc;
            } else {
                shrink = true;
            }
        } else {
            let xcc = comb(1.0 - psi, &xbar, psi, &worst);
            let fxcc = f(&xcc);
            if fxcc < fsim[N] {
                sim[N] = xcc;
                fsim[N] = fxcc;
            } else {
                shrink = true;
            }
        }
        if shrink {
            for j in 1..=N {
                sim[j] = comb(1.0 - sigma, &sim[0], sigma, &sim[j]);
                fsim[j] = f(&sim[j]);
            }
        }
        sort(&mut sim, &mut fsim);
    }
    sim[0]
}

/// Reads `{prefix}.target_counts.tsv` (columns `gene`, `target_len`, `reads`).
pub fn read_target_counts(path: &str) -> Result<HashMap<String, GeneReadCount>, Box<dyn std::error::Error>> {
    let reader = BufReader::new(std::fs::File::open(path)?);
    let mut lines = reader.lines();
    let header = lines.next().ok_or("empty target_counts file")??;
    let cols: Vec<&str> = header.split('\t').collect();
    let col = |name: &str| cols.iter().position(|c| *c == name).ok_or(format!("{} has no {} column", path, name));
    let (gi, li, ri) = (col("gene")?, col("target_len")?, col("reads")?);
    let mut out = HashMap::new();
    for line in lines {
        let line = line?;
        if line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        out.insert(
            f[gi].to_string(),
            GeneReadCount { reads: f[ri].parse()?, length: f[li].parse()? },
        );
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, rel: f64) -> bool {
        (a - b).abs() <= rel * b.abs().max(1e-300)
    }

    /// (n, mu, alpha, scipy logpmf, scipy two-sided p)
    const SCIPY: [(u64, f64, f64, f64, f64); 7] = [
        (0, 5.0, 0.5, -2.505525936990736, 0.16326530612244897),
        (7, 20.0, 0.01, -6.857850477303076, 0.0033263423850953154),
        (150, 100.0, 0.02, -7.402560619541305, 0.00970710442516155),
        (40, 100.0, 0.02, -12.382772468955324, 2.4414650004902656e-05),
        (1000, 800.0, 1e-6, -27.49696788863548, 1.122807960653041e-11),
        (3, 250.0, 0.05, -44.94513913688838, 7.03514333365091e-20),
        (900, 250.0, 0.05, -31.203511700408278, 1.023061578788083e-12),
    ];

    #[test]
    fn logpmf_matches_scipy() {
        for (n, mu, a, lp, _) in SCIPY {
            assert!(close(nb_logpmf(n, mu, a), lp, 1e-9), "n={} mu={} a={}: {}", n, mu, a, nb_logpmf(n, mu, a));
        }
    }

    #[test]
    fn two_sided_p_matches_scipy() {
        for (n, mu, a, _, p) in SCIPY {
            assert!(close(two_sided_p(n, mu, a), p, 1e-6), "n={} mu={} a={}: {}", n, mu, a, two_sided_p(n, mu, a));
        }
        // Count at the mean is not significant
        assert!(two_sided_p(100, 100.0, 0.02) > 0.9);
    }

    #[test]
    fn fit_matches_scipy() {
        let len = [29444, 31685, 45799, 57126, 4021, 10361, 49730, 57021, 16455, 20086, 52403, 26552, 17843,
            50006, 16905, 25733, 39342, 33876, 6972, 3598, 52204, 45703, 50597, 33212, 49417, 21124, 28255,
            47728, 9187, 19585];
        let reads = [152, 135, 254, 211, 47, 44, 167, 257, 107, 99, 269, 138, 102, 163, 94, 114, 120, 157,
            49, 41, 211, 239, 210, 153, 198, 64, 102, 212, 71, 101];
        let counts: Vec<GeneReadCount> =
            len.iter().zip(reads).map(|(&length, reads)| GeneReadCount { reads, length }).collect();
        let f = fit_nb(&counts).unwrap();
        // scipy: r 0.003618001506107363, l 7627.824108045929, alpha 0.02266289862054383
        assert!(close(f.r, 0.003618001506107363, 1e-3), "r {}", f.r);
        assert!(close(f.l, 7627.824108045929, 1e-2), "l {}", f.l);
        assert!(close(f.alpha, 0.02266289862054383, 1e-2), "alpha {}", f.alpha);
    }

    #[test]
    fn fit_needs_reads() {
        assert!(fit_nb(&[]).is_none());
        assert!(fit_nb(&[GeneReadCount { reads: 0, length: 1000 }]).is_none());
    }
}
