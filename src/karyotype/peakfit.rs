//! Beta-mixture MAF peak fitting -- the "k2" method.
//!
//! Port of `maf_peakfit_k2.py` (see `MAF_PEAKFIT_K2_METHOD.txt`). Replaces the
//! histogram local-maximum search in [`super::find_maf_peak`] with a fitted
//! model, so the reported peak is the crest of the largest-AREA bump rather
//! than whatever bin happened to be tallest.
//!
//! For one coverage level's folded MAF (`min(ref,alt)/depth`, in `(0, 0.5]`)
//! the pipeline is:
//!
//! 1. drop MAF below [`LOW_CUT`] (the near-homozygous noise tail),
//! 2. transform to `x = 2*MAF` so the full Beta support `(0,1)` is used and a
//!    balanced diploid peak lands at the `x=1` boundary,
//! 3. fit `p(x) = sum_k w_k Beta(x; a_k,b_k) + wbg Uniform(0,1)` by EM, with a
//!    concentration FLOOR (keeps `a>1` and `b>1`, so the mode stays interior
//!    instead of spiking on the 0.5 wall) and a uniform BACKGROUND (soaks up
//!    the low-MAF error tail),
//! 4. choose K by BIC over `1..=`[`MAX_K`],
//! 5. merge components closer than [`K_SD`] combined SDs and refit, twice,
//! 6. call the peak: crest of the total density inside the winning (largest
//!    weight = largest area) component's basin.
//!
//! Everything here is plain `f64` arithmetic. The one special function the
//! Beta PDF needs is `ln Γ`, implemented below, and since `x` is fixed for the
//! whole fit the per-point `ln x` / `ln(1-x)` are precomputed once -- the EM
//! inner loop is then two multiply-adds and an `exp` per point per component.

// --------------------------------------------------------------------------------
// Method parameters
//
// These are the KEY PARAMETERS block of MAF_PEAKFIT_K2_METHOD.txt. They define
// the method rather than tune it, so they live here instead of in the config;
// promoting one to a threshold later is a two-line change.
// --------------------------------------------------------------------------------

/// Drop folded MAF below this before fitting (the near-homozygous tail).
/// Note `parse_maf` already applies the same cut when reading the file.
pub const LOW_CUT: f64 = 0.10;
/// Largest number of Beta components BIC may choose. 2 is the "k2" method.
const MAX_K: usize = 2;
/// EM restarts; the highest-log-likelihood fit wins.
const N_RESTARTS: usize = 8;
/// Maximum Beta concentration `a+b`.
const C_CAP: f64 = 400.0;
/// Merge adjacent components when `gap < K_SD * (sd_i + sd_j)`.
const K_SD: f64 = 1.2;
/// Merge-then-refit passes.
const MERGE_PASSES: usize = 2;
/// Below this many loci a level's fit is flagged low confidence. Mirrors the
/// `min_maf_sites` threshold that already gates which levels get here at all.
const MIN_N: usize = 50;

/// Clamp applied to `x = 2*MAF` before fitting, keeping it off the `0`/`1`
/// boundary where the Beta PDF diverges.
const FIT_EPS: f64 = 1e-6;
/// The same idea for the readout grid. Deliberately a different (looser) value
/// than [`FIT_EPS`]; both mirror the reference implementation.
const GRID_EPS: f64 = 1e-9;
/// Readout grid: `linspace(0.10, 0.499, 700)` in MAF units.
const GRID_LO: f64 = 0.10;
const GRID_HI: f64 = 0.499;
const GRID_N: usize = 700;
/// EM iteration cap and the log-likelihood improvement below which we stop.
const MAX_ITERS: usize = 300;
const CONVERGE_TOL: f64 = 1e-6;
/// A level with fewer than this many loci after the low cut gets no fit at all
/// (distinct from [`MIN_N`], which only flags low confidence).
const MIN_N_FIT: usize = 5;

// --------------------------------------------------------------------------------
// ln Gamma
// --------------------------------------------------------------------------------

/// Lanczos coefficients, `g = 7`, `n = 9`. Kept at full published precision so
/// they can be checked against the reference by eye, even though `f64` rounds
/// away the last digit or two.
const LANCZOS_G: f64 = 7.0;
#[allow(clippy::excessive_precision)]
const LANCZOS_C: [f64; 9] = [
    0.999_999_999_999_809_93,
    676.520_368_121_885_1,
    -1_259.139_216_722_402_8,
    771.323_428_777_653_13,
    -176.615_029_162_140_59,
    12.507_343_278_686_905,
    -0.138_571_095_265_720_12,
    9.984_369_578_019_571_6e-6,
    1.505_632_735_149_311_6e-7,
];

/// Natural log of the Gamma function, Lanczos approximation (~15 significant
/// digits over the range we use, where `a` and `b` are in `(0, C_CAP]`).
///
/// Rust's `f64::ln_gamma` is still unstable, so we carry our own rather than
/// take a dependency for one function.
pub fn ln_gamma(x: f64) -> f64 {
    if x < 0.5 {
        // Reflection: lnΓ(x) = ln(π / sin(πx)) - lnΓ(1-x).
        std::f64::consts::PI.ln() - (std::f64::consts::PI * x).sin().abs().ln() - ln_gamma(1.0 - x)
    } else {
        let x = x - 1.0;
        let mut acc = LANCZOS_C[0];
        for (i, &c) in LANCZOS_C.iter().enumerate().skip(1) {
            acc += c / (x + i as f64);
        }
        let t = x + LANCZOS_G + 0.5;
        0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + acc.ln()
    }
}

/// `ln B(a,b) = lnΓ(a) + lnΓ(b) - lnΓ(a+b)`, the Beta PDF's normalizer.
fn ln_beta_fn(a: f64, b: f64) -> f64 {
    ln_gamma(a) + ln_gamma(b) - ln_gamma(a + b)
}

/// Beta(a,b) density at a single point. Only used on the readout grid; the EM
/// loop works from precomputed logs instead (see [`MafPrep`]).
fn beta_pdf(x: f64, a: f64, b: f64) -> f64 {
    ((a - 1.0) * x.ln() + (b - 1.0) * (1.0 - x).ln() - ln_beta_fn(a, b)).exp()
}

// --------------------------------------------------------------------------------
// Data prep
// --------------------------------------------------------------------------------

/// One level's folded MAF, transformed to the Beta domain with the per-point
/// logs the E-step needs cached alongside.
struct MafPrep {
    x: Vec<f64>,
    ln_x: Vec<f64>,
    ln_1mx: Vec<f64>,
}

impl MafPrep {
    fn new(maf: &[f64]) -> Self {
        let x: Vec<f64> = maf
            .iter()
            .map(|&m| (2.0 * m).clamp(FIT_EPS, 1.0 - FIT_EPS))
            .collect();
        let ln_x = x.iter().map(|v| v.ln()).collect();
        let ln_1mx = x.iter().map(|v| (1.0 - v).ln()).collect();
        Self { x, ln_x, ln_1mx }
    }

    fn len(&self) -> usize {
        self.x.len()
    }
}

// --------------------------------------------------------------------------------
// Model
// --------------------------------------------------------------------------------

/// A fitted K-component Beta mixture plus its uniform background weight.
#[derive(Debug, Clone)]
pub struct BetaMix {
    /// Per-component mixture weights (also the component's AREA, since each
    /// Beta integrates to 1). `w.sum() + wbg == 1`.
    pub w: Vec<f64>,
    pub a: Vec<f64>,
    pub b: Vec<f64>,
    /// Weight of the `Uniform(0,1)` background -- the estimated noise fraction.
    pub wbg: f64,
    pub log_l: f64,
    pub k: usize,
}

/// One component's readout, in MAF units.
#[derive(Debug, Clone, Copy)]
pub struct Component {
    /// Component mean, converted back to MAF (`mean / 2`).
    pub peak: f64,
    /// Component SD in MAF units.
    pub sd: f64,
    /// Mixture weight = area under this bump.
    pub weight: f64,
    /// `a + b`, i.e. how tight the bump is.
    pub concentration: f64,
}

/// Everything one level's fit produced.
#[derive(Debug, Clone)]
pub struct LevelFit {
    /// Loci surviving the low cut.
    pub n: usize,
    /// `n < MIN_N` -- the fit ran but should not be trusted on its own.
    pub low_confidence: bool,
    /// K chosen by BIC, before merging.
    pub bic_k: usize,
    /// K after the merge-then-refit passes.
    pub final_k: usize,
    /// Background weight of the final fit.
    pub noise_fraction: f64,
    /// Per-component readout, sorted by peak location.
    pub components: Vec<Component>,
    /// The CALLED peak: crest of the largest-area bump. `None` when no fit ran.
    pub overall_peak: Option<f64>,
    pub fit: Option<BetaMix>,
}

// --------------------------------------------------------------------------------
// EM
// --------------------------------------------------------------------------------

/// Deterministic replacement for the reference implementation's random
/// restarts: a Halton (low-discrepancy) point per restart, mapped into the
/// same `(0.2, 0.97)` band of starting means and sorted.
///
/// Halton spreads the restarts more evenly than sampling would, and -- more to
/// the point -- makes nasvar's peak call reproducible run to run.
fn init_means(k: usize, restart: usize) -> Vec<f64> {
    const BASES: [usize; 3] = [2, 3, 5];
    const LO: f64 = 0.2;
    const HI: f64 = 0.97;
    let mut mu: Vec<f64> = (0..k)
        .map(|j| {
            let base = BASES[j % BASES.len()];
            LO + (HI - LO) * van_der_corput(restart + 1, base)
        })
        .collect();
    mu.sort_by(f64::total_cmp);
    mu
}

/// Radical-inverse of `i` in the given base -- the 1-D building block of a
/// Halton sequence. Returns a value in `[0, 1)`.
fn van_der_corput(mut i: usize, base: usize) -> f64 {
    let mut f = 1.0;
    let mut out = 0.0;
    let bf = base as f64;
    while i > 0 {
        f /= bf;
        out += f * (i % base) as f64;
        i /= base;
    }
    out
}

/// Fit a K-component Beta mixture plus a uniform background by EM.
///
/// The M-step is weighted method-of-moments: from the responsibility-weighted
/// mean `m` and variance `v`, the implied concentration is `c = m(1-m)/v - 1`
/// and `(a, b) = (m*c, (1-m)*c)`. No per-component numerical optimization
/// needed, which is what keeps this fast enough to run inside nasvar.
///
/// Two guardrails, both from the method doc:
/// * FLOOR -- raise `c` to at least `1/min(m,1-m) + 2`, which forces `a>1` AND
///   `b>1`. Without it a component collapses into a J-shape spiking at the 0.5
///   boundary, which is an artifact rather than a peak.
/// * BACKGROUND -- a flat `Uniform(0,1)` component absorbs the low-MAF error
///   tail so it can't drag the real peaks down. Its weight is the noise
///   fraction we report.
///
/// EM only finds local optima, so we run [`N_RESTARTS`] initializations and
/// keep the highest log-likelihood.
fn fit_beta_mix(prep: &MafPrep, k: usize) -> BetaMix {
    let n = prep.len();
    let nf = n as f64;
    let n_comp = k + 1; // K Betas + the background
    let mut best: Option<BetaMix> = None;
    // Responsibilities, row-major `[component][point]`, reused across restarts.
    let mut resp = vec![0.0f64; n_comp * n];

    for restart in 0..N_RESTARTS {
        let mu = init_means(k, restart);
        // Start at a moderate concentration of 25, leaving 0.15 for the background.
        let mut a: Vec<f64> = mu.iter().map(|m| m * 25.0).collect();
        let mut b: Vec<f64> = mu.iter().map(|m| (1.0 - m) * 25.0).collect();
        let mut w = vec![0.85 / k as f64; k];
        let mut wbg = 0.15;
        let mut prev = f64::NEG_INFINITY;
        let mut ll = f64::NEG_INFINITY;

        for _ in 0..MAX_ITERS {
            // ---- E-step: unnormalized weight of each component at each point ----
            for kk in 0..k {
                let ln_w = w[kk].ln();
                let ln_norm = ln_beta_fn(a[kk], b[kk]);
                let (am1, bm1) = (a[kk] - 1.0, b[kk] - 1.0);
                let row = &mut resp[kk * n..(kk + 1) * n];
                for ((slot, &lx), &l1x) in row.iter_mut().zip(&prep.ln_x).zip(&prep.ln_1mx) {
                    *slot = (ln_w + am1 * lx + bm1 * l1x - ln_norm).exp();
                }
            }
            // Uniform(0,1) has density 1, so its contribution is just wbg.
            resp[k * n..(k + 1) * n].fill(wbg);

            ll = 0.0;
            for i in 0..n {
                let mut den = 0.0;
                for kk in 0..n_comp {
                    den += resp[kk * n + i];
                }
                den += 1e-300;
                ll += den.ln();
                for kk in 0..n_comp {
                    resp[kk * n + i] /= den;
                }
            }

            // ---- M-step (weights): mean responsibility, renormalized to 1 ----
            for kk in 0..k {
                w[kk] = (resp[kk * n..(kk + 1) * n].iter().sum::<f64>() / nf).max(1e-6);
            }
            wbg = (resp[k * n..(k + 1) * n].iter().sum::<f64>() / nf).max(1e-6);
            let s: f64 = w.iter().sum::<f64>() + wbg;
            for v in w.iter_mut() {
                *v /= s;
            }
            wbg /= s;

            // ---- M-step (shapes): weighted method of moments per component ----
            for kk in 0..k {
                let row = &resp[kk * n..(kk + 1) * n];
                let sk: f64 = row.iter().sum();
                if sk < 1e-6 {
                    continue; // dead component
                }
                let m: f64 = row.iter().zip(&prep.x).map(|(r, x)| r * x).sum::<f64>() / sk;
                let v: f64 = row
                    .iter()
                    .zip(&prep.x)
                    .map(|(r, x)| r * (x - m) * (x - m))
                    .sum::<f64>()
                    / sk;
                // FLOOR: smallest concentration keeping both a>1 and b>1.
                let c = (m * (1.0 - m) / v - 1.0).max(1.0 / m.min(1.0 - m) + 2.0);
                if c > 0.0 {
                    let c = c.min(C_CAP);
                    a[kk] = m * c;
                    b[kk] = (1.0 - m) * c;
                }
            }

            if ll - prev < CONVERGE_TOL {
                break;
            }
            prev = ll;
        }

        // NOTE: `ll` is the log-likelihood of the params as they stood at the
        // TOP of the final iteration, while `w/a/b` are one M-step newer. That
        // is what the reference implementation stores, and restart selection
        // compares these same slightly-stale values, so we keep the behaviour
        // rather than silently changing which restart wins.
        if best.as_ref().is_none_or(|bst| ll > bst.log_l) {
            best = Some(BetaMix {
                w: w.clone(),
                a: a.clone(),
                b: b.clone(),
                wbg,
                log_l: ll,
                k,
            });
        }
    }

    best.expect("N_RESTARTS >= 1 guarantees a fit")
}

// --------------------------------------------------------------------------------
// Model selection
// --------------------------------------------------------------------------------

/// Bayesian Information Criterion, lower is better.
///
/// Free parameters are counted as `3*K`: each component contributes a shape
/// pair `(a,b)` plus a weight, with the background weight and the sum-to-one
/// constraint roughly cancelling.
fn bic(fit: &BetaMix, n: usize) -> f64 {
    3.0 * fit.k as f64 * (n as f64).ln() - 2.0 * fit.log_l
}

/// Fit every `K in 1..=MAX_K` and return the one with the lowest BIC.
fn select_k(prep: &MafPrep) -> BetaMix {
    let n = prep.len();
    (1..=MAX_K)
        .map(|k| fit_beta_mix(prep, k))
        .min_by(|x, y| bic(x, n).total_cmp(&bic(y, n)))
        .expect("MAX_K >= 1 yields at least one fit")
}

// --------------------------------------------------------------------------------
// Readout and merging
// --------------------------------------------------------------------------------

/// Per-component `(peak, sd, weight, concentration)` in MAF units, sorted by
/// peak location.
///
/// For `Beta(a,b)`: `mean = a/(a+b)` and `var = ab/((a+b)^2 (a+b+1))`. Both the
/// location and the SD are halved to undo the `x = 2*MAF` transform.
pub fn components(fit: &BetaMix) -> Vec<Component> {
    let mut out: Vec<Component> = (0..fit.k)
        .map(|k| {
            let (a, b) = (fit.a[k], fit.b[k]);
            let mean = a / (a + b);
            Component {
                peak: 0.5 * mean,
                sd: 0.5 * (mean * (1.0 - mean) / (a + b + 1.0)).sqrt(),
                weight: fit.w[k],
                concentration: a + b,
            }
        })
        .collect();
    out.sort_by(|p, q| p.peak.total_cmp(&q.peak));
    out
}

/// How many peaks survive once components closer than their combined width are
/// treated as one.
///
/// Adjacent components (sorted by peak) are MERGED when
/// `gap = peak_j - peak_i < K_SD * (sd_i + sd_j)`, i.e. they sit within about
/// [`K_SD`] combined SDs and are really describing one bump. This only decides
/// the COUNT; the caller then refits at the reduced K so the shapes are
/// re-estimated cleanly rather than averaged.
///
/// Returns `(surviving_peaks, log)` where each log entry is
/// `(gap, threshold, merged)` for one adjacent pair.
fn n_after_merge(fit: &BetaMix) -> (usize, Vec<(f64, f64, bool)>) {
    let comps = components(fit);
    if comps.len() <= 1 {
        return (comps.len().max(1), Vec::new());
    }
    let mut groups = vec![comps[0]];
    let mut log = Vec::new();
    for next in &comps[1..] {
        let cur = *groups.last().expect("groups seeded with one component");
        let gap = next.peak - cur.peak;
        let thr = K_SD * (cur.sd + next.sd);
        let merged = gap < thr;
        log.push((gap, thr, merged));
        if merged {
            let ws = cur.weight + next.weight;
            *groups.last_mut().expect("groups is non-empty") = Component {
                peak: (cur.peak * cur.weight + next.peak * next.weight) / ws,
                sd: cur.sd.max(next.sd),
                weight: ws,
                concentration: cur.concentration.min(next.concentration),
            };
        } else {
            groups.push(*next);
        }
    }
    (groups.len(), log)
}

/// The readout grid in MAF units: `linspace(GRID_LO, GRID_HI, GRID_N)`.
fn readout_grid() -> Vec<f64> {
    let step = (GRID_HI - GRID_LO) / (GRID_N - 1) as f64;
    (0..GRID_N).map(|i| GRID_LO + step * i as f64).collect()
}

/// Index of the largest element, first one on a tie (matching `argmax`).
fn argmax(xs: &[f64]) -> usize {
    let mut best = 0;
    for (i, v) in xs.iter().enumerate() {
        if *v > xs[best] {
            best = i;
        }
    }
    best
}

/// Total mixture density in MAF space, for plotting or inspection.
///
/// The factor 2 is the Jacobian of `x = 2*MAF`; the background contributes a
/// flat `wbg * 2`.
pub fn density(fit: &BetaMix, grid_maf: &[f64]) -> Vec<f64> {
    grid_maf
        .iter()
        .map(|&g| {
            let x = (2.0 * g).clamp(GRID_EPS, 1.0 - GRID_EPS);
            (0..fit.k)
                .map(|k| fit.w[k] * 2.0 * beta_pdf(x, fit.a[k], fit.b[k]))
                .sum::<f64>()
                + fit.wbg * 2.0
        })
        .collect()
}

/// The CALLED peak: the largest-AREA bump, reported at its VISIBLE crest.
///
/// 1. winner = the highest-weight Beta component (weight is area, since each
///    Beta integrates to 1),
/// 2. basin = the contiguous stretch of grid where that component is the
///    largest contributor,
/// 3. peak = argmax of the TOTAL density inside that basin.
///
/// Using area rather than height is what lets a broad peak near 0.5 win over a
/// taller-but-narrower low-MAF spike; restricting to the basin then keeps the
/// reported crest on the winning bump instead of jumping to a neighbour.
pub fn overall_peak(fit: &BetaMix, grid: &[f64]) -> f64 {
    // Per-component weighted density across the grid, `[component][grid point]`.
    let comp: Vec<Vec<f64>> = (0..fit.k)
        .map(|k| {
            grid.iter()
                .map(|&g| {
                    let x = (2.0 * g).clamp(GRID_EPS, 1.0 - GRID_EPS);
                    fit.w[k] * 2.0 * beta_pdf(x, fit.a[k], fit.b[k])
                })
                .collect()
        })
        .collect();
    let total: Vec<f64> = (0..grid.len())
        .map(|j| comp.iter().map(|c| c[j]).sum::<f64>() + fit.wbg * 2.0)
        .collect();

    if fit.k == 1 {
        return grid[argmax(&total)];
    }

    let kwin = argmax(&fit.w);
    // Where does the winner out-contribute every other component?
    let mask: Vec<bool> = (0..grid.len())
        .map(|j| argmax(&comp.iter().map(|c| c[j]).collect::<Vec<_>>()) == kwin)
        .collect();

    let mean_maf = 0.5 * fit.a[kwin] / (fit.a[kwin] + fit.b[kwin]);
    let mut j0 = 0;
    for (j, &g) in grid.iter().enumerate() {
        if (g - mean_maf).abs() < (grid[j0] - mean_maf).abs() {
            j0 = j;
        }
    }
    if !mask[j0] {
        return mean_maf; // degenerate; fall back to the winner's own mean
    }

    // Grow the contiguous basin around j0.
    let (mut lo, mut hi) = (j0, j0);
    while lo > 0 && mask[lo - 1] {
        lo -= 1;
    }
    while hi < grid.len() - 1 && mask[hi + 1] {
        hi += 1;
    }
    grid[lo + argmax(&total[lo..=hi])]
}

// --------------------------------------------------------------------------------
// Orchestrator
// --------------------------------------------------------------------------------

/// Full pipeline for one level's folded MAF.
///
/// Note the merge loop runs [`MERGE_PASSES`] passes rather than one: the first
/// refit can leave two still-overlapping components that a single pass would
/// never re-check.
pub fn fit_level(maf: &[f64]) -> LevelFit {
    let kept: Vec<f64> = maf.iter().copied().filter(|&m| m >= LOW_CUT).collect();
    let n = kept.len();
    if n < MIN_N_FIT {
        return LevelFit {
            n,
            low_confidence: true,
            bic_k: 0,
            final_k: 0,
            noise_fraction: 0.0,
            components: Vec::new(),
            overall_peak: None,
            fit: None,
        };
    }

    let prep = MafPrep::new(&kept);
    let first = select_k(&prep);
    let bic_k = first.k;

    let mut fit = first;
    for _ in 0..MERGE_PASSES {
        let (k_merge, log) = n_after_merge(&fit);
        for (gap, thr, merged) in &log {
            log::debug!(
                "  peakfit merge check: gap={:.3} thr={:.3} -> {}",
                gap,
                thr,
                if *merged { "merge" } else { "keep" }
            );
        }
        if k_merge < fit.k {
            fit = fit_beta_mix(&prep, k_merge);
        } else {
            break; // stable
        }
    }

    let grid = readout_grid();
    let peak = overall_peak(&fit, &grid);
    LevelFit {
        n,
        low_confidence: n < MIN_N,
        bic_k,
        final_k: fit.k,
        noise_fraction: fit.wbg,
        components: components(&fit),
        overall_peak: Some(peak),
        fit: Some(fit),
    }
}

// --------------------------------------------------------------------------------
// Optional CN mapping (advisory only -- see the method doc's caveats)
// --------------------------------------------------------------------------------

/// Expected folded MAF of a clean CN state at tumor purity `rho`.
///
/// At low purity these all crowd toward 0.5 and stop being distinguishable,
/// which is exactly why the CN call is gated on COVERAGE and not on this.
pub fn expected_folded_maf(cn: u32, rho: f64) -> Option<f64> {
    match cn {
        1 => Some((1.0 - rho) / (2.0 - rho)),
        2 => Some(0.5),
        3 => Some(1.0 / (2.0 + rho)),
        4 => Some(1.0 / (2.0 + 2.0 * rho)),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ln_gamma_matches_known_values() {
        // lnΓ(0.5) = ln(sqrt(pi)); lnΓ(1) = 0; lnΓ(5) = ln(24).
        assert!((ln_gamma(0.5) - 0.572_364_942_924_700_1).abs() < 1e-12);
        assert!(ln_gamma(1.0).abs() < 1e-12);
        assert!((ln_gamma(5.0) - 24.0f64.ln()).abs() < 1e-12);
        assert!((ln_gamma(0.1) - 2.252_712_651_734_206).abs() < 1e-11);
        assert!((ln_gamma(400.0) - 1_994.509_233_436_133_4).abs() < 1e-9);
        assert!((ln_gamma(1.5) + 0.120_782_237_635_245_43).abs() < 1e-12);
    }

    #[test]
    fn beta_pdf_integrates_to_one() {
        // Crude midpoint rule over (0,1) -- enough to catch a wrong normalizer.
        for (a, b) in [(2.0, 5.0), (30.0, 8.0), (1.5, 1.5), (200.0, 20.0)] {
            let n = 200_000;
            let h = 1.0 / n as f64;
            let mass: f64 = (0..n)
                .map(|i| beta_pdf((i as f64 + 0.5) * h, a, b) * h)
                .sum();
            assert!((mass - 1.0).abs() < 1e-4, "Beta({a},{b}) mass = {mass}");
        }
    }

    #[test]
    fn van_der_corput_is_the_expected_sequence() {
        // Base 2: 1/2, 1/4, 3/4, 1/8, 5/8 ...
        let got: Vec<f64> = (1..=5).map(|i| van_der_corput(i, 2)).collect();
        let want = [0.5, 0.25, 0.75, 0.125, 0.625];
        for (g, w) in got.iter().zip(&want) {
            assert!((g - w).abs() < 1e-15);
        }
    }

    #[test]
    fn init_means_are_deterministic_sorted_and_in_range() {
        for restart in 0..N_RESTARTS {
            let mu = init_means(2, restart);
            assert_eq!(mu, init_means(2, restart));
            assert!(mu[0] <= mu[1]);
            assert!(mu.iter().all(|m| (0.2..=0.97).contains(m)));
        }
    }

    /// Build a sample whose histogram follows Beta(a,b), in MAF units. Not a
    /// random draw -- deterministic, so the assertions below are stable.
    fn sample_from_beta(a: f64, b: f64, count: usize) -> Vec<f64> {
        let grid = 2_000;
        let dens: Vec<f64> = (0..grid)
            .map(|i| beta_pdf((i as f64 + 0.5) / grid as f64, a, b))
            .collect();
        let total: f64 = dens.iter().sum();
        let mut out = Vec::with_capacity(count);
        for (i, d) in dens.iter().enumerate() {
            let reps = (d / total * count as f64).round() as usize;
            let x = (i as f64 + 0.5) / grid as f64;
            out.extend(std::iter::repeat_n(x * 0.5, reps)); // -> MAF units
        }
        out
    }

    #[test]
    fn recovers_a_single_balanced_peak() {
        // A tight bump at x=0.9 -> MAF 0.45.
        let maf = sample_from_beta(0.9 * 60.0, 0.1 * 60.0, 20_000);
        let r = fit_level(&maf);
        let peak = r.overall_peak.expect("fit produced a peak");
        assert_eq!(r.final_k, 1, "one bump should collapse to K=1");
        assert!((peak - 0.45).abs() < 0.01, "peak = {peak}");
        assert!(!r.low_confidence);
    }

    #[test]
    fn merges_two_overlapping_components() {
        // Two bumps well inside 1.2 combined SDs of each other.
        let mut maf = sample_from_beta(0.86 * 80.0, 0.14 * 80.0, 10_000);
        maf.extend(sample_from_beta(0.90 * 80.0, 0.10 * 80.0, 10_000));
        let r = fit_level(&maf);
        assert_eq!(r.final_k, 1, "overlapping bumps should merge to one peak");
    }

    #[test]
    fn keeps_two_separated_components() {
        // A low-MAF bump and a near-balanced bump, far apart.
        let mut maf = sample_from_beta(0.45 * 120.0, 0.55 * 120.0, 10_000);
        maf.extend(sample_from_beta(0.94 * 120.0, 0.06 * 120.0, 10_000));
        let r = fit_level(&maf);
        assert_eq!(r.final_k, 2, "separated bumps should stay separate");
        assert_eq!(r.components.len(), 2);
        assert!(r.components[0].peak < r.components[1].peak);
    }

    #[test]
    fn area_beats_height_for_the_called_peak() {
        // A narrow, TALL spike at low MAF plus a broad, larger-AREA bump near
        // 0.5. The call must land on the broad one.
        let mut maf = sample_from_beta(0.30 * 900.0, 0.70 * 900.0, 4_000); // narrow, MAF .15
        maf.extend(sample_from_beta(0.92 * 40.0, 0.08 * 40.0, 16_000)); // broad, MAF ~.46
        let r = fit_level(&maf);
        let peak = r.overall_peak.expect("fit produced a peak");
        assert!(peak > 0.35, "called peak {peak} should be the broad bump");
    }

    #[test]
    fn floor_keeps_the_diploid_peak_off_the_boundary() {
        let maf = sample_from_beta(0.95 * 200.0, 0.05 * 200.0, 20_000);
        let r = fit_level(&maf);
        let peak = r.overall_peak.expect("fit produced a peak");
        assert!(peak < 0.5, "peak {peak} must stay strictly below 0.5");
        assert!(peak > 0.44, "peak {peak} should still read as near-balanced");
    }

    #[test]
    fn too_few_loci_yields_no_fit() {
        let r = fit_level(&[0.42, 0.45, 0.47]);
        assert!(r.fit.is_none());
        assert!(r.overall_peak.is_none());
        assert!(r.low_confidence);
        assert_eq!(r.final_k, 0);
    }

    #[test]
    fn low_cut_drops_the_noise_tail() {
        let mut maf = vec![0.01, 0.02, 0.05, 0.09]; // below LOW_CUT
        maf.extend(sample_from_beta(0.9 * 60.0, 0.1 * 60.0, 1_000));
        let r = fit_level(&maf);
        assert_eq!(r.n, maf.iter().filter(|&&m| m >= LOW_CUT).count());
    }

    #[test]
    fn fit_is_reproducible() {
        let maf = sample_from_beta(0.88 * 50.0, 0.12 * 50.0, 5_000);
        let first = fit_level(&maf);
        let second = fit_level(&maf);
        assert_eq!(first.overall_peak, second.overall_peak);
        assert_eq!(first.final_k, second.final_k);
    }

    #[test]
    fn expected_folded_maf_is_sane() {
        assert_eq!(expected_folded_maf(2, 0.7), Some(0.5));
        // CN3 at full purity -> 1/3; CN1 at full purity -> 0 (pure LOH).
        assert!((expected_folded_maf(3, 1.0).unwrap() - 1.0 / 3.0).abs() < 1e-12);
        assert!(expected_folded_maf(1, 1.0).unwrap().abs() < 1e-12);
        assert_eq!(expected_folded_maf(5, 0.5), None);
    }
}
