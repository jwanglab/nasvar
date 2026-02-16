//! Affine-gap Needleman-Wunsch pairwise sequence alignment.

/// Alignment scoring parameters.
pub struct AlignScoring {
    pub match_score: i32,
    pub mismatch_score: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
}

impl Default for AlignScoring {
    fn default() -> Self {
        Self {
            match_score: 2,
            mismatch_score: -2,
            gap_open: -3,
            gap_extend: -1,
        }
    }
}

/// A single alignment operation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignOp {
    /// Bases aligned (match or mismatch) — consumes 1 ref + 1 query
    Match,
    /// Gap in query (deletion from query perspective) — consumes 1 ref
    InsRef,
    /// Gap in reference (insertion in query) — consumes 1 query
    InsQry,
}

/// Result of pairwise alignment.
pub struct PairwiseAlignment {
    pub score: i32,
    /// Run-length encoded alignment operations.
    pub ops: Vec<(AlignOp, usize)>,
}

const NEG_INF: i32 = i32::MIN / 2;

/// Traceback direction for the DP.
#[derive(Clone, Copy)]
enum Trace {
    M,
    Ix,
    Iy,
}

/// Align query to reference using affine-gap Needleman-Wunsch.
pub fn needleman_wunsch(
    reference: &[u8],
    query: &[u8],
    scoring: &AlignScoring,
) -> PairwiseAlignment {
    let n = reference.len();
    let m = query.len();
    let cols = m + 1;

    // Flat matrices: index = i * cols + j
    let size = (n + 1) * cols;
    let mut sm = vec![NEG_INF; size]; // Match/mismatch state
    let mut sx = vec![NEG_INF; size]; // Gap in query (consuming ref)
    let mut sy = vec![NEG_INF; size]; // Gap in ref (consuming query)

    // Traceback: which state produced the best score at each cell
    let mut tm = vec![Trace::M; size];
    let mut tx = vec![Trace::M; size];
    let mut ty = vec![Trace::M; size];

    // Initialization
    sm[0] = 0;
    for i in 1..=n {
        sx[i * cols] = scoring.gap_open + (i as i32) * scoring.gap_extend;
    }
    for (j, sy_val) in sy.iter_mut().enumerate().take(m + 1).skip(1) {
        *sy_val = scoring.gap_open + (j as i32) * scoring.gap_extend;
    }

    // Fill
    for i in 1..=n {
        for j in 1..=m {
            let idx = i * cols + j;
            let diag = (i - 1) * cols + (j - 1);
            let up = (i - 1) * cols + j;
            let left = i * cols + (j - 1);

            // Match/mismatch state
            let s = if reference[i - 1].eq_ignore_ascii_case(&query[j - 1]) {
                scoring.match_score
            } else {
                scoring.mismatch_score
            };

            let m_from_m = sm[diag];
            let m_from_x = sx[diag];
            let m_from_y = sy[diag];
            let best_m = m_from_m.max(m_from_x).max(m_from_y);
            sm[idx] = best_m + s;
            tm[idx] = if best_m == m_from_m {
                Trace::M
            } else if best_m == m_from_x {
                Trace::Ix
            } else {
                Trace::Iy
            };

            // Gap in query (InsRef) — consumes reference[i]
            let x_open = sm[up] + scoring.gap_open + scoring.gap_extend;
            let x_ext = sx[up] + scoring.gap_extend;
            if x_open >= x_ext {
                sx[idx] = x_open;
                tx[idx] = Trace::M;
            } else {
                sx[idx] = x_ext;
                tx[idx] = Trace::Ix;
            }

            // Gap in reference (InsQry) — consumes query[j]
            let y_open = sm[left] + scoring.gap_open + scoring.gap_extend;
            let y_ext = sy[left] + scoring.gap_extend;
            if y_open >= y_ext {
                sy[idx] = y_open;
                ty[idx] = Trace::M;
            } else {
                sy[idx] = y_ext;
                ty[idx] = Trace::Iy;
            }
        }
    }

    // Find best final state
    let end_idx = n * cols + m;
    let end_m = sm[end_idx];
    let end_x = sx[end_idx];
    let end_y = sy[end_idx];
    let best_end = end_m.max(end_x).max(end_y);
    let score = best_end;

    let mut state = if best_end == end_m {
        Trace::M
    } else if best_end == end_x {
        Trace::Ix
    } else {
        Trace::Iy
    };

    // Traceback
    let mut raw_ops: Vec<AlignOp> = Vec::new();
    let mut i = n;
    let mut j = m;

    while i > 0 || j > 0 {
        match state {
            Trace::M => {
                if i == 0 || j == 0 {
                    // Edge: should not happen in normal NW, but handle gracefully
                    if i > 0 {
                        raw_ops.push(AlignOp::InsRef);
                        i -= 1;
                    } else {
                        raw_ops.push(AlignOp::InsQry);
                        j -= 1;
                    }
                } else {
                    let prev = tm[i * cols + j];
                    raw_ops.push(AlignOp::Match);
                    i -= 1;
                    j -= 1;
                    state = prev;
                }
            }
            Trace::Ix => {
                let prev = tx[i * cols + j];
                raw_ops.push(AlignOp::InsRef);
                i -= 1;
                state = prev;
            }
            Trace::Iy => {
                let prev = ty[i * cols + j];
                raw_ops.push(AlignOp::InsQry);
                j -= 1;
                state = prev;
            }
        }
    }

    raw_ops.reverse();

    // Run-length encode
    let mut ops: Vec<(AlignOp, usize)> = Vec::new();
    for op in raw_ops {
        if let Some(last) = ops.last_mut()
            && last.0 == op {
                last.1 += 1;
                continue;
            }
        ops.push((op, 1));
    }

    PairwiseAlignment { score, ops }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn align(r: &[u8], q: &[u8]) -> PairwiseAlignment {
        needleman_wunsch(r, q, &AlignScoring::default())
    }

    #[test]
    fn test_nw_identical_sequences() {
        let a = align(b"ACGTACGT", b"ACGTACGT");
        assert_eq!(a.score, 16); // 8 * 2
        assert_eq!(a.ops, vec![(AlignOp::Match, 8)]);
    }

    #[test]
    fn test_nw_single_mismatch() {
        let a = align(b"ACGT", b"ACTT");
        // 3 matches (6) + 1 mismatch (-2) = 4
        assert_eq!(a.score, 4);
        assert_eq!(a.ops, vec![(AlignOp::Match, 4)]);
    }

    #[test]
    fn test_nw_single_insertion() {
        // Query has extra base
        let a = align(b"ACGT", b"ACGGT");
        // Expect: 4 matches + 1 insertion = 4*2 + (-3 + -1) = 4
        assert_eq!(a.score, 4);
        // Should have some Match ops surrounding an InsQry
        let total_match: usize = a.ops.iter().filter(|o| o.0 == AlignOp::Match).map(|o| o.1).sum();
        let total_ins: usize = a.ops.iter().filter(|o| o.0 == AlignOp::InsQry).map(|o| o.1).sum();
        assert_eq!(total_match, 4);
        assert_eq!(total_ins, 1);
    }

    #[test]
    fn test_nw_single_deletion() {
        // Reference has extra base
        let a = align(b"ACGGT", b"ACGT");
        let total_match: usize = a.ops.iter().filter(|o| o.0 == AlignOp::Match).map(|o| o.1).sum();
        let total_del: usize = a.ops.iter().filter(|o| o.0 == AlignOp::InsRef).map(|o| o.1).sum();
        assert_eq!(total_match, 4);
        assert_eq!(total_del, 1);
    }

    #[test]
    fn test_nw_homopolymer_length_diff() {
        // AAACCC vs AAAACCC — 1 extra A in query
        let a = align(b"AAACCC", b"AAAACCC");
        let total_match: usize = a.ops.iter().filter(|o| o.0 == AlignOp::Match).map(|o| o.1).sum();
        let total_ins: usize = a.ops.iter().filter(|o| o.0 == AlignOp::InsQry).map(|o| o.1).sum();
        assert_eq!(total_match, 6);
        assert_eq!(total_ins, 1);
    }

    #[test]
    fn test_nw_affine_gap_prefers_extension() {
        // One long gap should be preferred over two short gaps
        // Ref: AAACCCGGG, Query: AAAGGG
        // One 3-base gap: 6*2 + (-3 + 3*-1) = 12 - 6 = 6
        // Two gaps (e.g., 1+2): worse due to two gap opens
        let a = align(b"AAACCCGGG", b"AAAGGG");
        assert_eq!(a.score, 6);
        let total_del: usize = a.ops.iter().filter(|o| o.0 == AlignOp::InsRef).map(|o| o.1).sum();
        assert_eq!(total_del, 3);
        // Should be one contiguous gap, not split
        let gap_runs: usize = a.ops.iter().filter(|o| o.0 == AlignOp::InsRef).count();
        assert_eq!(gap_runs, 1);
    }

    #[test]
    fn test_nw_empty_reference() {
        let a = align(b"", b"ACGT");
        // All insertions: gap_open + 4 * gap_extend = -3 + 4*-1 = -7
        assert_eq!(a.score, -7);
        assert_eq!(a.ops, vec![(AlignOp::InsQry, 4)]);
    }

    #[test]
    fn test_nw_empty_query() {
        let a = align(b"ACGT", b"");
        assert_eq!(a.score, -7);
        assert_eq!(a.ops, vec![(AlignOp::InsRef, 4)]);
    }

    #[test]
    fn test_nw_both_empty() {
        let a = align(b"", b"");
        assert_eq!(a.score, 0);
        assert!(a.ops.is_empty());
    }

    #[test]
    fn test_nw_case_insensitive() {
        let a = align(b"acgt", b"ACGT");
        assert_eq!(a.score, 8);
        assert_eq!(a.ops, vec![(AlignOp::Match, 4)]);
    }
}
