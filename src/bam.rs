//! BAM-related utilities for split alignment handling and CIGAR operations.

// Re-export ContigMapper for backward compatibility with existing code
pub use crate::utils::contig::{ContigMapper, NamingConvention};

/// CIGAR operation enum for raw CIGAR decoding
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CigarOp {
    Match,
    Ins,
    Del,
    RefSkip,
    SoftClip,
    HardClip,
    Pad,
    Equal,
    Diff,
    Back,
}

impl CigarOp {
    pub fn from_u32(val: u32) -> (Self, u32) {
        let op = val & 0xF;
        let len = val >> 4;
        let kind = match op {
            0 => CigarOp::Match,
            1 => CigarOp::Ins,
            2 => CigarOp::Del,
            3 => CigarOp::RefSkip,
            4 => CigarOp::SoftClip,
            5 => CigarOp::HardClip,
            6 => CigarOp::Pad,
            7 => CigarOp::Equal,
            8 => CigarOp::Diff,
            9 => CigarOp::Back,
            _ => CigarOp::Match, // Default/Error?
        };
        (kind, len)
    }
}

#[derive(Debug, Clone)]
pub struct AlignmentCoords {
    pub qs: u32,
    pub qe: u32,
    pub ts: u32,
    pub te: u32,
    pub is_reverse: bool,
}

/// Convert BAM CIGAR coordinates to original strand coordinates.
///
/// PAF and BAM files store query coordinates differently - BAM coordinates
/// have to be converted back to the original strand coordinates before comparison.
pub fn fix_sam_coords(ts: u32, te: u32, cigar: &[u32], flags: u16) -> std::io::Result<AlignmentCoords> {
    let mut qlen = 0u32;
    let mut hard_clip_start = 0u32;
    let mut soft_clip_start = 0u32;
    let mut seen_non_clip = false;

    for val in cigar {
        let (op, len) = CigarOp::from_u32(*val);
        match op {
            CigarOp::Match | CigarOp::Equal | CigarOp::Diff | CigarOp::Ins => {
                qlen += len;
                seen_non_clip = true;
            }
            CigarOp::SoftClip => {
                qlen += len;
                if !seen_non_clip {
                    soft_clip_start = len;
                }
                seen_non_clip = true;
            }
            CigarOp::HardClip => {
                qlen += len;
                if !seen_non_clip {
                    hard_clip_start = len;
                }
            }
            CigarOp::Del | CigarOp::RefSkip => {
                seen_non_clip = true;
            },
            _ => { seen_non_clip = true; }
        }
    }

    // query start = (hard clip at start) + (soft clip at start)
    let mut qs = hard_clip_start + soft_clip_start;
    let mut qe = qs;

    // Calculate qe by summing aligned query bases (M, I, =, X)
    for val in cigar {
        let (op, len) = CigarOp::from_u32(*val);
        match op {
            CigarOp::Match | CigarOp::Equal | CigarOp::Diff | CigarOp::Ins => {
                qe += len;
            }
            _ => {}
        }
    }

    let is_reverse = (flags & 0x10) != 0;

    // Flip to original strand coordinates if reverse
    if is_reverse {
        let tmp = qlen - qe;
        qe = qlen - qs;
        qs = tmp;
    }

    Ok(AlignmentCoords {
        qs,
        qe,
        ts,
        te,
        is_reverse,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn encode_cigar(op: u32, len: u32) -> u32 {
        (len << 4) | op
    }

    #[test]
    fn test_cigar_op_decode() {
        // Match (0) with length 10
        let (op, len) = CigarOp::from_u32(encode_cigar(0, 10));
        assert_eq!(op, CigarOp::Match);
        assert_eq!(len, 10);

        // Insertion (1) with length 5
        let (op, len) = CigarOp::from_u32(encode_cigar(1, 5));
        assert_eq!(op, CigarOp::Ins);
        assert_eq!(len, 5);

        // SoftClip (4) with length 100
        let (op, len) = CigarOp::from_u32(encode_cigar(4, 100));
        assert_eq!(op, CigarOp::SoftClip);
        assert_eq!(len, 100);
    }

    #[test]
    fn test_fix_sam_coords_simple() {
        // Simple case: 100M alignment on forward strand
        let cigar = vec![encode_cigar(0, 100)]; // 100M
        let result = fix_sam_coords(1000, 1100, &cigar, 0).unwrap();

        assert_eq!(result.qs, 0);
        assert_eq!(result.qe, 100);
        assert_eq!(result.ts, 1000);
        assert_eq!(result.te, 1100);
        assert!(!result.is_reverse);
    }

    #[test]
    fn test_fix_sam_coords_with_soft_clip() {
        // 10S90M - soft clip at start
        let cigar = vec![
            encode_cigar(4, 10),  // 10S
            encode_cigar(0, 90),  // 90M
        ];
        let result = fix_sam_coords(1000, 1090, &cigar, 0).unwrap();

        assert_eq!(result.qs, 10);
        assert_eq!(result.qe, 100);
        assert!(!result.is_reverse);
    }

    #[test]
    fn test_fix_sam_coords_reverse() {
        // 100M on reverse strand
        let cigar = vec![encode_cigar(0, 100)]; // 100M
        let result = fix_sam_coords(1000, 1100, &cigar, 0x10).unwrap();

        // On reverse strand, coordinates are flipped
        assert_eq!(result.qs, 0);
        assert_eq!(result.qe, 100);
        assert!(result.is_reverse);
    }
}
