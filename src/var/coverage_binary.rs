//! Whole-genome 100 kb binary coverage track.
//!
//! Produces `<prefix>.coverage_100k.bin`, a compact fixed-grid record of
//! primary-alignment read starts per 100 kb bin across every contig in the
//! BAM header. Unlike `.coverage.tsv` (1 Mb, mask-adjusted, gap-filtered),
//! this track is:
//!   - **Raw counts**: primary alignment starts only, no mask filter, no
//!     normalization. Every bin position is present, in coordinate order.
//!   - **Fixed grid**: bins anchored at 0 per contig, size 100 kb.
//!   - **Small**: ~125 KB uncompressed for T2T-CHM13 (24 contigs, ~31k bins).
//!
//! Intended for client-side visualization (drop-zone viewer): the layout is
//! trivial to parse into a `Uint32Array` and index by
//! `chrom_offset + pos / bin_size`.
//!
//! Binary layout (little-endian throughout):
//!
//! ```text
//! offset  size  field
//! 0       8     magic         b"NVCOVBIN"
//! 8       4     version       u32 = 1
//! 12      4     bin_size_bp   u32 = 100_000
//! 16      2     n_contigs     u16
//! 18+           per contig:
//!                 1  name_len (u8)
//!                 N  name (utf-8)
//!                 4  length_bp (u32)
//!                 4  n_bins (u32)   // = ceil(length_bp / bin_size)
//! ...           data: u32[] counts (sum(n_bins) entries).
//!               contig i's counts start at cumulative sum of prior n_bins.
//! ```

use std::fs::File;
use std::io::{BufWriter, Write};

use log::info;

use crate::input::{AlignmentHeader, AlignmentRecord};

pub const BIN_SIZE: usize = 100_000;
const MAGIC: &[u8; 8] = b"NVCOVBIN";
const VERSION: u32 = 1;

/// Per-contig counts. Length is `ceil(length_bp / BIN_SIZE)`.
struct ContigBins {
    name: String,
    length_bp: u32,
    counts: Vec<u32>,
}

pub struct BinaryCoverage100k {
    contigs: Vec<ContigBins>,
}

impl BinaryCoverage100k {
    pub fn new(header: &AlignmentHeader) -> Self {
        let contigs = header.refs.iter().zip(header.lengths.iter())
            .map(|(name, &len_i32)| {
                let length_bp = len_i32.max(0) as u32;
                let n_bins = (length_bp as usize).div_ceil(BIN_SIZE).max(1);
                ContigBins {
                    name: name.clone(),
                    length_bp,
                    counts: vec![0u32; n_bins],
                }
            })
            .collect();
        Self { contigs }
    }

    /// Increment the bin containing this record's alignment start, if the
    /// record is a primary alignment (not unmapped, secondary, or
    /// supplementary). Uses 0-based `pos` to match `plot_focal_depth`'s
    /// bin edges (`[k * BIN_SIZE, (k+1) * BIN_SIZE)`).
    pub fn process(&mut self, record: &AlignmentRecord) {
        // 0x004 unmapped, 0x100 secondary, 0x800 supplementary
        if record.flag & 0x904 != 0 { return; }
        if record.ref_id < 0 || record.pos < 0 { return; }
        let id = record.ref_id as usize;
        let Some(chrom) = self.contigs.get_mut(id) else { return; };
        let bin = (record.pos as usize) / BIN_SIZE;
        if bin < chrom.counts.len() {
            chrom.counts[bin] += 1;
        }
    }

    pub fn write(&self, path: &str) -> std::io::Result<()> {
        let f = File::create(path)?;
        let mut w = BufWriter::new(f);

        w.write_all(MAGIC)?;
        w.write_all(&VERSION.to_le_bytes())?;
        w.write_all(&(BIN_SIZE as u32).to_le_bytes())?;
        w.write_all(&(self.contigs.len() as u16).to_le_bytes())?;

        for c in &self.contigs {
            let name_bytes = c.name.as_bytes();
            // The u8 name length caps ref names at 255 bytes; SAM permits up
            // to 254 (per the spec's regex), so we're safely inside.
            debug_assert!(name_bytes.len() <= u8::MAX as usize);
            w.write_all(&[name_bytes.len() as u8])?;
            w.write_all(name_bytes)?;
            w.write_all(&c.length_bp.to_le_bytes())?;
            w.write_all(&(c.counts.len() as u32).to_le_bytes())?;
        }

        for c in &self.contigs {
            for &count in &c.counts {
                w.write_all(&count.to_le_bytes())?;
            }
        }

        w.flush()?;
        let total_bins: usize = self.contigs.iter().map(|c| c.counts.len()).sum();
        info!(
            "[coverage-100k] wrote {} ({} contigs, {} bins, {} bytes data)",
            path, self.contigs.len(), total_bins, total_bins * 4,
        );
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::input::AlignmentRecord;

    fn make_header(refs: &[(&str, i32)]) -> AlignmentHeader {
        AlignmentHeader {
            text: String::new(),
            refs: refs.iter().map(|(n, _)| n.to_string()).collect(),
            lengths: refs.iter().map(|(_, l)| *l).collect(),
        }
    }

    fn rec(ref_id: i32, pos: i32, flag: u16) -> AlignmentRecord {
        AlignmentRecord {
            ref_id, pos, flag,
            ..AlignmentRecord::default()
        }
    }

    #[test]
    fn bins_count_primary_starts_and_skip_others() {
        let h = make_header(&[("chr1", 250_000), ("chr2", 100_000)]);
        let mut bc = BinaryCoverage100k::new(&h);

        // primary read at chr1:50 -> bin 0
        bc.process(&rec(0, 50, 0));
        // primary read at chr1:100_000 -> bin 1
        bc.process(&rec(0, 100_000, 0));
        // supplementary read at chr1:50 -> skipped
        bc.process(&rec(0, 50, 0x800));
        // secondary read at chr1:50 -> skipped
        bc.process(&rec(0, 50, 0x100));
        // unmapped -> skipped
        bc.process(&rec(0, 50, 0x004));
        // primary read at chr2:0 -> chr2 bin 0
        bc.process(&rec(1, 0, 0));

        assert_eq!(bc.contigs[0].counts, vec![1, 1, 0]);
        assert_eq!(bc.contigs[1].counts, vec![1]);
    }

    #[test]
    fn write_roundtrips_header_and_counts() {
        let h = make_header(&[("chrA", 100_000), ("chrB", 250_000)]);
        let mut bc = BinaryCoverage100k::new(&h);
        bc.process(&rec(0, 0, 0));
        bc.process(&rec(1, 200_000, 0));

        let tmp = std::env::temp_dir().join(format!("nasvar_covbin_{}.bin", std::process::id()));
        let path = tmp.to_str().unwrap();
        bc.write(path).unwrap();
        let bytes = std::fs::read(path).unwrap();
        std::fs::remove_file(path).ok();

        assert_eq!(&bytes[0..8], MAGIC);
        assert_eq!(u32::from_le_bytes(bytes[8..12].try_into().unwrap()), 1);
        assert_eq!(u32::from_le_bytes(bytes[12..16].try_into().unwrap()), 100_000);
        assert_eq!(u16::from_le_bytes(bytes[16..18].try_into().unwrap()), 2);

        // chrA header: len=4, "chrA", length_bp=100000, n_bins=1
        assert_eq!(bytes[18], 4);
        assert_eq!(&bytes[19..23], b"chrA");
        assert_eq!(u32::from_le_bytes(bytes[23..27].try_into().unwrap()), 100_000);
        assert_eq!(u32::from_le_bytes(bytes[27..31].try_into().unwrap()), 1);

        // chrB header: len=4, "chrB", length_bp=250000, n_bins=3
        assert_eq!(bytes[31], 4);
        assert_eq!(&bytes[32..36], b"chrB");
        assert_eq!(u32::from_le_bytes(bytes[36..40].try_into().unwrap()), 250_000);
        assert_eq!(u32::from_le_bytes(bytes[40..44].try_into().unwrap()), 3);

        // Data: chrA[0]=1, chrB[0]=0, chrB[1]=0, chrB[2]=1
        let data = &bytes[44..];
        assert_eq!(data.len(), 4 * 4);
        assert_eq!(u32::from_le_bytes(data[0..4].try_into().unwrap()), 1);
        assert_eq!(u32::from_le_bytes(data[4..8].try_into().unwrap()), 0);
        assert_eq!(u32::from_le_bytes(data[8..12].try_into().unwrap()), 0);
        assert_eq!(u32::from_le_bytes(data[12..16].try_into().unwrap()), 1);
    }
}
