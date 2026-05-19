//! Streaming deduplication module -- integrates into the pipeline
//!
//! Uses a FxHashSet of read hashes to track seen sequences. First occurrence
//! passes, subsequent occurrences are failed as PCR duplicates.
//!
//! Pair-aware: in paired-end mode, R1 and R2 are hashed together as a single
//! key via `process_pair()`. Both mates receive the same duplicate verdict,
//! eliminating asymmetric dedup and concordance-induced viral read loss.
//!
//! Uses full read sequence (skipping first 5bp for trim tolerance) rather than
//! a short prefix, to avoid false collisions between reads sharing conserved
//! domains.
//!
//! Memory: ~16 bytes per unique read/pair. 100M reads = ~1.6 GB.

use crate::config::DedupConfig;
use crate::modules::{AtomicStats, ModuleReport, QcModule};
use crate::pipeline::AnnotatedRecord;
use rustc_hash::FxHashSet;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Mutex;

/// Number of bases to skip from read start (trim tolerance)
const SKIP_BASES: usize = 5;

/// Streaming deduplication module
pub struct StreamingDedup {
    /// Thread-safe set of seen read hashes
    seen: Mutex<FxHashSet<u64>>,
    /// Standard stats
    stats: AtomicStats,
    /// Duplicate reads removed
    duplicates_removed: AtomicU64,
}

impl StreamingDedup {
    pub fn new(_config: &DedupConfig) -> Self {
        Self {
            seen: Mutex::new(FxHashSet::default()),
            stats: AtomicStats::new(),
            duplicates_removed: AtomicU64::new(0),
        }
    }
}

impl QcModule for StreamingDedup {
    /// Single-end dedup: hash full sequence
    fn process(&self, record: &mut AnnotatedRecord) {
        self.stats.record_processed();

        let seq = &record.record.sequence;
        if seq.len() < SKIP_BASES + 10 {
            return; // too short to hash meaningfully
        }

        let hash = hash_sequence(seq);

        let is_new = {
            let mut set = self.seen.lock().unwrap();
            set.insert(hash)
        };

        if !is_new {
            record.fail("pcr_duplicate");
            self.duplicates_removed.fetch_add(1, Ordering::Relaxed);
            self.stats.record_removed();
        }
    }

    /// Paired-end dedup: hash R1+R2 together so both get the same verdict
    fn process_pair(&self, r1: &mut AnnotatedRecord, r2: &mut AnnotatedRecord) {
        // Count both mates as processed
        self.stats.record_processed();
        self.stats.record_processed();

        let seq1 = &r1.record.sequence;
        let seq2 = &r2.record.sequence;

        if seq1.len() < SKIP_BASES + 10 || seq2.len() < SKIP_BASES + 10 {
            return;
        }

        // Combine R1 and R2 into a single hash
        let hash = hash_pair(seq1, seq2);

        let is_new = {
            let mut set = self.seen.lock().unwrap();
            set.insert(hash)
        };

        if !is_new {
            r1.fail("pcr_duplicate");
            r2.fail("pcr_duplicate");
            self.duplicates_removed.fetch_add(2, Ordering::Relaxed);
            self.stats.record_removed();
            self.stats.record_removed();
        }
    }

    fn report(&self) -> ModuleReport {
        let seen_count = self.seen.lock().unwrap().len() as u64;
        self.stats.to_report(
            self.name(),
            serde_json::json!({
                "duplicates_removed": self.duplicates_removed.load(Ordering::Relaxed),
                "unique_sequences": seen_count,
                "estimated_library_complexity": seen_count,
            }),
        )
    }

    fn name(&self) -> &str {
        "dedup"
    }
}

/// Hash full read sequence (skip first 5 bases for trim tolerance)
fn hash_sequence(sequence: &[u8]) -> u64 {
    let start = SKIP_BASES.min(sequence.len());
    let mut hash: u64 = 0xcbf29ce484222325; // FNV-1a offset basis
    for &b in &sequence[start..] {
        hash ^= b.to_ascii_uppercase() as u64;
        hash = hash.wrapping_mul(0x100000001b3); // FNV-1a prime
    }
    hash
}

/// Hash a read pair: combine R1 and R2 sequences into one hash.
/// Uses a different initial seed for R2 to avoid symmetric collisions
/// (where R1_a == R2_b and R2_a == R1_b would produce same combined hash).
fn hash_pair(seq1: &[u8], seq2: &[u8]) -> u64 {
    let start1 = SKIP_BASES.min(seq1.len());
    let start2 = SKIP_BASES.min(seq2.len());

    // Hash R1
    let mut h1: u64 = 0xcbf29ce484222325;
    for &b in &seq1[start1..] {
        h1 ^= b.to_ascii_uppercase() as u64;
        h1 = h1.wrapping_mul(0x100000001b3);
    }

    // Hash R2 with different seed
    let mut h2: u64 = 0x6c62272e07bb0142;
    for &b in &seq2[start2..] {
        h2 ^= b.to_ascii_uppercase() as u64;
        h2 = h2.wrapping_mul(0x100000001b3);
    }

    // Combine: XOR with rotation to avoid commutativity
    h1 ^ h2.rotate_left(31)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::DedupConfig;
    use biometal::FastqRecord;

    fn make_record(seq: &[u8]) -> AnnotatedRecord {
        let qual = vec![b'I'; seq.len()];
        AnnotatedRecord::new(FastqRecord::new("test".into(), seq.to_vec(), qual))
    }

    fn make_dedup() -> StreamingDedup {
        StreamingDedup::new(&DedupConfig {
            enabled: false,
            optical_distance: 2500,
            umi_aware: false,
        })
    }

    #[test]
    fn test_first_occurrence_passes() {
        let dedup = make_dedup();
        let seq: Vec<u8> = (0..100).map(|i| [b'A', b'T', b'G', b'C'][i % 4]).collect();
        let mut record = make_record(&seq);
        dedup.process(&mut record);
        assert!(!record.is_failed());
    }

    #[test]
    fn test_duplicate_fails() {
        let dedup = make_dedup();
        let seq: Vec<u8> = (0..100).map(|i| [b'A', b'T', b'G', b'C'][i % 4]).collect();

        let mut first = make_record(&seq);
        dedup.process(&mut first);
        assert!(!first.is_failed());

        let mut second = make_record(&seq);
        dedup.process(&mut second);
        assert!(second.is_failed());
    }

    #[test]
    fn test_different_reads_pass() {
        let dedup = make_dedup();
        let seq1: Vec<u8> = (0..100).map(|i| [b'A', b'T', b'G', b'C'][i % 4]).collect();
        let seq2: Vec<u8> = (0..100).map(|i| [b'G', b'C', b'T', b'A'][i % 4]).collect();

        let mut r1 = make_record(&seq1);
        let mut r2 = make_record(&seq2);
        dedup.process(&mut r1);
        dedup.process(&mut r2);
        assert!(!r1.is_failed());
        assert!(!r2.is_failed());
    }

    #[test]
    fn test_same_prefix_different_suffix_not_duplicate() {
        // Previously with 50-base prefix, these would collide.
        // With full-length hashing, they should be unique.
        let dedup = make_dedup();
        let mut seq1: Vec<u8> = (0..125).map(|i| [b'A', b'T', b'G', b'C'][i % 4]).collect();
        let mut seq2 = seq1.clone();
        // Change only the last 10 bases
        for i in 115..125 {
            seq2[i] = b'N';
        }

        let mut r1 = make_record(&seq1);
        let mut r2 = make_record(&seq2);
        dedup.process(&mut r1);
        dedup.process(&mut r2);
        assert!(!r1.is_failed());
        assert!(!r2.is_failed(), "Reads differing only in suffix should NOT be called duplicates");
    }

    #[test]
    fn test_pair_dedup_both_fail() {
        let dedup = make_dedup();
        let seq_r1: Vec<u8> = (0..100).map(|i| [b'A', b'T', b'G', b'C'][i % 4]).collect();
        let seq_r2: Vec<u8> = (0..100).map(|i| [b'C', b'G', b'A', b'T'][i % 4]).collect();

        // First pair passes
        let mut p1_r1 = make_record(&seq_r1);
        let mut p1_r2 = make_record(&seq_r2);
        dedup.process_pair(&mut p1_r1, &mut p1_r2);
        assert!(!p1_r1.is_failed());
        assert!(!p1_r2.is_failed());

        // Duplicate pair: both mates should fail
        let mut p2_r1 = make_record(&seq_r1);
        let mut p2_r2 = make_record(&seq_r2);
        dedup.process_pair(&mut p2_r1, &mut p2_r2);
        assert!(p2_r1.is_failed(), "R1 of duplicate pair should fail");
        assert!(p2_r2.is_failed(), "R2 of duplicate pair should fail");
    }

    #[test]
    fn test_pair_dedup_different_pairs_pass() {
        let dedup = make_dedup();
        let seq_a_r1: Vec<u8> = (0..100).map(|i| [b'A', b'T', b'G', b'C'][i % 4]).collect();
        let seq_a_r2: Vec<u8> = (0..100).map(|i| [b'C', b'G', b'A', b'T'][i % 4]).collect();
        let seq_b_r1: Vec<u8> = (0..100).map(|i| [b'G', b'A', b'T', b'C'][i % 4]).collect();
        let seq_b_r2: Vec<u8> = (0..100).map(|i| [b'T', b'C', b'G', b'A'][i % 4]).collect();

        let mut a_r1 = make_record(&seq_a_r1);
        let mut a_r2 = make_record(&seq_a_r2);
        dedup.process_pair(&mut a_r1, &mut a_r2);

        let mut b_r1 = make_record(&seq_b_r1);
        let mut b_r2 = make_record(&seq_b_r2);
        dedup.process_pair(&mut b_r1, &mut b_r2);

        assert!(!a_r1.is_failed());
        assert!(!b_r1.is_failed());
    }

    #[test]
    fn test_pair_not_commutative() {
        // Swapping R1 and R2 should produce a different hash
        let dedup = make_dedup();
        let seq_a: Vec<u8> = (0..100).map(|i| [b'A', b'T', b'G', b'C'][i % 4]).collect();
        let seq_b: Vec<u8> = (0..100).map(|i| [b'C', b'G', b'A', b'T'][i % 4]).collect();

        let mut p1_r1 = make_record(&seq_a);
        let mut p1_r2 = make_record(&seq_b);
        dedup.process_pair(&mut p1_r1, &mut p1_r2);

        // Swapped: R1=seq_b, R2=seq_a — should NOT collide
        let mut p2_r1 = make_record(&seq_b);
        let mut p2_r2 = make_record(&seq_a);
        dedup.process_pair(&mut p2_r1, &mut p2_r2);

        assert!(!p2_r1.is_failed(), "Swapped pair should not be duplicate");
    }

    #[test]
    fn test_report_counts() {
        let dedup = make_dedup();
        let seq: Vec<u8> = (0..100).map(|i| [b'A', b'T', b'G', b'C'][i % 4]).collect();

        for _ in 0..5 {
            let mut r = make_record(&seq);
            dedup.process(&mut r);
        }

        let report = dedup.report();
        assert_eq!(report.reads_processed, 5);
        assert_eq!(report.reads_removed, 4); // 1 unique + 4 duplicates
    }
}
