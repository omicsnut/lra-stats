use arrow;
use rust_htslib::bam::record::{Aux, Cigar};
use rust_htslib::{bam, bam::Read, htslib};
use serde::{Deserialize, Serialize};
use serde_json;

use std::fs::File;
use std::sync::Arc;

pub struct ReadLevelStats {
    pub read_ids: Vec<String>,
    pub chrom_ids: Vec<i32>,
    pub avg_quals: Vec<f64>,
    pub lengths: Vec<u64>,
    pub map_quals: Vec<u8>,
    pub identities: Vec<f64>,
    pub sam_flags: Vec<u16>,
    pub num_pass_duplex_reads: u64,
    pub pct_duplex_reads: f64,
}

pub fn extract(bam: &mut bam::Reader) -> ReadLevelStats {
    let mut result = ReadLevelStats {
        read_ids: vec![],
        chrom_ids: vec![],
        avg_quals: vec![],
        lengths: vec![],
        map_quals: vec![],
        identities: vec![],
        sam_flags: vec![],
        num_pass_duplex_reads: 0,
        pct_duplex_reads: 0.0,
    };

    let mut duplex_read_count: u64 = 0;
    let mut total_read_count: u64 = 0;

    for read in bam.rc_records().map(|r| r.unwrap()) {
        result.sam_flags.push(read.flags());
        if read.is_unmapped()
            || read.is_secondary()
            || read.is_quality_check_failed()
            || read.is_duplicate()
        {
            continue;
        }

        match read.aux(b"dx") {
            Ok(value) => {
                if value == Aux::I8(1)
                    || value == Aux::I16(1)
                    || value == Aux::I32(1)
                    || value == Aux::U8(1)
                    || value == Aux::U16(1)
                    || value == Aux::U32(1)
                    || value == Aux::Float(1.0)
                    || value == Aux::Double(1.0)
                {
                    duplex_read_count += 1;
                    total_read_count += 1;
                } else if value == Aux::I8(0)
                    || value == Aux::I16(0)
                    || value == Aux::I32(0)
                    || value == Aux::U8(0)
                    || value == Aux::U16(0)
                    || value == Aux::U32(0)
                    || value == Aux::Float(0.0)
                    || value == Aux::Double(0.0)
                {
                    total_read_count += 1;
                }
            }

            Err(e) => {
                panic!("Error reading dx aux field: {}", e)
            }
        }

        result.identities.push(gap_compressed_identity(&read));
        result.map_quals.push(read.mapq());
        result.lengths.push(read.seq_len() as u64);
        result.avg_quals.push(average_read_quality(&read));
        result.chrom_ids.push(read.tid());
        result.read_ids.push(
            String::from_utf8(read.qname().to_vec())
                .expect("could not parse read name in bam/cram"),
        );
    }

    result.num_pass_duplex_reads = duplex_read_count;
    result.pct_duplex_reads = 100.0 * (duplex_read_count as f64) / (total_read_count as f64);

    result
}

#[derive(Serialize, Deserialize)]
pub struct SamFlagCounts {
    pub total: usize,
    pub mapped: usize,
    pub qcfail: usize,
    pub duplicate: usize,
    pub unmapped: usize,
    pub primary: usize,
    pub secondary: usize,
    pub supplementary: usize,
}

#[derive(Serialize, Deserialize)]
pub struct AggregateStats {
    pub num_pass_reads: usize,
    pub num_pass_duplex_reads: u64,
    pub read_length_n50: u64,
    pub bases_yield: u64,
    pub mean_read_length: f64,
    pub median_read_length: f64,
    pub mean_pct_identity: f64,
    pub median_pct_identity: f64,
    pub pct_duplex_reads: f64,
    pub sam_flag_counts: SamFlagCounts,
}

impl ReadLevelStats {
    pub fn aggregate(&self) -> AggregateStats {
        if self.read_ids.len() < 2 {
            panic!("Not enough reads to aggregate read level metrics. Exiting.");
        }

        let identities_sum = self.identities.iter().sum::<f64>();
        let num_pass_reads = self.read_ids.len();
        let bases_yield = self.lengths.iter().sum::<u64>();
        let read_length_n50 = compute_n50(self.lengths.clone().as_mut(), bases_yield);
        let mean_read_length = round_to_digit(bases_yield as f64 / num_pass_reads as f64, 2);
        let median_read_length = round_to_digit(compute_median_u64(self.lengths.as_slice()), 2);
        let mean_pct_identity = round_to_digit(identities_sum / num_pass_reads as f64, 4);
        let median_pct_identity = round_to_digit(compute_median(self.identities.as_slice()), 4);
        let pct_duplex_reads = round_to_digit(self.pct_duplex_reads, 4);

        let mut sam_flag_counts = SamFlagCounts {
            total: 0,
            mapped: 0,
            qcfail: 0,
            duplicate: 0,
            unmapped: 0,
            primary: 0,
            secondary: 0,
            supplementary: 0,
        };

        for sam_flag in self.sam_flags.iter() {
            sam_flag_counts.total += 1;

            let is_qcfail = sam_flag & htslib::BAM_FQCFAIL as u16 != 0;
            let is_duplicate = sam_flag & htslib::BAM_FDUP as u16 != 0;
            let is_supplementary = sam_flag & htslib::BAM_FSUPPLEMENTARY as u16 != 0;
            let is_secondary = sam_flag & htslib::BAM_FSECONDARY as u16 != 0;
            let is_unmapped = sam_flag & htslib::BAM_FUNMAP as u16 != 0;

            if is_qcfail {
                sam_flag_counts.qcfail += 1;
            }

            if is_duplicate {
                sam_flag_counts.duplicate += 1;
            }

            if is_supplementary {
                sam_flag_counts.supplementary += 1;
            }

            if is_secondary {
                sam_flag_counts.secondary += 1;
            } else {
                sam_flag_counts.primary += 1;
            }

            if is_unmapped {
                sam_flag_counts.unmapped += 1
            } else {
                sam_flag_counts.mapped += 1;
            }
        }

        AggregateStats {
            num_pass_reads,
            num_pass_duplex_reads: self.num_pass_duplex_reads,
            read_length_n50,
            bases_yield,
            mean_read_length,
            median_read_length,
            mean_pct_identity,
            median_pct_identity,
            sam_flag_counts,
            pct_duplex_reads,
        }
    }

    pub fn write_as_arrow(self, out_path: String) {
        use arrow::array::{Float64Array, StringArray, UInt64Array, UInt8Array};
        let read_ids = Arc::new(StringArray::from(self.read_ids)) as _;
        let avg_quals = Arc::new(Float64Array::from(self.avg_quals)) as _;
        let lengths = Arc::new(UInt64Array::from(self.lengths)) as _;
        let map_quals = Arc::new(UInt8Array::from(self.map_quals)) as _;
        let identities = Arc::new(Float64Array::from(self.identities)) as _;

        let batch = arrow::record_batch::RecordBatch::try_from_iter([
            ("readIDs", read_ids),
            ("quals", avg_quals),
            ("lengths", lengths),
            ("mapQ", map_quals),
            ("identities", identities),
        ])
        .unwrap();

        use arrow::datatypes::{DataType, Field, Schema};
        let schema = Schema::new(vec![
            Field::new("readIDs", DataType::Utf8, false),
            Field::new("quals", DataType::Float64, false),
            Field::new("lengths", DataType::UInt64, false),
            Field::new("mapQ", DataType::UInt8, false),
            Field::new("identities", DataType::Float64, false),
        ]);

        use arrow::ipc::writer::FileWriter as ArrowFileWriter;
        let buffer = File::create(out_path).expect("could not create output file");

        let mut writer =
            ArrowFileWriter::try_new(buffer, &schema).expect("could not create arrow writer");

        writer
            .write(&batch)
            .expect("could not write arrow record batch");

        writer
            .finish()
            .expect("could not finish writing arrow record batch");
    }
}

impl AggregateStats {
    pub fn write_as_json(&self, out_path: String) {
        let buffer = File::create(out_path).expect("could not create output file");
        serde_json::to_writer_pretty(&buffer, &self).expect("msg");
        buffer.sync_all().expect("msg");
    }
}

fn average_read_quality(read: &bam::Record) -> f64 {
    let qual_sum = read.qual().iter().map(|&i| i as u64).sum::<u64>();
    (qual_sum as f64) / (read.seq_len() as f64)
}

/// de:f tag added by newer versions of minimap2 is the
/// gap-compressed per-base sequence divergence. It can
/// be converted into percent identity – 100 * (1 - de)
/// https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity
fn get_de(read: &bam::Record) -> Option<f64> {
    match read.aux(b"de") {
        Ok(value) => match value {
            Aux::Float(v) => Some(v as f64),
            _ => panic!("unexpected type for de tag {:?}", value),
        },

        Err(_) => None,
    }
}

fn get_nm(read: &bam::Record) -> u32 {
    match read.aux(b"NM") {
        Ok(value) => match value {
            Aux::U8(v) => u32::from(v),
            Aux::U16(v) => u32::from(v),
            Aux::U32(v) => v,
            _ => panic!("unexpected type for NM tag {:?}", value),
        },

        Err(e) => panic!("could not parse value from NM tag: {:?}", e),
    }
}

/// Compute gap compressed identity from cigar data and NM tag
/// https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity
fn compute_gap_compressed_identity(read: &bam::Record) -> f64 {
    let mut matches = 0;
    let mut gap_size = 0;
    let mut gap_count = 0;

    for item in read.cigar().iter() {
        match item {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => matches += *len,
            Cigar::Del(len) | Cigar::Ins(len) => {
                gap_size += *len;
                gap_count += 1;
            }
            _ => (),
        }
    }

    1.0 - ((get_nm(&read) - gap_size + gap_count) as f64 / (matches + gap_count) as f64)
}

/// Get the gap compressed identity of the read from either the
/// de:f tag is it is present or compute it from cigar & NM tag
/// https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity
fn gap_compressed_identity(read: &bam::Record) -> f64 {
    match get_de(read) {
        Some(de_value) => 100.0 * (1.0 - de_value),
        None => compute_gap_compressed_identity(read),
    }
}

fn compute_median<T: Into<f64> + Copy>(array: &[T]) -> f64 {
    if (array.len() % 2) == 0 {
        let ind_left = array.len() / 2 - 1;
        let ind_right = array.len() / 2;
        (array[ind_left].into() + array[ind_right].into()) / 2.0
    } else {
        array[array.len() / 2].into()
    }
}

fn compute_median_u64(array: &[u64]) -> f64 {
    if (array.len() % 2) == 0 {
        let ind_left = array.len() / 2 - 1;
        let ind_right = array.len() / 2;
        (array[ind_left] + array[ind_right]) as f64 / 2.0
    } else {
        array[array.len() / 2] as f64
    }
}

fn compute_n50(lengths: &mut Vec<u64>, total_num_bases: u64) -> u64 {
    lengths.sort_unstable();
    let mut counter = 0;
    for val in lengths.iter() {
        counter += *val;
        if counter > total_num_bases / 2 {
            return *val;
        }
    }

    lengths[lengths.len() - 1]
}

fn round_to_digit(val: f64, prec: usize) -> f64 {
    const BASE: f64 = 10.0;
    let scaling = BASE.powi(prec as i32);
    f64::trunc(val * scaling) / scaling
}
