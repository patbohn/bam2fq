use std::fs::{File, create_dir_all};
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::str::from_utf8;
//use std::io::{BufRead, BufReader};
//use flate2::read::MultiGzDecoder;

use flate2::write::GzEncoder;
use flate2::Compression;

use structopt::StructOpt;

use bam::BamReader;
use bam::record::Record;
use bam::record::tags::{TagValue};

///get phred (intergers) to error probability table
pub mod phred_int_to_prob;
use phred_int_to_prob::PHRED_TO_ERROR_PROB;

#[derive(Debug, StructOpt)]
#[structopt(name = "bam2fq", about = "Converts unmapped bam output of the dorado basecaller into separate simplex and duplex fastq.gz files and calculates QC metrics in parallel. Also outputs poly A tail estimate located in pt tag.")]
struct Config {
    #[structopt(parse(from_os_str), help = "Input fastq file")]
    input_file: PathBuf,

    #[structopt(short = "o", long, help = "Output prefix of simplex and duplex fastq.gz files", default_value ="fastq_output")]
    output_prefix: String,
    
    #[structopt(short="s", long, help = "Output prefix of per-read stats and qscore histogram", default_value ="stats_output")]
    output_stats_prefix: String,
    
    #[structopt(short="m", long, help = "[Optional] minimal read length to output filtered reads", default_value ="0")]
    min_read_length: usize,
    
    #[structopt(short="q", long, help = "[Optional] minimal mean qscore to output filtered reads", default_value ="0")]
    min_mean_qscore: f64,

}

const MAX_QSCORE: usize = 94;
const THREADS: u16 = 8;

fn main() -> std::io::Result<()> {
    
    let config = Config::from_args();
    
    //define output paths
    let per_read_filename = format!("{}_read_stats.csv.gz", config.output_stats_prefix);
    let hist_filename = format!("{}_all_phred_hist.csv", config.output_stats_prefix);
    let fq_simplex_filename = format!("{}_simplex.fastq.gz", config.output_prefix);
    
    
    let min_read_length: usize = config.min_read_length;
    let min_mean_qscore: f64 = config.min_mean_qscore;
    
    //create output files
    let per_read_statsfile = File::create(per_read_filename)?;
    let mut stats_writer = GzEncoder::new(per_read_statsfile, Compression::fast());
    writeln!(&mut stats_writer, "read_id,read_length,mean_phred,mean_error_rate,poly_a_estimate,duplex_state,filtering_passed")?;
    
    let fq_simplex_file = File::create(fq_simplex_filename)?;
    let mut simplex_writer = GzEncoder::new(fq_simplex_file, Compression::fast());
    

    let fq_duplex_filename = format!("{}_duplex.fastq.gz", config.output_prefix);
    let fq_duplex_file = File::create(fq_duplex_filename)?;
    let mut duplex_writer = GzEncoder::new(fq_duplex_file, Compression::fast());


    //read in bam file
    let bam_reader = BamReader::from_path(config.input_file, THREADS).unwrap();
    //read bam let reader = BufReader::new(MultiGzDecoder::new(input_file));
    
    //create array to count phred qscores
    let mut simplex_p_qscore_hist: [u64; MAX_QSCORE] = [0;MAX_QSCORE];
    let mut simplex_f_qscore_hist: [u64; MAX_QSCORE] = [0;MAX_QSCORE];

    let mut duplex_p_qscore_hist: [u64; MAX_QSCORE] = [0;MAX_QSCORE];
    let mut duplex_f_qscore_hist: [u64; MAX_QSCORE] = [0;MAX_QSCORE];
    
    for record in bam_reader {
        //extract read_id, sequence, quality scores per nt, simplex/duplex state and run_id
        let this_record = record.unwrap();
        
        let id = from_utf8(this_record.name()).unwrap_or("NA");
        let raw_seq = this_record.sequence().to_vec();
        let sequence = String::from_utf8_lossy(&raw_seq);
        let qualities = this_record.qualities().raw();
        let duplex_state = get_duplex_tag(&this_record);
        let raw_run_id = get_optional_string_BAM_tag(&this_record, b"RG", "NA".to_owned());
        let run_id = raw_run_id.split("_").next().unwrap_or("");
        let poly_a_estimate = get_optional_int_BAM_tag(&this_record, b"pt", -1);
        
        //(re-)calculate per read mean accuracy (to be seen whether we want to re-evaluate it after trimming?
        let (mean_error_prob, read_length) = calc_mean_median_error(&qualities);
        let mean_quality = error_prob_to_phred(mean_error_prob);
        
        let filtering_passed: bool = (sequence.len() >= min_read_length) & (mean_quality >= min_mean_qscore);
        //write per read stats to file
        writeln!(&mut stats_writer, "{},{},{:.1},{:1.2e},{},{},{}", id, read_length, mean_quality, mean_error_prob, poly_a_estimate, duplex_state, filtering_passed)?;
        
        //write read to simplex or duplex file depending on duplex state
        if filtering_passed {
            if duplex_state == 1 {
                writeln!(&mut duplex_writer, "@{} run_id={} duplex={} poly_A_length={}\n{}\n+\n{}", id, run_id, duplex_state, poly_a_estimate, sequence, phred_to_utf8(qualities))?;
                for c in qualities.iter() {
                    if c <= &(MAX_QSCORE as u8) {
                        duplex_p_qscore_hist[*c as usize] += 1
                    } else {
                        panic!("QScore outside of expected range!")
                    }
                }
                
            } else {
                writeln!(&mut simplex_writer, "@{} run_id={} duplex={} poly_A_length={}\n{}\n+\n{}", id, run_id, duplex_state, poly_a_estimate, sequence, phred_to_utf8(qualities))?;
                for c in qualities.iter() {
                    if c <= &(MAX_QSCORE as u8) {
                        simplex_p_qscore_hist[*c as usize] += 1
                    } else {
                        panic!("QScore outside of expected range!")
                    }
                }
                
            }
        } else {
            if duplex_state == 1 {
                for c in qualities.iter() {
                    if c <= &(MAX_QSCORE as u8) {
                        duplex_f_qscore_hist[*c as usize] += 1
                    } else {
                        panic!("QScore outside of expected range!")
                    }
                }
            } else {
                for c in qualities.iter() {
                    if c <= &(MAX_QSCORE as u8) {
                        simplex_f_qscore_hist[*c as usize] += 1
                    } else {
                        panic!("QScore outside of expected range!")
                    }
                }
            }
            
        }
    }
    
    //write phred qscore histogram
    let hist_output_file = File::create(hist_filename)?;
    let mut writer = BufWriter::new(hist_output_file);
    writeln!(&mut writer, "phred_score,simplex_pass_count,duplex_pass_count,simplex_fail_count,duplex_fail_count")?;

    for qscore in 0..MAX_QSCORE {
        writeln!(&mut writer, "{},{},{},{},{}", qscore,simplex_p_qscore_hist[qscore],duplex_p_qscore_hist[qscore],simplex_f_qscore_hist[qscore],duplex_f_qscore_hist[qscore])?;
        }
    Ok(())
}

fn get_duplex_tag(record : &Record) -> i8 {
    if let Some(tag_value) = record.tags().get(b"dx") {
        match tag_value {
            TagValue::Int(tag_int, _) => tag_int as i8,
            _ => {panic!("Unexpected tag type");}
            }
    } else {
        0
    }
}

fn get_optional_string_BAM_tag(record : &Record, tag: &[u8;2], na_val : String) -> String {
    if let Some(tag_value) = record.tags().get(tag) {
        match tag_value {
            TagValue::String(tag_string, _) => String::from_utf8_lossy(tag_string).to_string(),
            _ => panic!("Unexpected tag type"),
            }
    } else {
        na_val
    }
}

fn get_optional_int_BAM_tag(record : &Record, tag: &[u8;2], na_val : i32) -> i32 {
    if let Some(tag_value) = record.tags().get(tag) {
        match tag_value {
            TagValue::Int(tag_int, _) => tag_int as i32,
            _ => panic!("Unexpected tag type"),
            }
    } else {
        na_val
    }
}

fn calc_mean_median_error(quality_array: &[u8]) -> (f64, i64) {
    let mut total_prob = 0.0;
    let mut count = 0;
    
    for quality_score in quality_array.iter() {
        let phred = *quality_score as usize;
        let prob = PHRED_TO_ERROR_PROB[phred];
        total_prob += prob;
        count += 1;
    }
    return (total_prob / count as f64, count as i64);
}

fn error_prob_to_phred(prob: f64) -> f64 {
    return -10.0_f64 * prob.log10()
}

fn phred_to_utf8(quality_array: &[u8]) -> std::string::String {
    let mut offset_array: Vec<u8> = Vec::with_capacity(quality_array.len());
    
    for &byte in quality_array {
        offset_array.push(byte+33);
    }
    
    let quality_string = match String::from_utf8(offset_array) {
        Ok(s) => s,
        Err(_e) => {
            panic!("Error converting quality scores to String")
        }
    };
    return quality_string
}