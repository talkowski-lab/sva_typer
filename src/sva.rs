use std::{env, path::{Path, PathBuf}};

use crate::hmm::HMM;
use crate::builder::*;
use crate::utils::*;
use crate::reader::*;

const HEXAMER_REPEAT: &str = "CCCTCT";
const VNTR_REPEATS: &[&str] = &[
        "GCCTCTGCCCGGCCGCCCAGTCTGGGAAGTGAGGAGC",
        "GCCCGGCCAGCCGCCCCGTCCGGGAGGAGGTGGGGGGGTCAGCCCCC",
        "GCCGCCCCGACCGGGAAGTGAGGAGCCCCTCTGCCCG"
    ];

// TODO: Are these starts 0-based or 1-based? (They are 0-based from the MSA file)
const SVA_TYPES: &[(&str, &str, usize, usize, usize)] = &[

    ("SVA_A", "DF000001067.hmm", 70, 434, 886),
    ("SVA_B", "DF000001068.hmm", 70, 429, 882),
    ("SVA_C", "DF000001069.hmm", 70, 430, 883),
    ("SVA_D", "DF000001070.hmm", 70, 430, 884),
    ("SVA_E", "DF000001071.hmm", 70, 426, 879),
    ("SVA_F", "DF000001072.hmm", 70, 419, 872),

];


pub fn gen_sva_model(settings: &HMMBuildSettings) -> HMM {


    let hexamer_hmm = create_HMM_from_motifs(
        &[HEXAMER_REPEAT],
        &["hex"],
        settings,
        "hexamer_region"
    );

    let vntr_hmm = create_HMM_from_motifs(
        VNTR_REPEATS,
        &["VNTR_1", "VNTR_2", "VNTR_3"],
        settings,
        "VNTR_region"
    );

    let skip1 = create_skip_state(settings, Some("skip1"));
    let skip2 = create_skip_state(settings, Some("skip2"));
    let skip3 = create_skip_state(settings, Some("skip3"));

    append_HMM(vec![skip1, hexamer_hmm, skip2, vntr_hmm, skip3])
}

pub fn gen_sva_model_with_custom_hexseq(settings: &HMMBuildSettings, hex_motifs: &[String]) -> HMM {

    let motif_names = hex_motifs.iter()
        .enumerate()
        .map(|(i, _)| format!("hex_{}", i+1))
        .collect::<Vec<_>>();

    // I'm sure there's a less stupid way of doing this
    let motif_names_ref = motif_names.iter().map(|s| s.as_str()).collect::<Vec<_>>();

    let hexamer_hmm = create_HMM_from_motifs(
        &hex_motifs.iter().map(|s| s.as_str()).collect::<Vec<_>>(),
        &motif_names_ref,
        settings,
        "hexamer_region"
    );

    let vntr_hmm = create_HMM_from_motifs(
        VNTR_REPEATS,
        &["VNTR_1", "VNTR_2", "VNTR_3"],
        settings,
        "VNTR_region"
    );

    let skip1 = create_skip_state(settings, Some("skip1"));
    let skip2 = create_skip_state(settings, Some("skip2"));
    let skip3 = create_skip_state(settings, Some("skip3"));

    append_HMM(vec![skip1, hexamer_hmm, skip2, vntr_hmm, skip3])
}

pub fn sva_hmm_dir() -> PathBuf {
    if let Ok(dir) = env::var("SVA_HMM_PATH") {
        PathBuf::from(dir)
    } else {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("ref")
    }

}

pub fn gen_sva_model_with_innerseq_all_families(settings: &HMMBuildSettings) -> HMM {
    let hexamer_hmm = create_HMM_from_motifs(
        &[HEXAMER_REPEAT],
        &["hex"],
        settings,
        "hexamer_region"
    );

    let vntr_hmm = create_HMM_from_motifs(
        VNTR_REPEATS,
        &["VNTR_1", "VNTR_2", "VNTR_3"],
        settings,
        "VNTR_region"
    );

    let hmm_dir = sva_hmm_dir();

    let alu_region = parallelize_HMM(
        SVA_TYPES.iter().map(
            |(elem_type, path, alu_start, alu_end, _sine_start)| {
                read_hmm_file(
                    &hmm_dir.join(*path),
                    Some(format!("{elem_type}_ALU").as_str()), 
                    Some(*alu_start), 
                    Some(*alu_end))
            }).collect::<Result<Vec<_>,_>>().unwrap_or_else(|e| panic!("{e}")),
        "ALU"
    );

    let sine_region = parallelize_HMM(
        SVA_TYPES.iter().map(
            |(elem_type, path, _alu_start, _alu_end, sine_start)| {
                read_hmm_file(
                    &hmm_dir.join(*path), 
                    Some(format!("{elem_type}_SINE").as_str()), 
                    Some(*sine_start), 
                    None)
            }).collect::<Result<Vec<_>,_>>().unwrap_or_else(|e| panic!("{e}")),
        "SINE"
    );


    let skip1 = create_skip_state(settings, Some("skip1"));
    let skip2 = create_skip_state(settings, Some("skip2"));

    append_HMM(vec![skip1, hexamer_hmm, alu_region, vntr_hmm, sine_region, skip2])
}

pub fn gen_sva_model_with_innerseq(settings: &HMMBuildSettings) -> HMM {
    let hexamer_hmm = create_HMM_from_motifs(
        &[HEXAMER_REPEAT],
        &["hex"],
        settings,
        "hexamer_region"
    );

    let vntr_hmm = create_HMM_from_motifs(
        VNTR_REPEATS,
        &["VNTR_1", "VNTR_2", "VNTR_3"],
        settings,
        "VNTR_region"
    );

    let (elem_type, path, alu_start, alu_end, sine_start) = SVA_TYPES[0];

    let hmm_dir = sva_hmm_dir();

    let alu_region = read_hmm_file(
                    &hmm_dir.join(path),
                    Some(format!("{elem_type}_ALU").as_str()), 
                    Some(alu_start), 
                    Some(alu_end)).unwrap();

    let sine_region = read_hmm_file(
                    &hmm_dir.join(path), 
                    Some(format!("{elem_type}_SINE").as_str()), 
                    Some(sine_start), 
                    None
    ).unwrap();


    let skip1 = create_skip_state(settings, Some("skip1"));
    let skip2 = create_skip_state(settings, Some("skip2"));

    append_HMM(vec![skip1, hexamer_hmm, alu_region, vntr_hmm, sine_region, skip2])
}


pub fn trim_loop_intervals(final_intervals: &mut Vec<(&str, Interval)>) {
    // This is overly verbose because of borrowing rules
    let (_, hex_interval) = final_intervals.iter().find(|(n, _)| n == &"hexamer_region").unwrap();
    let hex_start = hex_interval.start;
    let hex_stop = hex_interval.stop;

    let (_, vntr_interval) = final_intervals.iter().find(|(n, _)| n == &"VNTR_region").unwrap();
    let vntr_start = vntr_interval.start;
    let vntr_stop = vntr_interval.stop;


    let mut new_hex_start = hex_start;
    let mut new_hex_stop = hex_stop;
    let mut new_vntr_start = vntr_start;
    let mut new_vntr_stop = vntr_stop;
    for (n, int) in final_intervals.iter_mut() {
        if n == &"hexamer_region_skip" {
            if int.start == hex_start {
                *n = "skip";
                new_hex_start = int.stop;
            }
            if int.stop == hex_stop {
                *n = "skip";
                new_hex_stop = int.start;
            }
        }
        if n == &"VNTR_region_skip" {
            if int.start == vntr_start {
                *n = "skip";
                new_vntr_start = int.stop;
            }
            if int.stop == vntr_stop {
                *n = "skip";
                new_vntr_stop = int.start;
            }
        }
    }
    for (n, int) in final_intervals.iter_mut() {
        if n == &"hexamer_region" {
            int.start = new_hex_start;
            int.stop = new_hex_stop;
        }
        if n == &"VNTR_region" {
            int.start = new_vntr_start;
            int.stop = new_vntr_stop;
        }
    }

}



#[cfg(test)]
mod tests {
    use super::*;
    // use crate::utils::*;

    const SVA_F_SEQ: &str = "CTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCCTCTTTCCACGGTCTCCCTCTCATGCGGAGCCGAAGCTGGACTGTACTGCTGCCATCTCGGCTCACTGCAACCTCCCTGCCTGATTCTCCTGCCTCAGCCTGCCGAGTGCCTGCGATTGCAGGCACGCGCCGCCACGCCTGACTGGTTTTGGTGGAGACGGGGTTTCGCTGTGTTGGCCGGGCCGGTCTCCAGCCCCTAACCGCGAGTGATCCGCCAGCCTCGGCCTCCCGAGGTGCCGGGATTGCAGACGGAGTCTCGTTCACTCAGTGCTCAATGGTGCCCAGGCTGGAGTGCAGTGGCGTGATCTCGGCTCGCTACAACCTACACCTCCCAGCCGCCTGCCTTGGCCTCCCAAAGTGCCGAGATTGCAGCCTCTGCCCGGCCGCCGCCCCGTCTGGGAGGTGAGGAGCGCCTCTGCCCGGCCGCCCATCGTCTGGGANGTGAGGAGCCCCTCTGCCCGGCCGCCCCGTCTGGGAGGTGAGGAGCGCCTCCGCCCGGCCGCCGCCCCGTCCGGGAGGTGAGGAGCGTCTCCGCCCGGCCGCCCNCCGTCCGGGANGTGAGGAGCGCCTCCGCCCGGCCGCCCCGTCCGGGANGTGAGGAGCGCCTCCGCCCGGCCAGCCGCCCCGTCCGGGAGGTGGGGGGGTCAGCCCCCCGCCCGGCCAGCCGCCCCGTCCGGGAGGAGGTGGGGGGGTCAGCCCCCCGCCCGGCCAGCCGCCCCGTCCGGGAGGTGAGGGGCGCCTCTGCCCGGCCGCCCCTACTGGGAAGTGAGGAGCCCCTCTGCCCGGCCACCGCCCCGTCCGGGAGGTGTGCCCAACAGCTCATTGAGAACGGGCCAGGATGACAATGGCGGCTTTGTGGAATAGAAAGGCGGGAAAGGTGGGGAAAAGATTGAGAAATCGGATGGTTGCCGTGTCTGTGTAGAAAGAAGTAGACATGGGAGACTTTTCATTTTGTTCTGCACTAAGAAAAATTCCTCTGCCTTGGGATCCTGTTGATCTGTGACCTTACCCCCAACCCTGTGCTCTCTGAAACATGTGCTGTGTCCACTCAGGGTTAAATGGATTAAGGGCGGTGCAAGATGTGCTTTGTTAAACAGATGCTTGAAGGCAGCATGCTCGTTAAGAGTCATCACCAATCCCTAATCTCAAGTAATCAGGGACACAAACACTGCGGAAGGCCGCAGGGTCCTCTGCCTAGGAAAACCAGAGACCTTTGTTCACTTGTTTATCTGCTGACCTTCCCTCCACTATTGTCCCATGACCCTGCCAAATCCCCCTCTGTGAGAAACACCCAAGAATTATCAATAAAAAAAATNAAAAAAAAAA";

    #[test]
    fn sva_test() {
        let settings = HMMBuildSettings::default();
        let hmm = gen_sva_model(&settings);
        hmm.check_valid();
        let mut result = hmm.query(&sequence_to_bytes(SVA_F_SEQ));
        // trim_loop_intervals(&mut result);
        let mut writer = std::io::stdout();
        // pprint_intervals(&mut writer, result);
        // panic!();
    }

    #[test]
    fn complex_sva_test() {
        let settings = HMMBuildSettings::default();
        let hmm = gen_sva_model_with_innerseq_all_families(&settings);
        hmm.check_valid();
    }
}
