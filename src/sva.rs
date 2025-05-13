use crate::hmm::HMM;
use crate::builder::*;
use crate::utils::*;

pub fn gen_sva_model(settings: &HMMBuildSettings) -> HMM {

    let hexamer_repeat = "CCCTCT";

    let vntr_repeats = vec![
        "GCCTCTGCCCGGCCGCCCAGTCTGGGAAGTGAGGAGC",
        "GCCCGGCCAGCCGCCCCGTCCGGGAGGAGGTGGGGGGGTCAGCCCCC",
        "GCCGCCCCGACCGGGAAGTGAGGAGCCCCTCTGCCCG"
    ];

    let hexamer_hmm = create_HMM_from_motifs(
        vec![hexamer_repeat],
        vec!["hex"],
        settings,
        "hexamer_region"
    );

    let vntr_hmm = create_HMM_from_motifs(
        vntr_repeats,
        vec!["VNTR_1", "VNTR_2", "VNTR_3"],
        settings,
        "VNTR_region"
    );
    append_HMM(vec![hexamer_hmm, vntr_hmm])
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
        trim_loop_intervals(&mut result);
        let mut writer = std::io::stdout();
        pprint_intervals(&mut writer, result);
        panic!();
    }
}
