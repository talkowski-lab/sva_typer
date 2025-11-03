#[cfg(test)]
mod tests {
    use std::time;
    use std::path::Path;
    use std::time::Duration;
    use bio::io::fasta;
    use crate::sva;
    use crate::utils::*;
    use crate::builder::*;

    // TODO: Can I make a separate function that takes in a closure?
    
    #[test]
    fn bench_sva_seqs() {
        let count = 10;
        let settings = HMMBuildSettings::default();
        let sva_simple = sva::gen_sva_model(&settings);
        let sva_complex = sva::gen_sva_model_with_innerseq(&settings);
        let sva_complex_fam = sva::gen_sva_model_with_innerseq_all_families(&settings);

        let seq_file = Path::new(env!("CARGO_MANIFEST_DIR")).join("test/SVA_ref_core.fa");

        for (label, model) in [
            ("Simple", sva_simple),
            ("Complex", sva_complex),
            ("Complex_All_Families", sva_complex_fam)
        ] {
            let mut total_duration = Duration::from_secs(0);
            let mut total_seq_count = 0;
            for _ in 0..count {
                let mut reader = fasta::Reader::from_file(&seq_file).unwrap();
                
                let t = time::Instant::now();
                for record in reader.records() {
                    let record = record.unwrap();
                    let query = std::str::from_utf8(record.seq()).unwrap().to_uppercase();
                    let result = model.query(&sequence_to_bytes(&query));
                    total_seq_count += 1;
                }
                let d = t.elapsed();
                total_duration += d;
            }
            let avg_duration = total_duration / total_seq_count;
            eprintln!("{}: {}s {}us", label, avg_duration.as_secs(), avg_duration.subsec_micros());
        }
        panic!()


    }
    

}
