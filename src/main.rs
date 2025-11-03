use anyhow::Result;
use sva_typer::{
    builder::HMMBuildSettings,
    cli::{Args, SVAModelType}, 
    utils::*,
    sva
};
use clap::Parser;
use bio::{self, io::fasta::{self, FastaRead}};
use rayon::prelude::*;

fn run(args: Args) -> Result<()> {
    let mut reader = fasta::Reader::from_file(&args.file)?;
    let settings = HMMBuildSettings::try_from(&args)?;
    let hmm = match args.sva_model {
        SVAModelType::Simple => sva::gen_sva_model(&settings),
        SVAModelType::Complex => sva::gen_sva_model_with_innerseq(&settings),
        SVAModelType::ComplexAllFamilies => sva::gen_sva_model_with_innerseq_all_families(&settings)
    };
    sva::gen_sva_model(&settings);

    if args.cores == 1 {
        let mut writer = open_write(args.output_file.as_deref())?;
        writeln!(writer, "ID\tregion\tstart\tend")?;
        // TODO: Turn this into a parallel loop
        for (i, record) in reader.records().enumerate() {
            eprint!("Record {}\r", i);
            let record = record?;
            let query = std::str::from_utf8(record.seq()).unwrap().to_uppercase();
            let result = hmm.query(&sequence_to_bytes(&query));
            // sva::trim_loop_intervals(&mut result);
            tsvprint_intervals(&mut writer,record.id(), result)?;
        }
    } else {
        let mut writer = open_write(args.output_file.as_deref())?;
        writeln!(writer, "ID\tregion\tstart\tend")?;
        let mut total_i = 0;

        let mut record = fasta::Record::new();

        let batch_size = 1000;
        loop {
            let mut batch = Vec::with_capacity(batch_size);
            for i in 0..batch_size {
                // eprintln!("{i}");
                reader.read(&mut record)?;
                if record.is_empty() {
                    break
                }
                
                batch.push(record.clone());
                total_i += 1;
                eprint!("Record {}\r", total_i);
            }
            if batch.is_empty() {
                break
            }

            let results = batch.par_chunks(75).map(|vr| {
                vr.iter().map(|r| {
                    let query = std::str::from_utf8(r.seq()).unwrap().to_uppercase();
                    let result = hmm.query(&sequence_to_bytes(&query));
                    result
                }).collect::<Vec<_>>()
            }).flatten().collect::<Vec<_>>();


            for (record, result )in std::iter::zip(batch, results) {
                tsvprint_intervals(&mut writer,record.id(), result)?;
            }
        }
    }
    Ok(())



}

fn main() {
    if let Err(e) = run(Args::parse()) {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    }
}
