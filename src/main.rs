use anyhow::Result;
use sva_typer::{
    builder::HMMBuildSettings,
    cli::Args, 
    utils::*,
    sva
};
use clap::Parser;
use bio::{self, io::fasta};

fn run(args: Args) -> Result<()> {
    let reader = fasta::Reader::from_file(&args.file)?;
    let settings = HMMBuildSettings::try_from(&args)?;
    let hmm = sva::gen_sva_model(&settings);

    let mut writer = open_write(args.output_file.as_deref())?;
    writeln!(writer, "ID\tregion\tstart\tend")?;
    // TODO: Turn this into a parallel loop
    for (i, record) in reader.records().enumerate() {
        eprint!("Record {}\r", i);
        let record = record?;
        let query = std::str::from_utf8(record.seq()).unwrap().to_uppercase();
        let mut result = hmm.query(&sequence_to_bytes(&query));
        sva::trim_loop_intervals(&mut result);
        tsvprint_intervals(&mut writer,record.id(), result)?;
    }
    Ok(())
}

fn main() {
    if let Err(e) = run(Args::parse()) {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    }
}
