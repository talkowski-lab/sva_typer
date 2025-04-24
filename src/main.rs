use anyhow::Result;
use sva_typer::builder::{create_pHMM, HMMBuildSettings};
use sva_typer::utils::*;
use sva_typer::cli::Args;
use clap::Parser;


fn run(args: Args) -> Result<()> {
    let seq = sequence_to_bytes("ACGT");
    let settings: HMMBuildSettings = args.try_into()?;
    let mut hmm = create_pHMM(&seq, settings, None);
    let query = sequence_to_bytes("AGTTTGT");
    let result = hmm.query(&query);
    println!("Result: {:?}", result);

    Ok(())
}

fn main() {
    if let Err(e) = run(Args::parse()) {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    }
}
