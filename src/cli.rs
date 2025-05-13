use std::path::PathBuf;

use crate::builder::{HMMBuildSettings, HMMBuildError};
use anyhow::Result;
use clap::Parser;

fn between_0_1_parser(s: &str) -> Result<f64> {
    let val: f64 = s.parse()?;
    if (0.0..=1.0).contains(&val) {
        Ok(val)
    } else {
        Err(anyhow::anyhow!("Value must be between 0 and 1"))
    }
}
#[derive(Parser)]
pub struct Args {

    
    #[arg(value_name = "FILE")]
    pub file: PathBuf,
    /// Output file
    #[arg(short, long="output")]
    pub output_file: Option<PathBuf>,
    /// Probability of match state to match
    #[arg(
        long,
        default_value_t = HMMBuildSettings::default().match_to_match,
        value_parser=between_0_1_parser,
        help_heading = "HMM Build Parameters",
    )]
    pub match_to_match: f64,

    /// Probability of match state to insertion
    #[arg(
        long,
        default_value_t = HMMBuildSettings::default().match_to_ins,
        value_parser=between_0_1_parser,
        help_heading = "HMM Build Parameters",

    )]
    pub match_to_ins: f64,

    /// Probability of insertion state extension
    #[arg(
        long,
        default_value_t = HMMBuildSettings::default().ins_extend,
        value_parser=between_0_1_parser,
        help_heading = "HMM Build Parameters",

    )]
    pub ins_extend: f64,

    /// Probability of deletion state extension
    #[arg(
        long,
        default_value_t = HMMBuildSettings::default().del_extend,
        value_parser=between_0_1_parser,
        help_heading = "HMM Build Parameters",

    )]
    pub del_extend: f64,

    /// Probability of loop repeating
    #[arg(
        long,
        default_value_t = HMMBuildSettings::default().loop_prob,
        value_parser=between_0_1_parser,
        help_heading = "HMM Build Parameters",

    )]
    pub loop_prob: f64,
    /// Probability of skip state continuing
    #[arg(
        long,
        default_value_t = HMMBuildSettings::default().skip_to_skip,
        value_parser=between_0_1_parser,
        help_heading = "HMM Build Parameters",

    )]
    pub skip_to_skip: f64,

    /// Probability of match emission correctness
    #[arg(
        long,
        default_value_t = HMMBuildSettings::default().match_emit_correct,
        value_parser=between_0_1_parser,
        help_heading = "HMM Build Parameters",

    )]
    pub match_emit_correct: f64,

}

impl TryFrom<&Args> for HMMBuildSettings {
    type Error = HMMBuildError;

    fn try_from(value: &Args) -> std::result::Result<Self, Self::Error> {
        HMMBuildSettings::new(
            value.match_to_match,
            value.match_to_ins,
            value.ins_extend,
            value.del_extend,
            value.loop_prob,
            value.skip_to_skip,
            value.match_emit_correct,
        )
    }
}

