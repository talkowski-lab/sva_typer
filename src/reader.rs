use std::{io::{self, BufRead, BufReader}, num::{ParseFloatError, ParseIntError}, ops::Neg, path::Path};
use crate::hmm::{HMM, HMMState};
use thiserror::Error;

#[derive(Error, Debug)]
pub enum ReaderError {
    #[error("IO error: {0}")]
    Io(#[from] io::Error),
    #[error("Parse error: {0}")]
    Parse(String),
    #[error("Parse float error: Could not parse \"{str}\" as float")]
    ParseFloat {
        str: String,
        source: ParseFloatError
    },
    #[error("Parse int error: {0}")]
    ParseInt(#[from] ParseIntError),
    #[error("File terminates early")]
    TruncatedFile

}

fn neg_vec<T: Neg>(x: Vec<T>) -> Vec<<T as Neg>::Output> {
    x.into_iter().map(|x| -x).collect::<Vec<_>>()
}

#[allow(clippy::too_many_arguments)]
fn from_vectors(
    prefix: Option<&str>,
    match_to_match_probs: Vec<f64>,
    match_to_ins_probs: Vec<f64>,
    match_to_del_probs: Vec<f64>,
    ins_to_match_probs: Vec<f64>,
    ins_to_ins_probs: Vec<f64>,
    del_to_match_probs: Vec<f64>,
    del_to_del_probs: Vec<f64>,
                                               
    match_emit_vecs: Vec<Vec<f64>>,
    ins_emit_vecs: Vec<Vec<f64>>
) -> HMM {


    // Make them all negative
    let match_to_match_probs = neg_vec(match_to_match_probs);
    let match_to_ins_probs = neg_vec(match_to_ins_probs);
    let match_to_del_probs = neg_vec(match_to_del_probs);
    let ins_to_match_probs = neg_vec(ins_to_match_probs);
    let ins_to_ins_probs = neg_vec(ins_to_ins_probs);
    let del_to_match_probs = neg_vec(del_to_match_probs);
    let del_to_del_probs = neg_vec(del_to_del_probs);
    let match_emit_vecs = match_emit_vecs.into_iter().map(neg_vec).collect::<Vec<_>>();
    let ins_emit_vecs = ins_emit_vecs.into_iter().map(neg_vec).collect::<Vec<_>>();



    let mut hmm = HMM::new();

    let prefix = match prefix {
        Some(p) => {
            if !p.ends_with('_') {
                format!("{}_", p)
            } else {
                p.to_string()
            }
        }
        None => "".to_string()
    };

    hmm.add_state(HMMState::empty_state(format!("{prefix}start")));

    hmm.add_state(HMMState::new_with_logprob(
        format!("{prefix}M0"),
        Some(match_emit_vecs[0].clone()),
        vec![format!("{prefix}start")],
        vec![0.0],
    ));

    hmm.add_state(HMMState::new_with_logprob(
        format!("{prefix}I0"),
        Some(ins_emit_vecs[0].clone()),
        vec![format!("{prefix}M0"), format!("{prefix}I0")],
        vec![match_to_ins_probs[0], ins_to_ins_probs[0]]
    ));

    let hmm_length = match_to_match_probs.len();

    for i in 1..hmm_length {
        // Add matches
        let mut match_prev_states = vec![
            format!("{prefix}M{}", i - 1),
            format!("{prefix}I{}", i - 1),
        ];

        let mut match_prev_state_probs = vec![
            match_to_match_probs[i-1],
            ins_to_match_probs[i-1],
        ];

        if i > 1 {
            match_prev_states.push(format!("{prefix}D{}", i-1));
            match_prev_state_probs.push(del_to_match_probs[i-1]);
        }

        hmm.add_state(HMMState::new_with_logprob(
            format!("{prefix}M{i}"),
            Some(match_emit_vecs[i].clone()),
            match_prev_states,
            match_prev_state_probs
        ));

        // Add Insertion State
        hmm.add_state(HMMState::new_with_logprob(
            format!("{prefix}I{i}"),
            Some(ins_emit_vecs[i].clone()),
            vec![format!("{prefix}M{i}"), format!("{prefix}I{i}")],
            vec![match_to_ins_probs[i], ins_to_ins_probs[i]],
        ));
        

        let mut del_prev_states = vec![format!("{prefix}M{}", i - 1)];
        let mut del_prev_state_probs =
            vec![match_to_del_probs[i-1]];

        if i > 1 {
            // Only add previous deletion state if i > 1
            del_prev_states.push(format!("{prefix}D{}", i - 1));
            del_prev_state_probs.push(del_to_del_probs[i-1]);
        }

        // Add deletion state
        hmm.add_state(HMMState::new_with_logprob(
            format!("{prefix}D{}", i),
            None,
            del_prev_states,
            del_prev_state_probs,
        ));
    }
    hmm.add_state(HMMState::new_with_logprob(
        format!("{prefix}end"),
        None,
        vec![
            format!("{prefix}M{}", hmm_length-1),
            format!("{prefix}I{}", hmm_length-1),
            format!("{prefix}D{}", hmm_length-1)
        ],
        // Technically this is in the hmm file, but if the vectors have been subset, the transition
        // probabilities won't add to 1
        vec![
            (1.0 - match_to_ins_probs[hmm_length-1].exp()).ln(),
            ins_to_match_probs[hmm_length-1],
            0.0
        ]
    ));
    hmm


}
 
fn read_lines(reader: &mut impl BufRead, prefix: Option<&str>, start_pos: Option<usize>, end_pos: Option<usize>) -> Result<HMM, ReaderError> {
    let mut match_to_match_probs = Vec::new();
    let mut match_to_ins_probs = Vec::new();
    let mut match_to_del_probs = Vec::new();
    let mut ins_to_match_probs = Vec::new();
    let mut ins_to_ins_probs = Vec::new();
    let mut del_to_match_probs = Vec::new();
    let mut del_to_del_probs = Vec::new();

    let mut match_emit_vecs = Vec::new();
    let mut ins_emit_vecs = Vec::new();

    // Read until HMM Start

    let mut lines = reader.lines();

    for line in lines.by_ref() {
        let line = line?;
        let first_sec = line.split_whitespace().collect::<Vec<_>>()[0];
        if first_sec == "HMM" {
            break
        }
    }
    let _s = match lines.next() {
        Some(s) => s?,
        None => return Err(ReaderError::TruncatedFile),
    };

    let mut line1: String;
    let mut line2: String;
    let mut line3: String;
    loop {
        line1 = lines.next().ok_or(ReaderError::TruncatedFile)??;

        if line1.starts_with("//") {
            break
        }
        line2 = lines.next().ok_or(ReaderError::TruncatedFile)??;
        line3 = lines.next().ok_or(ReaderError::TruncatedFile)??;

        let splits1 = line1.split_whitespace().collect::<Vec<_>>();
        let splits2 = line2.split_whitespace().collect::<Vec<_>>();
        let splits3 = line3.split_whitespace().collect::<Vec<_>>();

        if splits1[0] == "COMPO" {
            continue
        }

        // let pos: usize = splits1[0].parse()?;
        let emission_probs: Vec<f64> = splits1[1..5].iter().map(|s| s.parse().map_err(
            |e| ReaderError::ParseFloat { str: String::from(*s), source: e}
        )).collect::<Result<Vec<_>, _>>()?;
        let insertion_emission_probs: Vec<f64> = splits2[0..4].iter().map(|s| s.parse().map_err(
            |e| ReaderError::ParseFloat { str: String::from(*s), source: e}
        )).collect::<Result<Vec<_>, _>>()?;
        let transition_probs: Vec<f64> = splits3[0..7].iter().map(|s| match s {
            &"*" => Ok(f64::NEG_INFINITY),
            _ => s.parse()
        }.map_err(
            |e| ReaderError::ParseFloat { str: String::from(*s), source: e}
        )).collect::<Result<Vec<_>, _>>()?;

        match_emit_vecs.push(emission_probs);
        ins_emit_vecs.push(insertion_emission_probs);

        match_to_match_probs.push(transition_probs[0]);
        match_to_ins_probs.push(transition_probs[1]);
        match_to_del_probs.push(transition_probs[2]);
        ins_to_match_probs.push(transition_probs[3]);
        ins_to_ins_probs.push(transition_probs[4]);
        del_to_match_probs.push(transition_probs[5]);
        del_to_del_probs.push(transition_probs[6]);
        // break
    }
    // eprintln!("{}", match_to_match_probs.len());
    // panic!();
    let start = start_pos.unwrap_or(0);
    let end = end_pos.unwrap_or(match_to_match_probs.len());


    match_to_match_probs = match_to_match_probs.drain(start..end).collect();
    match_to_ins_probs = match_to_ins_probs.drain(start..end).collect();
    match_to_del_probs = match_to_del_probs.drain(start..end).collect();
    ins_to_match_probs = ins_to_match_probs.drain(start..end).collect();
    ins_to_ins_probs= ins_to_ins_probs.drain(start..end).collect();
    del_to_match_probs= del_to_match_probs.drain(start..end).collect();
    del_to_del_probs = del_to_del_probs.drain(start..end).collect();
    match_emit_vecs = match_emit_vecs.drain(start..end).collect();
    ins_emit_vecs= ins_emit_vecs.drain(start..end).collect();

    let hmm = from_vectors(
        prefix,
        match_to_match_probs,
        match_to_ins_probs,
        match_to_del_probs,
        ins_to_match_probs,
        ins_to_ins_probs,
        del_to_match_probs,
        del_to_del_probs,
        match_emit_vecs,
        ins_emit_vecs
    );
    Ok(hmm)
}

pub fn read_hmm_file(file: &Path, prefix: Option<&str>, start_pos: Option<usize>, end_pos: Option<usize>) -> Result<HMM, ReaderError> {
    let file = std::fs::File::open(file)?;
    let mut b = BufReader::new(file);
    read_lines(&mut b, prefix, start_pos, end_pos)
}


mod test {

    use super::*;
    #[test]
    fn hmm_reader() {
        let file = Path::new("test/DF000001067.hmm");
        match read_hmm_file(file, Some("Test"), None, None) {
            Ok(hmm) => {
                hmm.check_valid()
            },
            Err(e) => panic!("{e}"),
        }
    }

    #[test]
    fn hmm_reader_subset() {
        let file = Path::new("test/DF000001067.hmm");
        match read_hmm_file(file, Some("Test"), Some(5), Some(52)) {
            Ok(hmm) => {
                hmm.check_valid()
            },
            Err(e) => panic!("{e}"),
        }
    }

}
