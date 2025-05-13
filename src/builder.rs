use crate::utils::*;
use crate::hmm::*;
use thiserror::Error;
use std::iter::zip;

#[derive(Error, Debug)]
pub enum HMMBuildError {
    #[error("{0} must be between 0 and 1")]
    BuildParameterError(String),
}

#[derive(Clone, Copy)]
pub struct HMMBuildSettings {
    pub match_to_match: f64,
    pub match_to_ins: f64,
    pub ins_extend: f64,
    pub del_extend: f64,
    pub loop_prob: f64,
    pub enter_skip_loop: f64,
    pub skip_to_skip: f64,
    pub match_emit_correct: f64,
}

impl HMMBuildSettings {
    pub fn new(
        match_to_match: f64,
        match_to_ins: f64,
        ins_extend: f64,
        del_extend: f64,
        loop_prob: f64,
        enter_skip_loop: f64,
        skip_to_skip: f64,
        match_emit_correct: f64,
    ) -> Result<Self, HMMBuildError> {
        if !(0.0..=1.0).contains(&match_to_match) {
            return Err(HMMBuildError::BuildParameterError(
                "match_to_match".to_string(),
            ));
        }
        if !(0.0..=1.0).contains(&match_to_ins) {
            return Err(HMMBuildError::BuildParameterError(
                "match_to_ins".to_string(),
            ));
        }
        if !(0.0..=1.0).contains(&ins_extend) {
            return Err(HMMBuildError::BuildParameterError("ins_extend".to_string()));
        }
        if !(0.0..=1.0).contains(&del_extend) {
            return Err(HMMBuildError::BuildParameterError("del_extend".to_string()));
        }
        if !(0.0..=1.0).contains(&match_emit_correct) {
            return Err(HMMBuildError::BuildParameterError(
                "match_emit_correct".to_string(),
            ));
        }
        if !(0.0..=1.0).contains(&(match_to_match + match_to_ins)) {
            return Err(HMMBuildError::BuildParameterError(
                "match_to_match + match_to_ins".to_string(),
            ));
        }
        Ok(HMMBuildSettings {
            match_to_match,
            match_to_ins,
            ins_extend,
            del_extend,
            loop_prob,
            enter_skip_loop,
            skip_to_skip,
            match_emit_correct,
        })
    }

    pub fn match_emission_vec(&self, match_index: u8) -> Vec<f64> {
        let match_emit_incorrect = (1.0 - self.match_emit_correct) / 3_f64;
        let mut emission = vec![match_emit_incorrect; 4];
        emission[match_index as usize] = self.match_emit_correct;
        emission
    }
}

impl Default for HMMBuildSettings {
    fn default() -> Self {
        HMMBuildSettings::new(
            0.9, 
            0.04,
            0.05, 
            0.1,
            0.9, 
            0.05,
            0.9, 
            0.9
        ).unwrap()
    }
}

#[allow(non_snake_case)]
pub fn create_HMM_from_motifs(motifs: Vec<&str>, motifnames: Vec<&str>, settings: &HMMBuildSettings, loop_name: &str) -> HMM {
    
    let motif_hmms = zip(motifs, motifnames)
        .map(|(s, m)| create_pHMM(
                &sequence_to_bytes(s), settings, Some(m)
        ))
        .collect::<Vec<_>>();
    // Now add skip state
    // motif_hmms.push(create_skip_state(settings, Some(loop_name)));
    let mut hmm = parallelize_HMM(motif_hmms, &format!("{loop_name}_loop"));
    loop_HMM(&mut hmm, loop_name, settings, true);
    hmm
}

#[allow(non_snake_case)]
pub fn create_pHMM(seq: &[u8], settings: &HMMBuildSettings, prefix: Option<&str>) -> HMM {
    let mut hmm = HMM::new();

    let prefix = match prefix {
        Some(p) => {
            if !p.ends_with('_') {
                format!("{}_", p)
            } else {
                p.to_string()
            }
        }
        None => "".to_string(),
    };

    let match_seed_prob =
        2.0 * (1.0 - settings.match_to_match) / (seq.len() * (seq.len() - 1)) as f64;

    // Add start state
    hmm.add_state(HMMState::empty_state(format!("{prefix}start")));

    // Add first match and insertion states
    let seq_chars = seq;
    hmm.add_state(HMMState::new(
        format!("{prefix}M0"),
        Some(settings.match_emission_vec(seq_chars[0])),
        vec![format!("{prefix}start")],
        vec![settings.match_to_match],
    ));
    hmm.add_state(HMMState::new(
        format!("{prefix}I0"),
        Some(vec![0.25, 0.25, 0.25, 0.25]),
        vec![format!("{prefix}M0"), format!("{prefix}I0")],
        vec![settings.match_to_ins, settings.ins_extend],
    ));

    for (i, c) in seq_chars.iter().enumerate().skip(1) {
        // Add match
        //
        let mismatch_prob = match_seed_prob * ((seq_chars.len() - i) as f64);
        let mut match_prev_states = vec![
            format!("{prefix}start"),
            format!("{prefix}M{}", i - 1),
            format!("{prefix}I{}", i - 1),
        ];
        let mut match_prev_state_probs = vec![
            mismatch_prob,
            settings.match_to_match,
            (1.0 - settings.ins_extend),
        ];
        if i > 1 {
            // Only add previous deletion state if i > 1
            match_prev_states.push(format!("{prefix}D{}", i - 1));
            match_prev_state_probs.push(1.0 - settings.del_extend);
        }

        hmm.add_state(HMMState::new(
            format!("{prefix}M{i}"),
            Some(settings.match_emission_vec(*c)),
            match_prev_states,
            match_prev_state_probs,
        ));

        // Add Insertion State
        hmm.add_state(HMMState::new(
            format!("{prefix}I{i}"),
            Some(vec![0.25, 0.25, 0.25, 0.25]),
            vec![format!("{prefix}M{i}"), format!("{prefix}I{i}")],
            vec![settings.match_to_ins, settings.ins_extend],
        ));

        let mut del_prev_states = vec![format!("{prefix}M{}", i - 1)];
        let mut del_prev_state_probs =
            vec![(1.0 - settings.match_to_ins - settings.match_to_match)];

        if i > 1 {
            // Only add previous deletion state if i > 1
            del_prev_states.push(format!("{prefix}D{}", i - 1));
            del_prev_state_probs.push(settings.del_extend);
        }

        // Add deletion state
        hmm.add_state(HMMState::new(
            format!("{prefix}D{}", i),
            None,
            del_prev_states,
            del_prev_state_probs,
        ));
    }
    // Add final final state
    hmm.add_state(HMMState::new(
        format!("{prefix}end"),
        None,
        vec![
            format!("{prefix}M{}", seq_chars.len() - 1),
            format!("{prefix}I{}", seq_chars.len() - 1),
            format!("{prefix}D{}", seq_chars.len() - 1),
        ],
        vec![1.0 - settings.match_to_ins, 1.0 - settings.ins_extend, 1.0],
    ));
    hmm.order_states();
    hmm
}

pub fn create_skip_state(settings: &HMMBuildSettings, prefix: Option<&str>) -> HMM {
    let mut hmm = HMM::new();

    let prefix = match prefix {
        Some(p) => {
            if !p.ends_with('_') {
                format!("{}_", p)
            } else {
                p.to_string()
            }
        }
        None => "".to_string(),
    };
    hmm.add_state(HMMState::empty_state(format!("{prefix}skip_start")));

    hmm.add_state(HMMState::new(
            format!("{prefix}skip_state"),
            Some(vec![0.25, 0.25, 0.25, 0.25]),
            vec![format!("{prefix}skip_start"), format!("{prefix}skip_state")],
            vec![1.0, settings.skip_to_skip]
    ));

    hmm.add_state(HMMState::new(
            format!("{prefix}skip_end"),
            None,
            vec![format!("{prefix}skip_state")],
            vec![1.0-settings.skip_to_skip]
    ));
    hmm.order_states();

    hmm

}

/// Creates new HMM that is made up of a sequnce of HMMs
/// * `hmms`: A vector of HMM objects, in order that they should be stitched together
#[allow(non_snake_case)]
pub fn append_HMM(hmms: Vec<HMM>) -> HMM {

    let mut all_states: Vec<&str> = Vec::new();

    
    let mut new_hmm = HMM::new();
    let mut end_state_old: &str = "";
    let mut start_state_new: &str = "";
    for (hmm_i, hmm) in hmms.iter().enumerate() {
        if hmm_i >= 1 {
            end_state_old = hmms[hmm_i - 1].get_end_states()[0];
            start_state_new = hmm.get_start_states()[0];
        }
        for state in hmm.states.iter() {
            if all_states.contains(&state.identifier.as_str()) {
                panic!("There is a repeat state name {}", state.identifier)
            }
            all_states.push(state.identifier.as_str());
            let mut new_state = state.clone();
            if state.identifier == start_state_new {
                new_state.set_transitions(vec![end_state_old.to_string()], vec![1.0]);
            }
            new_hmm.add_state(new_state);
        }
    }
    new_hmm.order_states();
    new_hmm
}

#[allow(non_snake_case)]
/// Creates new HMM that is made up of a sequence of HMMs
/// * `hmms`: A vector of HMM objects
/// * `region_prefix`: A prefix to indicate the region, for the new starts and ends
pub fn parallelize_HMM(hmms: Vec<HMM>, region_prefix: &str) -> HMM {
    let region_prefix = if region_prefix.ends_with('_') {
        region_prefix.to_string()
    } else {
        format!("{region_prefix}_")
    };

    let mut new_hmm = HMM::new();
    let mut all_states = Vec::new();
    let all_starts = hmms.iter()
        .map(|s| s.get_start_states()[0])
        .collect::<Vec<_>>();

    let all_ends = hmms.iter()
        .map(|s| s.get_end_states()[0])
        .collect::<Vec<_>>();

    let transition_prob = 1.0 / hmms.len() as f64;

    new_hmm.add_state(HMMState::empty_state(format!("{region_prefix}start")));

    for hmm in hmms.iter() {
        for state in hmm.states.iter() {
            if all_states.contains(&state.identifier.as_str()) {
                panic!("There is a repeat state name {}", state.identifier)
            }
            all_states.push(state.identifier.as_str());
            let mut new_state = state.clone();
            if all_starts.contains(&state.identifier.as_str()) {
                new_state.set_transitions(vec![format!("{region_prefix}start")], vec![transition_prob])
            }
            new_hmm.add_state(new_state);
        }
    }
    new_hmm.add_state(HMMState::new(
        format!("{region_prefix}end"),
        None,
        all_ends.iter().map(|s| s.to_string()).collect::<Vec<_>>(),
        vec![1.0; all_ends.len()]
    ));
    new_hmm.order_states();
    new_hmm
}

#[allow(non_snake_case)]
pub fn loop_HMM(hmm: &mut HMM, loop_prefix: &str, settings: &HMMBuildSettings, skip_loop: bool) {


    //TODO: Actually incorporate skip_loop
    let loop_prefix = if loop_prefix.ends_with('_') {
        loop_prefix.to_string()
    } else {
        format!("{loop_prefix}_")
    };
    // I had to change it to String instead of &str because the borrowing got complex
    let start_state_id = hmm.get_start_states()[0].to_string();
    let end_state_id = hmm.get_end_states()[0].to_string();

    let mut skip = create_skip_state(settings, Some(loop_prefix.as_str()));
    let skip_start_id = skip.get_start_states()[0].to_string();
    let skip_end_id = skip.get_end_states()[0].to_string();


    hmm.add_state(HMMState::empty_state(
        format!("{loop_prefix}start")
    ));

    hmm.add_state(HMMState::new(
        format!("{loop_prefix}end"),
        None,
        vec![end_state_id.clone()],
        vec![1.0-settings.loop_prob]
    ));
    let start_state = hmm.states.iter_mut().find(|s|
        s.identifier == start_state_id
    ).expect("Could not find start state");

    start_state.set_transitions(
        vec![
            format!("{loop_prefix}start"),
            end_state_id.clone(),
            skip_end_id.clone()
        ],
        vec![
            1.0,
            settings.loop_prob * (1.0- settings.enter_skip_loop),
            1.0
        ]
    );

    // Now add in the skip state
    let skip_start = skip.states.iter_mut().find(|s|
        s.identifier == skip_start_id
    ).expect("Could not find loop start state");

    skip_start.set_transitions(
        vec![end_state_id.clone()],
        vec![settings.loop_prob * settings.enter_skip_loop]
    );

    for state in skip.states {
        hmm.states.push(state);
    }

    hmm.order_states();
}

#[cfg(test)]
mod tests {
    use super::*;


    #[test]
    fn build_test() {
        let seq = sequence_to_bytes("AC");
        let settings = HMMBuildSettings::default();
        let hmm = create_pHMM(&seq, &settings, Some("test"));
        // panic!("States: {:?}", hmm.states);
        assert_eq!(hmm.states.len(), (3 * seq.len() - 1) + 2);
    }
    #[test]
    fn valid_hmm() {
        let seq = sequence_to_bytes("ACGT");
        let settings = HMMBuildSettings::default();
        let hmm = create_pHMM(&seq, &settings, Some("test"));

        hmm.check_valid();
    }

    #[test]
    fn skip_state_valid() {
        let settings = HMMBuildSettings::default();
        let hmm = create_skip_state(&settings, Some("test"));
        hmm.check_valid();
    }

    #[test]
    fn valid_hmm_complex() {
        let settings = HMMBuildSettings::default();

        let hmm = create_HMM_from_motifs(
            vec!["ACGTG", "GTAAG", "GAACT"],
            vec!["Rep1", "Rep2", "Rep3"],
            &settings,
            "test"
        );
        hmm.check_valid();
    }

    #[test]
    fn query_test() {
        let seq = sequence_to_bytes("ACGT");
        let settings = HMMBuildSettings::default();
        let hmm = create_pHMM(&seq, &settings, Some("test"));
        let query = sequence_to_bytes("AGTTTGT");
        let result = hmm.query(&query);
        let mut writer = std::io::stdout();
        pprint_intervals(&mut writer, result);
        // panic!();
    }

    #[test]
    fn complex_query_test() {
        let settings = HMMBuildSettings::default();

        let hmm = create_HMM_from_motifs(
            vec!["ACG", "GTA", "TCC"],
            vec!["Rep1", "Rep2", "Rep3"],
            &settings,
            "test"
        );

        let query = sequence_to_bytes("ACGACGGTAACGTCCTCCTTCC");
        let result = hmm.query(&query);
        let mut writer = std::io::stdout();
        pprint_intervals(&mut writer, result);
        // panic!();

    }

    #[test]
    fn skip_test() {
        let settings = HMMBuildSettings::default();

        let motifs = vec!["ACGTGCGAT", "GTAACGAG", "GAAGCTACT"];

        let hmm = create_HMM_from_motifs(
            motifs.clone(),
            vec!["Rep1", "Rep2", "Rep3"],
            &settings,
            "test"
        );
        let query = sequence_to_bytes(&format!("{}{}{}{}{}{}{}{}{}", motifs[0], motifs[0], "ATGATCGATTTGTAAACTACTGGGACCCTGT", motifs[0], motifs[1], motifs[2], motifs[1], motifs[2], motifs[1]));
        let result = hmm.query(&query);
        let mut writer = std::io::stdout();
        pprint_intervals(&mut writer, result);
        // panic!();

    }
}
