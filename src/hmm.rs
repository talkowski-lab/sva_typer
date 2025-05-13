use std::collections::{HashMap, HashSet};
use std::iter::zip;
use crate::utils::Interval;

#[derive(Clone, Debug)]
pub enum HmmEmission {
    NoEmit,
    Emission(Vec<f64>),
}

#[derive(Clone, Debug)]
pub struct HMMState {
    pub identifier: String,
    emission: HmmEmission,
    prev_states: Vec<String>,
    prev_state_transitions: Vec<f64>,
}

impl HMMState {
    pub fn new(
        identifier: String,
        emission: Option<Vec<f64>>,
        prev_states: Vec<String>,
        prev_state_transitions: Vec<f64>,
    ) -> Self {
        HMMState {
            identifier,
            emission: match emission {
                None => HmmEmission::NoEmit,
                Some(v) => HmmEmission::Emission(v.iter().map(|p| p.ln()).collect()),
            },
            prev_states,
            prev_state_transitions: prev_state_transitions.iter().map(|p| p.ln()).collect(),
        }
    }
    pub fn empty_state(identifier: String) -> Self {
        HMMState {
            identifier,
            emission: HmmEmission::NoEmit,
            prev_states: Vec::new(),
            prev_state_transitions: Vec::new(),
        }
    }

    pub fn set_transitions(&mut self, prev_states: Vec<String>, prev_state_transitions: Vec<f64>) {
        self.prev_states = prev_states;
        self.prev_state_transitions = prev_state_transitions.iter().map(|p| p.ln()).collect();
    }
}

pub struct HMM {
    pub states: Vec<HMMState>,
}

impl HMM {
    pub fn new() -> Self {
        HMM { states: Vec::new() }
    }

    pub fn add_state(&mut self, state: HMMState) {
        self.states.push(state);
    }

    pub fn get_index_map(&self) -> HashMap<String, usize> {
        let mut map = HashMap::new();
        for (i, state) in self.states.iter().enumerate() {
            map.insert(state.identifier.clone(), i);
        }
        map
    }

    // TODO: Change to result output
    pub fn check_valid(&self) {
        // Check that there is only one start
        let start_states = self.get_start_states();
        if start_states.len() != 1 {
            panic!(
                "There should be exactly one start state, found {}",
                start_states.len()
            );
        }
        let end_states = self.get_end_states();
        if end_states.len() != 1 {
            panic!(
                "There should be exactly one end state, found {}",
                end_states.len()
            );
        }

        // First check that all states are constructed correctly
        for state in self.states.iter() {
            if state.prev_states.len() != state.prev_state_transitions.len() {
                panic!("State {} does not have the same number of previous states as previous state transitions", state.identifier);
            }
        }
        let all_identifiers: Vec<String> =
            self.states.iter().map(|s| s.identifier.clone()).collect();
        let all_identifier_set: HashSet<String> = all_identifiers.iter().cloned().collect();

        if all_identifiers.len() != all_identifier_set.len() {
            panic!("State identifiers are not unique");
        }

        self.states.iter().for_each(|state| {
            for prev_state in state.prev_states.iter() {
                if !all_identifier_set.contains(prev_state) {
                    panic!(
                        "State {} has a previous state {} that does not exist",
                        state.identifier, prev_state
                    );
                }
            }
        });

        self.check_valid_emissions();
        self.check_valid_transitions();
    }

    /// Properly orders the nonemitting states of the HMM so that their viterbi scores are
    /// calculated after their previous states
    ///
    /// Make sure this is called after creating an HMM or the probabilities might get messed up
    pub fn order_states(&mut self) {
        let pos_map = self.get_index_map();
        let copy = self.states.clone();
        let mut normal_states: Vec<String> = Vec::new();
        let mut silent_states: Vec<String> = Vec::new();

        for state in self.states.iter() {
            match state.emission {
                HmmEmission::NoEmit => silent_states.push(state.identifier.clone()),
                HmmEmission::Emission(_) => normal_states.push(state.identifier.clone()),
            }
        }

        //
        let mut sorted = Vec::new();

        while !silent_states.is_empty() {
            let mut unused = Vec::new();

            for state in &silent_states {
                let has_incoming_silent = self.states[pos_map[state]]
                    .prev_states
                    .iter()
                    .any(|s| silent_states.contains(s));
                if !has_incoming_silent {
                    sorted.push(state.to_string())
                } else {
                    unused.push(state.to_string());
                }
            }

            assert!(unused.len() < silent_states.len());
            silent_states = unused;
        }
        normal_states.extend(sorted);

        self.states.clear();
        for state in normal_states.iter() {
            self.states.push(copy[pos_map[state]].clone())
        }
    }

    /// 
    /// States MUST be in order using self.order_states() beforehand. Otherwise the log
    /// probabilities for the nonemitting states will get messed up. This is so self can be passed
    /// as a immutable reference making it easier for parallelization
    /// * `query`: 
    pub fn query(&self, query: &[u8]) -> Vec<(&str, Interval)> {
    // pub fn query(&mut self, query: &[u8]) -> (Vec<&str>, Vec<usize>) {
        let start_state = self.get_start_states()[0];
        let end_state = self.get_end_states()[0];
        let (_lp_mat, trace_mat) = self.gen_viterbi_mats(query, start_state);
        // panic!("Trace mat: {:?}", _lp_mat);
        let index_map = self.get_index_map();
        let (path, query_indexes) = self.traceback(&trace_mat, &index_map, start_state, end_state);

        // (path, query_indexes)

        convert_to_intervals(path, query_indexes)
    }

    fn gen_viterbi_mats(
        &self,
        query: &[u8],
        start_state: &str,
    ) -> (Vec<Vec<f64>>, Vec<Vec<Option<&str>>>) {
        let mut lp_mat = vec![vec![f64::NEG_INFINITY; query.len() + 1]; self.states.len()];
        let mut trace_mat = vec![vec![None; query.len() + 1]; self.states.len()];
        let map = self.get_index_map();
        lp_mat[map[start_state]][0] = 0.0;
        for query_i in 0..(query.len() + 1) {
            for (state_i, state) in self.states.iter().enumerate() {
                if let Some((score, trace_state)) =
                    self.viterbi_score(query, query_i, state, &lp_mat, &map)
                {
                    lp_mat[state_i][query_i] = score;
                    trace_mat[state_i][query_i] = Some(trace_state);
                }
            }
        }
        (lp_mat, trace_mat)
    }
    fn viterbi_score<'a>(
        &self,
        query: &[u8],
        index: usize,
        state: &'a HMMState,
        lp_mat: &[Vec<f64>],
        index_map: &HashMap<String, usize>,
    ) -> Option<(f64, &'a str)> {
        let (ln_em, traceback): (f64, usize) = match &state.emission {
            HmmEmission::NoEmit => (0.0, 0),
            HmmEmission::Emission(emit_probs) => {
                if index == 0 {
                    // index = 0 means before first character
                    return None;
                }
                let base = query[index - 1];
                match base {
                    4 => (0.0, 1), // 4 corresponds to N
                    _ => (emit_probs[base as usize], 1),
                }
            }
        };
        let mut best_state = None;
        let mut best_score = f64::NEG_INFINITY;

        for (i, prev_state) in state.prev_states.iter().enumerate() {
            let prev_state_i = index_map[prev_state];
            let prev_score = lp_mat[prev_state_i][index - traceback];
            let trans_lp = state.prev_state_transitions[i];
            let score = prev_score + trans_lp + ln_em;
            if score > best_score {
                best_score = score;
                best_state = Some(prev_state);
            }
        }
        best_state.map(|s| (best_score, s.as_str()))
    }

    fn traceback<'a>(
        &self,
        trace_mat: &[Vec<Option<&'a str>>],
        index_map: &HashMap<String, usize>,
        start_state: &'a str,
        end_state: &'a str
    ) -> (Vec<&'a str>, Vec<usize>) {
        let mut state = end_state;
        let mut state_i = index_map[state];
        let mut index = trace_mat[0].len() - 1;
        let mut path = Vec::new();
        let mut query_indexes = Vec::new();


        while state != start_state {
            path.push(state);
            query_indexes.push(index);
            let prev_state = trace_mat[state_i][index]
                .unwrap_or_else(|| panic!("No previous state found for {}", state));
            if let HmmEmission::Emission(_) = self.states[state_i].emission {
                index -= 1;
            }
            state = prev_state;
            state_i = index_map[state];
        }
        path.push(start_state);
        query_indexes.push(index);
        path.reverse();
        query_indexes.reverse();
        (path, query_indexes)
    }

    /// Checks if all states have valid emissions
    fn check_valid_emissions(&self) {
        let mut sum;
        for state in self.states.iter() {
            if let HmmEmission::Emission(d) = &state.emission {
                sum = (d.iter().map(|p| p.exp()).sum::<f64>() * 1000.0).round() / 1000.0;
                if sum != 1.0 {
                    panic!(
                        "Emission probabilities for state {} sum to {}, not 1",
                        state.identifier, sum
                    );
                }
            }
        }
    }

    /// Checks if all transition probabilities sum to 1
    fn check_valid_transitions(&self) {
        let prob_matrix = self.get_transition_matrix();
        let mut sum;
        let end_state = self.get_end_states()[0];
        for (i, row) in prob_matrix.iter().enumerate() {
            //TODO: This is kinda lazy, change to a better check for end states
            if end_state == self.states[i].identifier {
                continue;
            }
            sum = (row.iter().map(|p| p.exp()).sum::<f64>() * 1000.0).round() / 1000.0;
            if sum != 1.0 {
                panic!(
                    "Transition probabilities for state {} sum to {}, not 1",
                    self.states[i].identifier, sum
                );
            }
        }
    }
    pub fn get_start_states(&self) -> Vec<&str> {
        self.states
            .iter()
            .filter(|s| s.prev_states.is_empty())
            .map(|s| s.identifier.as_str())
            .collect::<Vec<&str>>()
    }

    pub fn get_end_states(&self) -> Vec<&str> {
        let mut states_out = HashSet::new();

        for state in self.states.iter() {
            for prev_state in state.prev_states.iter() {
                states_out.insert(prev_state);
            }
        }

        let end_states = self
            .states
            .iter()
            .filter(|s| !states_out.contains(&s.identifier))
            .map(|s| s.identifier.as_str())
            .collect::<Vec<&str>>();

        end_states
    }
    pub fn get_transition_matrix(&self) -> Vec<Vec<f64>> {
        let map = self.get_index_map();
        // panic!("{:?}", map);
        let mut matrix = vec![vec![f64::NEG_INFINITY; self.states.len()]; self.states.len()];
        for (i, state) in self.states.iter().enumerate() {
            for (prev_state, prob) in zip(&state.prev_states, &state.prev_state_transitions) {
                let x = map[prev_state];
                matrix[x][i] = *prob;
            }
        }
        matrix
    }
}

impl Default for HMM {
    fn default() -> Self {
        Self::new()
    }
}

/// Outputs 0-based intervals for the motif starts and stops
/// FYI: This functional heavily assumes the naming schema for states used throughout the building process
/// * `state_names`: 
/// * `state_pos`: 
fn convert_to_intervals(state_names: Vec<&str>, state_pos: Vec<usize>) -> Vec<(&str, Interval)> {
    // This removes loop states by just assuming that "loop_start" and "loop_end" are in the names,
    // which might not always be true, just fyi
    let mut intervals = vec![];

    for (s, i) in zip(
        state_names, 
        state_pos.iter().map(|i| match i { 
            // There's some funkiness because we had to add a fake start state to make the
            // traceback work, but 0 and 1 are the same position, and everything else is minus 1
            0 => 0,
            _ => i-1
        })
    )  {
            if s.contains("_start") && !s.contains("_loop_start") {
                intervals.push((s.strip_suffix("_start").unwrap(), Interval {start: i, stop: usize::MAX}))
            }
            if s.contains("_end") && !s.contains("_loop_end") {
                for (n, interval) in intervals.iter_mut() {
                    if n == &s.strip_suffix("_end").unwrap() && interval.stop == usize::MAX {
                        interval.stop = i
                    }
                }
            }
    }
    intervals
}
