#[derive(Debug)]
pub struct Interval {
    pub start: usize,
    pub stop: usize
}
fn char_to_index(c: char) -> u8 {
    match c {
        'A' => 0,
        'C' => 1,
        'G' => 2,
        'T' => 3,
        'N' => 4,
        _ => panic!("Invalid character"),
    }
}

pub fn sequence_to_bytes(seq: &str) -> Vec<u8> {
    seq.chars().map(char_to_index).collect()
}
pub fn pprint_intervals(intervals: Vec<(&str, Interval)>) {
    for (s, interval) in intervals {
        println!("{s}: {} - {}", interval.start, interval.stop)

    }
}
