use std::{
    fs::File, 
    io::{self, BufWriter, Write}, 
    path::Path
};

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

pub fn pprint_intervals<T: Write>(writer: &mut T, intervals: Vec<(&str, Interval)>) -> io::Result<()> {
    for (s, interval) in intervals {
        writeln!(writer, "{s}: {} - {}", interval.start, interval.stop)?
    }
    Ok(())
}

pub fn tsvprint_intervals<T: Write>(writer: &mut T, seqname: &str, intervals: Vec<(&str, Interval)>) -> io::Result<()> {
    for (s, interval) in intervals {
        writeln!(writer, "{seqname}\t{s}\t{}\t{}", interval.start, interval.stop)?
    }
    Ok(())
}

pub fn open_write(f: Option<&Path>) -> io::Result<Box<dyn Write>> {
    match f {
        Some(fname) => Ok(Box::new(BufWriter::new(File::create(fname)?))),
        None => Ok(Box::new(BufWriter::new(io::stdout())))
    }

}
