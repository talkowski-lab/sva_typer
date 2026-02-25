use std::{
    fs::File, 
    io::{self, BufWriter, Write}, 
    path::Path,
    iter::zip
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

pub fn write_header(writer: &mut impl Write, write_hmm_state: bool, write_query_seq: bool) -> io::Result<()> {
    match write_hmm_state {
        true => writeln!(writer, "ID\tstate\tquery_i\tquery_base"),
        false => match write_query_seq {
            true => writeln!(writer, "ID\tregion\tstart\tend\tseq"),
            false => writeln!(writer, "ID\tregion\tstart\tend")
        }
    }
}

pub fn tsvprint_intervals(writer: &mut impl Write, seqname: &str, intervals: Vec<(&str, Interval)>) -> io::Result<()> {
    for (s, interval) in intervals {
        writeln!(writer, "{seqname}\t{s}\t{}\t{}", interval.start, interval.stop)?
    }
    Ok(())
}

pub fn tsvprint_intervals_withseq(writer: &mut impl Write, seqname: &str, query: &str, intervals: Vec<(&str, Interval)>) -> io::Result<()> {
    for (s, interval) in intervals {
        writeln!(writer, "{seqname}\t{s}\t{}\t{}\t{}", interval.start, interval.stop, &query[interval.start..interval.stop])?
    }
    Ok(())
}

pub fn tsvprint_hmmstates(writer: &mut impl Write, seqname: &str, query: &str, state_names: Vec<&str>, state_pos: Vec<usize>) -> io::Result<()> {
    for (s, i) in zip(
        state_names, 
        state_pos.iter().map(|i| match i { 
            // There's some funkiness because we had to add a fake start state to make the
            // traceback work, but 0 and 1 are the same position, and everything else is minus 1
            0 => 0,
            _ => i-1
        })
    )  {
            writeln!(writer, "{seqname}\t{s}\t{i}\t{}", query.chars().nth(i).unwrap())?
    }
    Ok(())


}

pub fn open_write(f: Option<&Path>) -> io::Result<Box<dyn Write>> {
    match f {
        Some(fname) => Ok(Box::new(BufWriter::new(File::create(fname)?))),
        None => Ok(Box::new(BufWriter::new(io::stdout())))
    }

}
