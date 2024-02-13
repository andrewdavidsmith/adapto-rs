/* MIT License
 *
 * Copyright (c) 2023-2024 Andrew Smith
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use, copy,
 * modify, merge, publish, distribute, sublicense, and/or sell copies
 * of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

use rayon::prelude::*;
use std::cmp::{max, min};
use std::error::Error;
use std::io::{Read, Write};
use std::ptr;

// the rust_htslib crate is not ideal for our purpose
use rust_htslib::bgzf;
use rust_htslib::bgzf::CompressionLevel as CompLvl;
use rust_htslib::tpool::ThreadPool;

/// Just the naive algorithm for string matching with bounded
/// mismatches.
fn naive_matching(adaptor: &[u8], read: &[u8], min_frac: f64, min_ltrs: usize) -> usize {
    let n = adaptor.len();
    let m = read.len();
    let d_delta = 1f64 - min_frac;

    let i_lim1 = if n > m { 0 } else { m + 1 - n };
    // let min_match = (min_frac * n as f64).ceil() as usize;
    let mut d_max = n as f64 * d_delta;
    for i in 0..i_lim1 {
        let mut d: usize = 0;
        let mut j: usize = 0;
        while d <= d_max as usize && j < n {
            if read[i + j] != adaptor[j] {
                d += 1;
            }
            j += 1;
        }
        if d <= d_max as usize {
            return i;
        }
    }

    let i_lim2 = m + 1 - min_ltrs;
    let mut j_lim = m - i_lim1;
    for i in i_lim1..i_lim2 {
        d_max -= d_delta;
        let mut d: usize = 0;
        let mut j: usize = 0;
        while d <= d_max as usize && j < j_lim {
            if read[i + j] != adaptor[j] {
                d += 1;
            }
            j += 1;
        }
        if d <= d_max as usize {
            return i;
        }
        j_lim -= 1;
    }
    m
}

/// Find the positions in the read of the first non-N and last non-N.
fn trim_n_ends(read: &[u8]) -> (usize, usize) {
    (
        match read.iter().position(|&x| x != b'N') {
            Some(x) => x,
            _ => 0,
        },
        match read.iter().rposition(|&x| x != b'N') {
            Some(x) => x + 1,
            _ => 0,
        },
    )
}

/// Find the positions in the read where quality scores indicate the
/// read should be trimmed. This is copied from cutadapt source.
fn qual_trim(qual: &[u8], cut_front: i32, cut_back: i32) -> (usize, usize) {
    const QUAL_BASE: i32 = 33; // assumes base quality starts at 33

    /* ADS: COPIED FROM cutadapt SOURCE */
    let n = qual.len();

    //  find trim position for 5' end
    let mut start: usize = 0;
    let mut s: i32 = 0;
    let mut max_qual: i32 = 0;

    if cut_front > 0 {
        let cut_front = cut_front + QUAL_BASE;
        for i in 0..n {
            s += (cut_front + QUAL_BASE) - qual[i] as i32;
            if s < 0 {
                break;
            }
            if s > max_qual {
                max_qual = s;
                start = i + 1;
            }
        }
    }
    // same for 3' end
    let mut stop: usize = n;
    max_qual = 0;
    s = 0;
    let cut_back = cut_back + QUAL_BASE;
    for i in (0..n).rev() {
        s += cut_back - qual[i] as i32;
        if s < 0 {
            break;
        }
        if s > max_qual {
            max_qual = s;
            stop = i;
        }
    }
    if start >= stop {
        (start, stop) = (0, 0)
    }
    (start as usize, stop as usize)
}

fn shift(buf: &mut [u8], cursor: &mut usize, filled: &mut usize) {
    let mut j = 0;
    for i in *cursor..*filled {
        buf[j] = buf[i];
        j += 1;
    }
    *filled = j;
    *cursor = 0;
}

/// FQRec is a FASTQ record that represents the position of the start
/// of the name (n), the start of the read sequence (r), the start of
/// the other name, the one with the "+" (o), and the start of the
/// quality scores (q). The end of the record (e) is the start of the
/// next record, but having it in each FQRec allows these to be
/// processed by independent threads.
#[derive(Default)]
struct FQRec {
    n: usize, // start of "name"
    r: usize, // start of "read"
    o: usize, // start of "other"
    q: usize, // start of "quality" scores
    e: usize, // end of the record
}

impl std::fmt::Display for FQRec {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "@={}, seq={}, +={}, qual={}, end={}",
            self.n, self.r, self.o, self.q, self.e
        )
    }
}

impl FQRec {
    fn process(
        &mut self,
        adaptor: &[u8],
        cutoff: u8,
        min_frac: f64,
        min_ltrs: usize,
        buf: &Vec<u8>,
    ) {
        let seqlen = if self.r < self.o {
            self.o - self.r - 1
        } else {
            0
        };
        let (qstart, qstop) = qual_trim(&buf[self.q..self.q + seqlen], 0, cutoff as i32);
        // consecutive N values at both ends
        let (nstart, nstop) = trim_n_ends(&buf[self.r..self.r + seqlen]);
        // so no N or low qual bases can interfere with adaptor
        let mut stop = min(qstop, nstop);

        // find the adaptor at the 3' end
        let adaptor_start =
            naive_matching(adaptor, &buf[self.r..self.r + stop], min_frac, min_ltrs);

        stop = min(stop, adaptor_start);
        let (_, nstop) = trim_n_ends(&buf[self.r..self.r + stop]);
        stop = min(stop, nstop);
        let start = min(max(qstart, nstart), stop);

        /* ADS: Removing the comments in the next two lines breaks up
         * this function, which would allow the work to be done in two
         * loops, but that would mean waiting for slower threads. */

        // }
        // fn compress(&mut self, buf: &Vec<u8>) {

        /* ADS: below here, the instance variables other than n and e
         * become invalidated
         */

        let b = buf.as_ptr() as *mut u8;
        let r_sz = stop - start;

        // Set cursor just at the end of the read name (first
        // space/tab/newline)
        let mut cursor = self.n
            + match &buf[self.n..self.r].iter().position(|&x| x == b' ') {
                Some(x) => x,
                _ => &self.r,
            };
        unsafe {
            *b.add(cursor) = b'\n';
        }
        cursor += 1;
        unsafe {
            ptr::copy(b.add(self.r + start), b.add(cursor), r_sz);
        }
        cursor += r_sz;
        unsafe {
            *b.add(cursor) = b'\n';
        }
        cursor += 1;
        unsafe {
            *b.add(cursor) = b'+';
        }
        cursor += 1;
        unsafe {
            *b.add(cursor) = b'\n';
        }
        cursor += 1;
        /* ADS: the code above removes the "header" after the '+'
         * symbol on lines numbered 2 (mod 3)
         */
        unsafe {
            ptr::copy(b.add(self.q + start), b.add(cursor), r_sz);
        }
        cursor += r_sz;
        unsafe {
            *b.add(cursor) = b'\n';
        }
        cursor += 1;
        self.e = cursor;
        // postcondition of this function: self.n and self.e are valid
    }
    fn write<W: Write>(&self, buf: &Vec<u8>, writer: &mut W) {
        writer.write(&buf[self.n..self.e]).unwrap();
    }
}

#[inline(always)]
fn next_line(buf: &mut [u8], filled: usize, offset: usize) -> usize {
    for i in offset..filled {
        if buf[i] == b'\n' {
            return i + 1;
        }
    }
    usize::MAX
}

#[inline(always)]
fn get_next_record(buf: &mut [u8], cursor: &mut usize, filled: usize) -> FQRec {
    // ADS: here is where we should detect malformed records
    let n = *cursor;
    let r = next_line(buf, filled, n);
    let o = next_line(buf, filled, r);
    let q = next_line(buf, filled, o);
    let e = next_line(buf, filled, q);
    if e != usize::MAX {
        *cursor = e;
        debug_assert!(buf[n] == b'@');
    }
    FQRec { n, r, o, q, e }
}

fn process_reads<R: Read, W: Write>(
    buffer_size: usize,
    adaptor: &[u8],
    reader: &mut R,
    mut writer: &mut W,
    cutoff: u8,
    min_frac: f64,
    min_ltrs: usize,
) -> Result<(), Box<dyn Error>> {
    let mut buf: Vec<u8> = vec![b'\0'; buffer_size];
    let mut filled = 0usize;
    let mut cursor = 0usize;

    let mut recs: Vec<FQRec> = Vec::new();

    loop {
        // move any unused data to start of buffer
        shift(&mut buf, &mut cursor, &mut filled);

        // read the input to fill the buffer
        filled += reader.read(&mut buf[filled..])?;

        // find the sequenced read records
        recs.clear(); // keep capacity
        loop {
            let fq = get_next_record(&mut buf, &mut cursor, filled);
            if fq.e == usize::MAX {
                break;
            }
            recs.push(fq);
        }

        // find end-points of trimmed reads
        recs.par_iter_mut()
            .for_each(|fq_rec| fq_rec.process(&adaptor, cutoff, min_frac, min_ltrs, &buf));

        /* ADS: could do separately: make record a contiguous chunk */
        // recs.iter_mut().for_each(|x| x.compress(&buf));

        // write all records to output file
        recs.iter_mut().for_each(|x| x.write(&mut buf, &mut writer));

        // exit if previous read hit end of file
        if filled < buf.len() {
            break;
        }
    }

    Ok(())
}

pub fn remove_adaptors(
    zip: bool,
    n_threads: u32,
    buf_sz: usize,
    adaptor: &[u8],
    input: &String,
    output: &String,
    cutoff: u8,
    min_frac: f64,
    min_ltrs: usize,
) -> Result<(), Box<dyn Error>> {
    let lvl = match zip {
        true => CompLvl::Default,
        false => CompLvl::NoCompression,
    };
    let mut reader = bgzf::Reader::from_path(input)?;
    let mut writer = bgzf::Writer::from_path_with_level(output, lvl)?;

    let tpool = ThreadPool::new(n_threads - 1)?;
    if n_threads > 1 {
        reader.set_thread_pool(&tpool)?;
        writer.set_thread_pool(&tpool)?;
    }
    process_reads(
        buf_sz,
        adaptor,
        &mut reader,
        &mut writer,
        cutoff,
        min_frac,
        min_ltrs,
    )
}
