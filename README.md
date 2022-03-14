# README - In silico read generator

Includes two scripts, ReadGen.py and align_seq.py.

- **ReadGen.py** takes a reference DNA sequence in fasta format and generates
randomly produced reads in a FASTQ format. The program can generate single-end
or paired-end reads.

- **align_seq.py** aligns the generated reads to a reference sequence.

## Installation
The scripts are written in python, and use the numpy package.
- Python v3.8.8
- Numpy v1.20.1

## Reads generator - ReadGen.py
### Usage
The scripts are run through the command line.

 ```
 ReadsGen.py [-h] -i REFSEQ [--length READLEN] [--coverage COVERAGE] [--pair-end] [-o OUT]
 ```

 Flag options:
   -h, --help:           show help message and exit
   -i:           Input reference sequence in fasta format. Reads will be generated from this sequence.
   --length:     Length of the reads to be generated. Default is 100.
   --coverage:  Sequencing depth or coverage. Default is 10.
   --pair-end:       Generate paired-end reads. If this argument is not included, reads will be single-end.
   -o:               Basename for output files. Default is 'reads' for single-end, 'reads_1' and 'reads_2' for paired-end.

 Example usage:
```
 python ReadsGen.py -i rCRS.fa --pair-end --coverage 20 --readlen 150 -o reads
 ```
 This command generates paired-end reads of length 150, stored in two outputfiles
that will be named reads_1.fastq and read_2.fastq. Since the specified coverage
is 20, the sequence will be 'read' 20 times.

## Output
The outputs are FASTQ files; for the single-end reads one file is generated and
for paired-end each group of reads is printed to a different file, so 2 output
files numbered 1 and 2 are generated.

## Aligner - align_seq.py
This script splits the reads into k-mers and searches for their position in
the reference sequence, then reconstructs the sequence from the found positions.
Once a kmer matches the sequence, the rest of the read is compared so there
can be some mismatches on the reads, but at least one of the k-mers must match
for the read to be aligned.

### Usage

 ```
 align_reads.py [-h] -ref REF_SEQ [-i READS] [-p1 READS1] [-p2 READS2] [-o OUT] [-o OUT]
 ```
 optional arguments:
   -h, --help:   show help message and exit
   -ref:  Input reference sequence in fasta format
   -i:      Input reads in FASTQ format. 1 file, for single-end reads.
   -p1:    FASTQ file 1 for paired end reads
   -p2:   FASTQ file 2 for paired end reads
   --k: Length of the k-mers
   -o:        Basename for assembly file. Default is 'assembly'

 Example usage:
```
 python align_reads.py -ref rCRS.fa -p1 reads1.fastq -p2 reads2.fastq
 ```
This command takes two files containing paired-end reads and will align them
to the reference sequence.

## Output
This script generates two files:

- 'assembly' file is a fasta file that contains one sequence obtained from the
reads, according to their positions in the reference sequence. Missing positions
are written as 'N', and the sequence will have the same length as the reference.

- 'contigs' is a fasta file that contains the sequences that were mapped
without missing positions between them.

Also a percentage identity between the generated sequence and the can be printed to the screen if it specified in the
command.
