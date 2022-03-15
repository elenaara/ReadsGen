# -*- coding: utf-8 -*-
"""
-- In silico read generator --

name: ReadsGen.py

Created on Mon Mar 7 11:16:57 2022

@author: Elena Aramendía

This script generates reads from a reference sequence, with 
desired read length and coverage. Output is a FASTQ file.


Procedure:
    1. Input file and parameters are read. 
    Input ssequence can be randomly mutated if specified.
    2. Reads are generated by splitting the sequence:
        
        - For single-end reads: number of reads according to readlength are 
        calculated, and a random number of positions in the sequence are 
        generated. From each of this ranodm positions bases are read until
        reaching the read length, and then stored in a list as a read. 
        This procedure (splitting and reading all the positions) is repeated 
        as many times as specified coverage. For each reading, start and end
        of the sequence are fixed.
        
        -  For paired-end reads: similar procdure as single-end, but half of 
        the positions are calculated and for each position 2 reads are taken:
        one starting at that position and 1 starting 10 positions ahead. 
        The second read is converted into the reverse complement.
        Each group of reads is stored in a separate list.
        
    3. Quality line is generated: for each read a quality line of the same 
    length is generated. The quality line can show 'perfect' scores (all bases
    show a Q40 score or  0.0001 error rate), or, if the --q flag is included,
    more 'realistic' scores are generated. For the first 90% length of the read
    error values are taken from an uniform distribution taking values from 
    0.0001 to 0.001, based on error rates generated from the main Illumina 
    sequencing platforms.
    For the end of the read quality scores tend to drop, so the values for the
    last 10% of the read are generated from a uniform distribution with 
    slightly higher values (from 0.0004 to 0.0015).
    This values are converted to ASCII characters by calling the qscore function.

    4. Print sequence to output files. For each sequence 4 lines are printed,
    following the FASTQ format:
        header, nucleotide sequence, spacer line and quality line.
    For the single-end reads one file is generated and for paired-end each 
    group of reads is printed to a different file, so 2 output files numbered
    1 and 2 are generated.
    
Functions:
    fasta1line:
        Takes fasta file and returns a dictionary, used to get the whole 
        sequence in one line in case it is in multiple lines.
        
    qscore:
        Takes an error rate between 0 and 1 and converts it to phred score in ASCII.
    rev_comp:
        Takes a DNA sequence and generates the reverse complement sequence.
Modules:
    argparse
    random
    numpy
    
Input:
    - ref sequence: required, sequence to get reads from
    - read length: lenght of the reads, if not specified 100 is default
    - coverage: reading depth or coverage, if not specified 10 is default
    - pair-end: if this argument is included the script generates paired-end reads
    If not, default is single-end
    
Output:
    FASTQ file (2 files for paired-end reads)
    
Example usage:
    
    python ReadsGen.py -i rCRS.fa
    
    - With optional arguments:
        
    python ReadsGen.py -i rCRS.fa --coverage 20 --readlen 150 -o reads
    
"""
# %% Modules
import argparse
import random
import numpy as np

# %% Functions

def fasta1line(fasta_file):
    """
    Function to make single-line sequences from a fasta file, in case they are in multiple lines
    
    - Input: fasta_file
        A fasta file with headers starting with '>'
        If the first line it's not a header the funciton will give an error.
    
    - Return: fasta_dict
        Dictionary with header as key and sequence as value.
    
    """
    # Define lists of sequences, sequence fragments and headers
    seq = []
    seq_fragments = []
    header_list = []
    fasta_dict = {}
    index = 0 # counter
    fasta = False # Check that is a valid FASTA file
    for line in fasta_file:
    # Search for the header line
        if line.startswith('>'):
            fasta = True
            header = line.rstrip() # Remove newline character
            header_list.append(header)
            #If this is NOT the first sequence
            if seq_fragments:
                #Add sequence to a list of sequences
                curr_seq = ''.join(seq_fragments)
                seq.append(curr_seq)
            seq_fragments = []
        else:
            # Found more of existing sequence
            curr_seq = line.rstrip() # remove new line character
            seq_fragments.append(curr_seq) # Add fragment to current sequence   
    if fasta == False:
       raise Exception('Not a valid format. Please try again with a FASTA file.')
    # Appending the last fasta sequence 
    if seq_fragments:
        #if the file is not empty
        curr_seq = ''.join(seq_fragments)
        seq.append(curr_seq)
    for i in header_list: # Create dictionary with all the sequences
        fasta_dict[i] = seq[index]
        index += 1
    return fasta_dict

def rev_comp(sequence):
    """
    Function to get the reverse complement to a DNA sequence

    Parameters
    ----------
    sequence : str
        DNA sequence to reverse.

    Returns
    -------
    str
        Returns reverse complement of the sequence (ACs are turned into TGs and viceversa, and the sequence is backwards).

    """
    reverse = {'A':'T','T':'A','C':'G','G':'C','N':'N'}
    reverse_seq = ''
    for i in sequence.upper():
        reverse_seq = reverse_seq + reverse[i]
    return reverse_seq[::-1]

def mutate(position, sequence):
    """
    Function to randomly mutate a desired position of a DNA sequence.
    All possible nucleotides have the same probability.

    Parameters
    ----------
    sequence: str
        Sequence to mutate.
        
    position : int
        Position in the sequence to mutate.
    
    Returns
    -------
    mut_seq: str
        Mutated sequence.

    """
    Nts = ['A','C','T','G']
    new = random.choice(Nts)
    old = sequence[position] 
    while new == old: # If the chosen nucleotide is the same as the old one
    # Generate another one until they are different
        new = random.choice(Nts)
    mut_seq = seq[0:position] + new + seq[position+1:]
    return(mut_seq)

def qscore(error):
    """
    Function to turn error rate value into Phred score in ASCII

    Parameters
    ----------
    error : float
        Error rate, number between 0 and 1.

    Returns
    -------
    int
        ASCII character representing the input error rate.

    """
    #  Q = -10log10(e)
    import math
    Q = -10*math.log(error,10)
    asc = int(Q) + 33
    return chr(asc)
# %% Argument parser

parser = argparse.ArgumentParser(description='This program takes a reference sequence and generates reads from this sequence with a desired length and coverage.',
                                 epilog='Example usage: python ReadsGen.py -i rCRS.fa --pair-end --coverage 20 --readlen 150 -o reads')


parser.add_argument(                                                  # fasta file
    '-i',
    type = argparse.FileType('r'),
    dest = 'refseq',
    required = True,
    help = 'Input reference sequence in fasta format'
    )

parser.add_argument(                                                  # length of the reads 
    '--length',
    type = int,
    dest = 'readlen',
    required = False,
    default = 100,
    help = 'Length of the reads'
    )

parser.add_argument(                                                  # sequencing depth or coverage
    '--coverage',
    type = float,
    dest = 'coverage',
    required = False,
    default = 10,
    help = 'Sequencing depth (how many times a given nucleotide is included in unique reads)'
    )

parser.add_argument(                                                  # paired reads
    '--pair-end',
    action='store_true',
    dest = 'paired',
    help = 'Generate paired-end reads'
    )

parser.add_argument(                                                  # quality
    '--q',
    action='store_true',
    dest = 'q',
    help = 'Generate different qualities. By default all positions show a quality of I (Q40)'
    )

parser.add_argument(                                                  # quality
    '--mut',
    type= int,
    dest = 'mut',
    required = False,
    default = 0,
    help = 'Randomly mutate the specified number of positions'
    )

parser.add_argument(                                                  # output 
    '-o',
    type = str,
    dest = 'out',
    required = False,
    default = 'reads', ## CHECK THIS
    help = 'Basename for output files'
    )


args = parser.parse_args() 

#%%

# Import reference sequence 

# Reference seq
fasta = fasta1line(args.refseq)
# Get header and seq
for key in fasta:
    header = key[1:]
    seq = fasta[key]

# Mutate reference sequence if specified
if args.mut > 0:
    mut_positions = [int(random.uniform(0,len(seq))) for i in range(args.mut)]
    for pos in mut_positions:
        seq = mutate(pos,seq)

# %% Generate reads

# Read Parameters
readlen = args.readlen
coverage = args.coverage
paired = args.paired


# %% SINGLE END READS
# Split reference seq 
if args.paired == False:
    reads = [] # list for generated reads
    count = 0 # counter for the number of times the sequence is "read", as many as desired coverage
    while count <= coverage:    
        
        # Number of reads
        Nreads = int(len(seq)/readlen)
        # Generate splitting positions
        #positions = np.random.randint(0,len(seq),Nreads)
        positions = [int(random.uniform(0,len(seq))) for i in range(Nreads)]
        positions.append(0)
        positions.append(len(seq)-100)
        for i in positions: # Start reading in each of the positions
        # Read from position i until reaching read length
            read = [seq[i:j] for j in range(i, len(seq) + 1) if len(seq[i:j]) == readlen]
            # Add current read to the reads list
            reads += read 
                
                
        count += 1 # add counter
    
# Quality line
    if args.q:
        # generate 'realistic' qualities
        qualities = []
        for read in reads:       
            quality_values = 0.9*readlen 
            # Median error = 0.1026% -> 0.0010
            randQ = random.uniform(0.0004,0.001)
            error_rates1 = np.random.uniform(0.0001,randQ,size=int(quality_values)) 
            #statistics.median(error_rates1)
            error_rates2 = np.random.uniform(randQ,0.0015,size=int(readlen - quality_values))
            #statistics.median(error_rates2)
            error_ascii1 = ''
            error_ascii2 = ''
            #  Q = -10log10(e)
            # 
            for i in error_rates1:
                Q = qscore(i)
                error_ascii1 += Q
            for i in error_rates2:
                Q = qscore(i)
                error_ascii2 += Q
            qual = error_ascii1 + error_ascii2
            qualities.append(qual)
    else:
        qualities = []
        qual = readlen*'I'
        for i in range(len(reads)):
            qualities.append(qual)
            
    
    # Print header, seq, spacer, quality line
    # Output files
    output = args.out + '.fastq'
    with open(output, 'w') as out1:
        for i in range(0,len(reads)):
            numseq = i + 1
            print('@' + str(numseq) + '\t' + header, file = out1) # header
            print(reads[i], file = out1) # read sequence
            print('+', file = out1) # spacer line
            print(qualities[i], file = out1) # quality line


# %% PAIRED END READS
if args.paired:
     # lists for generated reads
     reads1 = []
     reads2 = []
     count = 0 # counter for the number of times the sequence is "read", as many as desired coverage
     while count <= coverage:    
         
         # Number of reads
         Nreads = int(len(seq)/(2*readlen))
         # Generate splitting positions
         #positions = np.random.randint(0,len(seq),Nreads)
         positions = [int(random.uniform(0,len(seq))) for i in range(Nreads)]
         positions.append(0)
         positions.append(len(seq)-100)
         for i in positions: # Start reading in each of the positions
         # Read from position i until reaching read length
             read1 = [seq[i:j] for j in range(i, len(seq) + 1) if len(seq[i:j]) == readlen]
             # Move for paired-end read
             i2 = i + 10
             read2 = [seq[i2:j] for j in range(i2,len(seq) + 1) if len(seq[i2:j]) == readlen]
             if read2:
                 read2_rev = rev_comp(read2[0])
                 reads2.append(read2_rev)
             # Add current read to the reads list
             reads1 += read1
            
                
         count += 1 # add counter

     # Quality line
     if args.q:
         # For reads1
         qualities1 = []
         for read in reads1:       
             quality_values = 0.9*readlen 
             
             randQ = random.uniform(0.0004,0.001)
             error_rates1 = np.random.uniform(0.0001,randQ,size=int(quality_values)) 
             
             error_rates2 = np.random.uniform(randQ,0.0015,size=int(readlen - quality_values))
             
             error_ascii1 = ''
             error_ascii2 = ''
        
             for i in error_rates1:
                 Q = qscore(i)
                 error_ascii1 += Q
             for i in error_rates2:
                 Q = qscore(i)
                 error_ascii2 += Q
             qual = error_ascii1 + error_ascii2
             qualities1.append(qual)
    
        # For reads2
         qualities2 = []
         for read in reads2:       
             quality_values = 0.9*readlen 
             
             randQ = random.uniform(0.0004,0.001)
             error_rates1 = np.random.uniform(0.0001,randQ,size=int(quality_values)) 
             
             error_rates2 = np.random.uniform(randQ,0.0015,size=int(readlen - quality_values))
             
             error_ascii1 = ''
             error_ascii2 = ''
             
             # 
             for i in error_rates1:
                 Q = qscore(i)
                 error_ascii1 += Q
             for i in error_rates2:
                 Q = qscore(i)
                 error_ascii2 += Q
             qual = error_ascii1 + error_ascii2
             qualities2.append(qual)
     else:
         qualities1 = []
         qualities2 = []
         qual = readlen*'I'
         for i in range(len(reads1)):
             qualities1.append(qual)
             qualities2.append(qual)
             
# Print FASTQ

# print header, seq, spacer, quality line
# One file for each reads

    # READS 1
    # Output file
     output1 = args.out + '_1.fastq'
     with open(output1, 'w') as out1:
         for i in range(0,len(reads1)):
             numseq = i + 1
             print('@' + str(numseq) + '_1\t' + header, file = out1) # header
             print(reads1[i], file = out1) # read sequence
             print('+', file = out1) # spacer line
             print(qualities1[i], file = out1) # quality line
         
    # READS 2
    # Output file 
     output2 = args.out + '_2.fastq'
     with open(output2, 'w') as out2:
         for i in range(0,len(reads2)):
             numseq = i + 1
             print('@' + str(numseq) + '_2\t' + header, file = out2) # header
             print(reads2[i], file = out2) # read sequence
             print('+', file = out2) # spacer line
             print(qualities2[i], file = out2) # quality line


