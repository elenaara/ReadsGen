# -*- coding: utf-8 -*-
"""

-- Assembler --

Created on Thu Mar 10 12:55:05 2022

@author: Elena Aramendía

This script takes reads in a FASTQ format and alignes them to a reference 
sequence, generating an assembled sequence using the found positions of the 
reads. Input reads should all have the same length.

The script splits the reads into k-mers and only searches for exact matches
to these k-mers, so it will not map reads that do not contain exact maches
of k length. This could be improved generating a list of similar k-mers for
each search.



Procedure:
    1. Input files are read. If the input is paired-end reads, the sequences
    in the second file are converted to their reverse complement and stored
    together with the sequences in the first file.
    
    2. Search for the position of each read in the reference sequence:
        - Read is split into kmers of k length (can be chosen by the user, 
          default is 11). 
        - Search for each k-mer in the reference sequence. 
        - If a match is found, the whole read is compared to the reference 
        for that position, and the percentage identity is calculated. 
        - Reads with a percentage of 90% or higher are mapped to that position. 
        If a read has already been mapped (has a percentage of 90% or higher)
        no more k-mers for this read are searched, go to the next one.
        - When we have a read that has matches with high identity, check that
        we do not have another read mapped to that position. If there is
        already a read, compare idenity and keep the best one.
        
    3. Once we have mapped the reads, generate the sequence taking the 
    nucleotides we have mapped to each position. If a position has no 
    nts mapped to it, and 'N' will be written.
    
    4. Print the whole sequence (including the Ns) and the contings (sequences
    wihout the Ns) to output files. If specified with a flag, percentage
    identity between the generated sequence and the reference can be calculated
    and printed.
        

Functions:
    fasta1line:
        Takes fasta file and returns a dictionary, used to get the whole 
        sequence in one line in case it is in multiple lines.
    
    seq_fastq:
        Takes a fastq file and stores only header and sequence lines in a dictionary.
        
    rev_comp:
        Takes a DNA sequence and generates the reverse complement sequence.
    
    k_mers:
        Split sequence into k-mers.
    
    identity:
        Compares two sequences between defined positions and calculates
        percentage identity.
    
Modules:
    argparse
    re


Input:
    - ref sequence (-ref)
    - reads (-i for single-end, -p1 and -p2 for paired-end)
    - k (--k, k-mer length)
    - identity (--identity, flag to print percentage identity between 
                generated sequence and reference)
    - o basename for outputfiles

Output:
    - assembly file (whole sequence with missing nts as 'N')
    - contigs file
Example usage:
    python align_reads.py -ref rCRS.fa -p1 reads1.fastq -p2 reads2.fastq
"""

# %% Modules
import argparse
import re

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

def seq_fastq(fastq): 
    """
    Function to get headers and sequences from a fastq file.
    
    - Input: fastq
        A fastq file.
        In this case headers start with @, and are followed but
        alphanumeric characters and a tab.
        
    - Return: fastq_dict
        Dictionary with header as key and sequence as value.
    
    """
    fastq_dict = {}
    for line in fastq:
        if re.match('@[\w]*\t[\w]*',line):
            header = line.rstrip()
            seq = fastq.readline().rstrip()
            fastq_dict[header] = seq
    return(fastq_dict)

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


def k_mers(seq, k):
    """
    Function to split a sequence into k-mers of k length
    
    Input: 
        - seq (str)
        Nucleotide sequence.
        - k (int)
        Length of the k-mers. Should be an odd number.
    
    - Return: kmer-list
        List with all possible kmers in the provided sequence.
    
    """
    kmer_list = []
    seq = seq.upper()
    for i in range(len(seq) - k):
        for j in range(i, len(seq)):
            if len(seq[i:j]) == k:
                kmer = seq[i:j]
                kmer_list.append(kmer)
    return(kmer_list)

def identity(seq,ref_seq,start,end):
    """
    Compares two sequences between the 'start' and 'end' positions
    and calculates idenitty percetage (percentage of residues that are the same).
    
    Parameters
    ----------
    seq : str
        Shorter sequence
    ref_seq : str
        Longer sequence
    index : int
        Position to start comparing the sequences
    k: int
        Length of the k-mers
    
    Returns
    -------
    Identity percentage

    """
    length = len(seq)
    ident = 0
    j = 0 # counter for shorter sequence
   
    if end < len(ref_seq):
        for i in range(start,end):
            if ref_seq[i] == seq[j]:
                ident += 1
            
            j += 1 # add to counter
        iden_per = ident*100/length
        return iden_per
    else:
        for i in range(start,len(ref_seq)):
            if ref_seq[i] == seq[j]:
                ident += 1
            
            j += 1 # add to counter
        iden_per = ident*100/length
        return iden_per
    

    
# %% Argument parser

parser = argparse.ArgumentParser(description='This program takes a reference sequence and a fastq file (or 2 in case they are paired-end reads) and aligns the reads to the reference if there is enoguh similarity',
                                 epilog='Example usage: python align_reads.py -ref rCRS.fa -p1 reads1.fastq -p2 reads2.fastq')


parser.add_argument(                                                  # Reference fasta file
    '-ref',
    type = argparse.FileType('r'),
    dest = 'ref_seq',
    required = True,
    help = 'Input reference sequence in fasta format'
    )

parser.add_argument(                                                  # Input single end 
    '-i',
    type = argparse.FileType('r'),
    dest = 'reads',
    required = False,
    help = 'Input reads in FASTQ format'
    )

parser.add_argument(                                                  # Input 1 paired-end
    '-p1',
    type = argparse.FileType('r'),
    dest = 'reads1',
    required = False,
    help = 'FASTQ file 1 for paired end reads'
    )

parser.add_argument(                                                  # Input 2 paired end
    '-p2',
    type = argparse.FileType('r'),
    dest = 'reads2',
    required = False,
    help = 'FASTQ file 2 for paired end reads'
    )


parser.add_argument(                                                  # k-mer length
    '--k',
    type = int,
    dest = 'k',
    required = False,
    default=11,
    help = 'K-mer length, must be an integer'
    )

parser.add_argument(                                                  # calculate identities
    '--identity',
    action='store_true',
    dest = 'identity',
    help = 'Print percentage idenity between aligned sequence and reference'
    )

parser.add_argument(                                                  # output 
    '-o',
    type = str,
    dest = 'out',
    required = False,
    default = 'assembly',
    help = 'Basename for output files'
    )


args = parser.parse_args() 
#%% Import sequences

# Reference sequence
#ref = open('C:/Users/earam/Desktop/cole/BINP29/Project/rCRS.fa', 'r')
ref = args.ref_seq
ref_dict = fasta1line(ref)
ref_seq = list(ref_dict.values())[0]
ref_header = list(ref_dict.keys())[0][1:]

# Read sequence
#reads = open('C:/Users/earam/Desktop/cole/BINP29/Project/mit_reads100_10.fastq', 'r')
# Single
if args.reads:
    fastq_dict = seq_fastq(args.reads)
    readlen = len(list(fastq_dict.values())[0])

# Paired reads
elif args.reads1 and args.reads2:
    # merge files into one?
    reads1_dict = seq_fastq(args.reads1)
    reads2_dict = seq_fastq(args.reads2)
    reads2_rev= {} # Paired reads are reversed, get the forward strand
    for key in reads2_dict:
        seq = reads2_dict[key]
        rev_seq = rev_comp(seq)
        reads2_rev[key] = rev_seq
    fastq_dict = {**reads1_dict, **reads2_rev}
    # if onlye one exist raise exception, need both reads
else:
    raise(Exception('Wrong input. Please provide one file with the flag -i or 2 files using the flags -p1 and -p2.'))

# Calculate length of the reads
readlen = len(list(fastq_dict.values())[0])


# %% Align reads 

# Find position of the kmers in the sequence
read_index = {} # dictionary for reads and positions
#k = args.k # k-mer length
k = 11
map_reads = 0 # counter for mapped reads, to calculate how many of the reads were mapped

for key in fastq_dict: # for each read
    seq = fastq_dict[key] # get sequence
    kmers = k_mers(seq,k)
    compared = False
    for kmer in kmers: 
        if compared == False:         
                if kmer in ref_seq: # if a kmer matches the ref sequence
        
                # Find position of all matches
                    ind = [index.start() for index in re.finditer(seq[:k],ref_seq)] 
                    for i in ind: # For each position where there is a match
                    # Compare the whole read with ref and get and identity percentage
                        pos_read = seq.index(kmer)
                        start = i - pos_read
                        end = i + (len(seq) - pos_read)
                    # Compare the whole read
                        iden_per = identity(seq,ref_seq,start,end)
                        
                        # Establish an identity percentage cutoff and compare
                        # If the identity is bigger or equal than our cutoff
                        if iden_per >= 90:
                            # Already compared this sequence
                            compared = True 
                            map_reads += 1
                            # Check if we already have a sequence in this position
                            if start not in read_index: # If not, add sequence and position to a dictionary
                                read_index[start] = seq
                            if start in read_index: # If there is already a sequence for this position
                            # Compare identity percentages and keep sequence with higher identity
                                if iden_per > identity(read_index[start],ref_seq,start,end):
                                    read_index[start] = seq 
                    

#%% Assemble sequence 
sequence = ''
#last_pos = sorted(list(read_index.keys()))[-1]
#assembly_length = last_pos + len(read_index[last_pos])
for i in range(len(ref_seq)): # For each position fo the reference sequence
    if i in read_index: # Search if there is a read mapped to this positions
    # if so, add first nt of this read to the assembled sequence
        position = 0
        index = i
        sequence += read_index[i][position]
    else: # No reads mapped to this position
        if not sequence: # If it is the first position
            sequence += '-'
        elif sequence.count('-') == len(sequence): # First positions
            sequence += '-'
        elif position >= readlen - 1: # If we have aligned all nts from last read
            sequence += '-'
        else: # If we have a previous reads that maps here
            position += 1
            sequence += read_index[index][position]

# Output files

# First file it's whole sequence with the unknown nts (N)
out1=open(args.out + '.fa', 'w')
N = sequence.count('-')
# Show mising positions as 'N'
trans = sequence.maketrans('-','N') 
assembly = sequence.translate(trans)
print('>assembly\t'+ 'L:' + str(len(sequence)) + '\tN:' + str(N) + '\tref:' + ref_header,file=out1)
print(assembly,file=out1)

# Second file, contigs
out2 = open(args.out+'_contigs.fa', 'w')
numseq = 1
for fragment in re.split('-+',sequence):
    if fragment:
        print('>contig_'+ str(numseq) + '\tL:' + str(len(fragment)),file=out2)
        print(fragment,file=out2)
        numseq += 1

# Calculate percentage identity between generated sequence and ref
if args.identity:
    per_i = identity(assembly,ref_seq,0,len(ref_seq))
    per_map = map_reads*100/len(fastq_dict)
    print('Percentage identity: {0:.2f}%'.format(per_i))
    print('Mapped reads: {0:.2f}%'.format(per_map))

