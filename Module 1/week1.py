#Assignment 1, part 2 

import argparse, gzip



def checkFileType(fileName):
    """
    Returns True if the file is a fasta file, and False if the file is fastq
    """
    if fileName.endswith('.gz'): #if it is gzipped
        fileObj = gzip.open(fileName, 'rb')
    else: #else, just open the file regularly
        fileObj = open(fileName, 'r+')
    i = 0
    for line in fileObj:
        if type(line) == bytes:
            line = line.decode()
        if i == 2:
            print(line.strip('\n'))
            if line.strip('\n').endswith('+'):
                return False
            return True
        i+=1

#parse arguments
parser = argparse.ArgumentParser(description="Generates reverse compliment of a fasta/fastq file and counts for a sequence in each")
parser.add_argument('-s', '--sequence', type=str, help='the sequence that is searched for', default="ATGTTG")
parser.add_argument('-f', '--file', type=str, help='fasta/fastq file name')
args = parser.parse_args()

seq = args.sequence
name = args.file

#test files
fastq = 'brca.example.illumina.0.1.fastq'
fasta = 'GRCH38p2_chr22.fasta'
fastagz = 'GRCH38p2_chr22.fasta.gz'

#test print statement
print(checkFileType(fastagz))