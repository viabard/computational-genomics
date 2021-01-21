#Assignment 1, part 2 

import argparse, gzip, re, tqdm

def checkFileType(fileName):
    """
    Returns a boolean touple:
        First - True is fasta, False is fastq
        Second - True is zipped, False is not zipped
    """
    if fileName.endswith('.gz'): #if it is gzipped
        fileObj = gzip.open(fileName, 'rb')
        zipped = True
    else: #else, just open the file regularly
        fileObj = open(fileName, 'r+')
        zipped = False
    i = 0
    for line in fileObj:
        if type(line) == bytes:
            line = line.decode()
        if i == 2:
            if line.strip('\n').endswith('+'):
                return False, zipped
            return True, zipped
        i+=1

def reverseCompliment(sequence): 
    sequence = sequence.upper()
    revCompliment = ''
    rcs = {}
    rcs['A'] = 'T'
    rcs['T'] = 'A'
    rcs['U'] = 'A'
    rcs['G'] = 'C'
    rcs['C'] = 'G'
    rcs['Y'] = 'R'
    rcs['R'] = 'Y'
    rcs['S'] = 'S'
    rcs['W'] = 'W'
    rcs['K'] = 'M'
    rcs['M'] = 'K'
    rcs['B'] = 'V'
    rcs['D'] = 'H'
    rcs['H'] = 'D'
    rcs['V'] = 'B'
    rcs['N'] = 'N'
    for character in sequence:
        if character in rcs.keys():
            revCompliment += rcs.get(character) 
        else:
            revCompliment += character
    return revCompliment[::-1]

def rcCount(sequence, fileName):
    attributes = checkFileType(fileName)
    choppedSeqList = [[]] #list to hold all of the sequences
    seqList = []
    rcList = [] #list to hold all of the reverse compliments
    counts = [0, 0]
    regex = re.compile(sequence.upper())
    if attributes[1] == True: #it is zipped
        fileObj = gzip.open(fileName, 'rb')
    else: #it is not zipped
        fileObj = open(fileName, 'r')
    if attributes[0] == True: #It is a fasta file
        count = -1
        for line in fileObj:
            if type(line) == bytes: #if it was gzipped
                line = line.decode() #decode
            if line[0] == '>':
                count+=1
            else:
                choppedSeqList[count].append(line.strip('\n').upper())
        for seq in choppedSeqList:
            seqList.append(''.join(seq))
    else: #It is a fastq file
        iterator = 0
        count = 0
        for line in fileObj:
            if iterator%4 == 1:
                if type(line) == bytes: #if it was gzipped
                    line = line.decode() #decode\
                seqList.append(line.strip('\n').upper())
                count += 1
            iterator += 1
    iterator = 0
    for seq in seqList:
        counts[0] += len(re.findall(regex, seq))
        rcList.append(reverseCompliment(seq))
        counts[1] += len(re.findall(regex, rcList[iterator]))
        iterator+=1
    print(seqList, rcList, counts)
        
#parse arguments
parser = argparse.ArgumentParser(description="Generates reverse compliment of a fasta/fastq file and counts for a sequence in each")
parser.add_argument('-s', '--sequence', type=str, help='the sequence that is searched for', default="ATGTTG")
parser.add_argument('-f', '--file', type=str, help='fasta/fastq file name')
args = parser.parse_args()

seq = args.sequence.upper()
name = args.file

#test files
fastq = 'brca.example.illumina.0.1.fastq'
fasta = 'GRCH38p2_chr22.fasta'
dromel = 'Dromel_Adh.fasta'
p53 = 'ensembl_nucleotide_aligned_noDash.fasta'
p53more = 'p53_homologous_unedited_nt.fasta'
fastagz = 'GRCH38p2_chr22.fasta.gz'
testing = 'testing.fasta'
testingq = 'testing.fastq'

#test print statement
#print(reverseComplimentAndCount("", fastagz))
rcCount("ATGTTG", fastq)