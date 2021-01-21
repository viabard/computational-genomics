#Assignment 1, part 2 

import argparse, gzip, re, tqdm #imports

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
    """
    Creating the reverse compliment of a sequence(string)...
    Returns the reverse compliment. If there are any dashes or non-nucleotide (non IUPAC), it is just re-added.
    """
    sequence = sequence.upper()
    revCompliment = ''
    temporary = []
    rcs = {} #dictionary of IUPAC nucleotides
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
            temporary.append(rcs.get(character)) 
        else:
            temporary.append(character)
    revCompliment = ''.join(temporary) #join is faster than adding/math00    
    return revCompliment[::-1]

def rcCount(sequence, fileName):
    """
    Finds information from the file and then finds the sequence matches of that sequence, and the reverse compliment.
    """
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
    return [seqList, rcList, counts]
        
#parse arguments
parser = argparse.ArgumentParser(description="Generates reverse compliment of a fasta/fastq file and counts for a sequence in each")
parser.add_argument('-s', '--sequence', type=str, help='the sequence that is searched for', default="ATGTTG")
parser.add_argument('-f', '--file', type=str, help='fasta/fastq file name', default='brca.example.illumina.0.1.fastq')
args = parser.parse_args()
seq = args.sequence.upper()
name = args.file


#oppening file variable
fileThing = open("output.txt", 'w')

print("Checking", name, "and it's reverse compliment for", seq, "...")
fileThing.write("Checking "+ name+ " and it's reverse compliment for " + seq + "...\n")

#extracting the information from the function return
tempList = rcCount(seq, name)
seqList = tempList[0]
rcList = tempList[1]
counts = tempList[2]

#printing and writing to file the counts
print(name, "counts:", counts[0]) 
print(name, "reverse compliment counts:", counts[1])
fileThing.write(''+ name + " counts: " + str(counts[0]) + '\n')
fileThing.write(''+ name + " reverse compliment counts: " + str(counts[1]))

print("Output in output.txt") #output.txt
fileThing.close()