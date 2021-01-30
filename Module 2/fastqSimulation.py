"""
Assignment 2, fastq data simulation

Works with argument parsing
"""
import argparse, random, subprocess

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

def getReads(inFile, readLength, depth):
    reads = []
    fileThing = open(inFile, 'r')
    fastaRegion = fileThing.readlines()
    nucleotides = ''
    fileThing.close() 
    sp = subprocess.Popen("grep ['>'] -v " + inFile + " | grep [ATCGN] -o | wc -l", shell=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    numberOfBases = int(sp.stdout.read())

    numberOfReads = int(numberOfBases*depth/readLength)
    for i in range(len(fastaRegion)):
        if i != 0:
            nucleotides += fastaRegion[i].strip('\n')
    rcNucleotides = reverseCompliment(nucleotides)
    
    print("Making " + str(numberOfReads) + " reads...")
    for i in range(numberOfReads):
        rand = random.randrange(0, numberOfBases-readLength)
        if (random.randrange(0,1) == 1):
            reads.append(nucleotides[rand:rand+readLength])
        else:
            reads.append(rcNucleotides[rand:rand+readLength])
    return reads

def makeFastq(reads, outFile, inFile, readLength):
    print("Making fastq file: " + outFile)
    i = 0
    out = open(outFile, 'w')
    qual = ''
    for i in range(readLength):
        qual += 'I'
    i = 0
    for read in reads:
        out.write("@Read_" + str(i) + " from " + inFile + "\n")
        out.write(str(read) + '\n')
        out.write("+\n")
        out.write(qual + '\n')
        i+=1

#parse arguments
try:
    parser = argparse.ArgumentParser(description="A program that simulates fastq file short-read data from 100,000 bp regions of Venter's 21st chromosomes.")
    parser.add_argument('-i', '--file', type = str, help = "input file")
    parser.add_argument('-o', '--output', type = str, help = "output file")
    parser.add_argument('-l', '--readlengths', type = int, help = "length of reads in fastq file")
    parser.add_argument('-d', '--depth', type = int, help = "depth of coverage of regions")
    args = parser.parse_args()

    #making 
    inFile = args.file 
    outFile = args.output
    readLength = args.readlengths
    depth = args.depth
        
    reads = getReads(inFile, readLength, depth)
    makeFastq(reads, outFile, inFile, readLength)
except:
    print("python3 fastqSimulation.py -h")