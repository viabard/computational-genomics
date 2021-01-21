"""
week1 bonus, reads a fastq file and generates a file that lists the error probabilties at each base positions
"""
import argparse, gzip, re

def getQualityScore(asciiQualities): 
    qualities = []
    for character in asciiQualities:
        qualities.append(10**(-(ord(character)/10)))
    return qualities

def getProbabilities(fileName):
    """
    opens the fastq file and gets the quality score, and converts it to the probabilities of sequence error
    """
    fileObj = open(fileName, 'r')
    seqList = []
    qualList = []
    qualities = []
    iterator = 0
    seqCount = 0
    qualCount = 0
    for line in fileObj:
        if iterator%4 == 3:
            if type(line) == bytes:
                line = line.decode()
            qualList.append(line.strip('\n').upper())
            qualCount += 1
        if iterator%4 == 1:
            if type(line) == bytes: #if it was gzipped
                line = line.decode() #decode
            seqList.append(line.strip('\n').upper())
            seqCount += 1
        iterator += 1
    for i in range(len(qualList)):
        qualities.append(getQualityScore(qualList[i]))
    return seqList, qualities

parser = argparse.ArgumentParser(description="Reads a fastq file and generates a file that lists the error probabilities at each base position")
parser.add_argument('-f', '--file', type=str, help='fastq file name', default='brca.example.illumina.0.1.fastq')
args = parser.parse_args()
name = args.file

#opening file variable
fileThing = open('bonusOutput.txt', 'w')
probabilities = getProbabilities(name)[1]
print("Writing to bonusOutput.txt...")
for i in range(len(probabilities)):
    fileThing.writelines("Read " + str(i) + ": " + str(probabilities[i]) + "\n")