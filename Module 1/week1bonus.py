"""
week1 bonus, reads a fastq file and generates a file that lists the error probabilties at each base positions
"""
import argparse, gzip #imports (argparse for argument parsing, gzip for unzipping)

def getQualityScore(asciiQualities): 
    """
    Converts the illumina quality ascii scores into the error probability
    Takes in a string and returns a list
    """
    qualities = []
    for character in asciiQualities:
        qualities.append(10**(-(ord(character)-33/10)))
    return qualities

def getProbabilities(fileName):
    """
    opens the fastq file and gets the quality score, and converts it to the probabilities of sequence error
    """
    fileObj = open(fileName, 'r') #opening file variable
    seqList = [] #holds the sequences
    qualList = [] #holds the ascii qualities
    qualities = [] #holds the error probabilities 
    iterator = 0
    seqCount = 0 #how many sequences
    qualCount = 0 #how many sets of ascii qualities (should be the same as number of sequences)
    for line in fileObj:
        if iterator%4 == 3: #every 4th line per read (ascii qualities)
            if type(line) == bytes: #if it's gzipped
                line = line.decode() #decode
            qualList.append(line.strip('\n').upper())
            qualCount += 1
        if iterator%4 == 1: #every second line per read (sequence)
            if type(line) == bytes: #if it was gzipped
                line = line.decode() #decode
            seqList.append(line.strip('\n').upper())
            seqCount += 1
        iterator += 1
    fileObj.close() #closing file variable
    for i in range(len(qualList)): #feeding strings of ascii qualities to get probabilities
        qualities.append(getQualityScore(qualList[i]))
    return seqList, qualities #returns seqlist (list of strings), and qualities (list of lists of floats)

#parsing command line information
parser = argparse.ArgumentParser(description="Reads a fastq file and generates a file that lists the error probabilities at each base position")
parser.add_argument('-f', '--file', type=str, help='fastq file name', default='brca.example.illumina.0.1.fastq')
args = parser.parse_args()
name = args.file

#writing per base probabilities to bonusOutput.txt
fileThing = open('bonusOutput.txt', 'w') #opening file variable
probabilities = getProbabilities(name)[1] #getting just the probabilities, not sequences
print("Writing to bonusOutput.txt...")
for i in range(len(probabilities)):
    fileThing.writelines("Read " + str(i) + ": " + str(probabilities[i]) + "\n")
fileThing.close() #closing file variable