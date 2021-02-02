"""
Assignment 2, simulate fasta data, but make a mistakes with a given error rate
"""
import argparse, random, sys

def getReads(error):
    nucleotides = {'A', 'T', 'C', 'G', 'N'}
    print("Reading stdin, if no input, cancel program and pipe stdin (like cat file.fasta | python3 week2bonus.py -e 0.0001)")
    f = sys.stdin.read()
    reads = f.splitlines(keepends=True)
    newReads = []
    sys.stdout.write(reads[0] + "w/ error rate: " + str(error) + '\n')
    for i in range(len(reads)-1):
        newReads.append('')
        for nucleotide in reads[i+1]:
            if nucleotide != '\n':
                if random.random() <= error:
                    nucleotide = nucleotide.upper()
                    nucleotides.remove(nucleotide)
                    newReads[i] += str(random.sample(nucleotides, 1)[0])
                    nucleotides.add(nucleotide)
                else:
                    newReads[i] += str(nucleotide)
            else:
                newReads[i] += '\n'
        sys.stdout.write(newReads[i])
    return newReads

def stdout(reads):
    for line in reads:
        sys.stdout.write(line)

#parse arguments
try:
    parser = argparse.ArgumentParser(description="A program that simulates fasta short-read data with an error rate in stdout from fasta data given through stdin")
    parser.add_argument('-e', '--errorRate', type = float, help = "the percent error for each base (less than 1)")
    args = parser.parse_args()

    error = args.errorRate
    reads = getReads(error)
except:
    print("Try <python3 week2bonus.py -h>")
