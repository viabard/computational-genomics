"""
week4 bonus, calculates Nx of a genome assembly, where x is 50, or any other number between 10 and 90)
"""
import argparse, sys, re

parser = argparse.ArgumentParser(description="Reads a fasta contigs file and calculates the N50 (or specified N)")
parser.add_argument('-i', '--infile', type=str, help='fasta file name', default="contigs.fa")
parser.add_argument('-n', '--nvalue', type=int, help='n value', default=50)
parser.add_argument('-s', '--size', type=int, help='size of genome', default=100000)
args = parser.parse_args()
infile = args.infile
nvalue = args.nvalue
genomeSize = args.size

try:
    fileThing = open(infile, 'r')
    if(nvalue <= 10 or nvalue >= 90 or genomeSize <= 0):
        sys.exit()
except(FileNotFoundError):
    sys.exit("File not found... (Default file is contigs.fa). Try --help.")
except(SystemExit):
    fileThing.close()
    sys.exit("Use an N value that is between 10 and 90, and a positive genome size...")
except:
    fileThing.close()
    sys.exit("Something went horribly wrong...")

i = fileThing.read()
fileThing.close()
i = re.split(r"(^>.*$)", i, flags=re.M)
contigs, lengths = [], []
for line in i:
    if(line.startswith('\n')):
        contigs.append(line.replace('\n', ''))
while contigs:
    lengths.append(len(contigs[0]))
    contigs = contigs[1:]
lengths = sorted(lengths, reverse=True)   
total, countingTo = 0, genomeSize*(nvalue/100)
iterator = 0
while(total < countingTo):
    total += lengths[iterator]
    current = lengths[iterator]
    iterator+=1
    if(iterator == len(lengths)-1):
        sys.exit("N value does not exist... Contigs do not cover " + str(nvalue) + r"% of the genome.")
print("N" + str(nvalue) + " = " + str(current) + " bp\nFor genome size: " + str(genomeSize) + " bp")