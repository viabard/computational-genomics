"""
    From known barcodes, create fastq files for each barcode with all the reads that have that barcode
"""
barcodes = ['ATGAGATCTT', 'AGCTCATTTC', 'TGAAAATCTT', 'TATCCAGCCA', 'AGGCAGGCAG', 'CTTGTTACTA', 'AAGGCACAAG', 'TGCTCGCTGA', 'GTACCGCCGT', 'CCTCACCAGC']

inFile = "week3_seqs_bc.fq"
allLines, reads = [], []

with open(inFile, 'r') as f:
    """get all of the lines in the file"""
    for line in f:
        allLines.append(line)

while allLines:
    """empty out allLines, and in the process, make a 2D array that holds each read as an array of size 4"""
    reads.append(allLines[0:4])
    allLines = allLines[4:]

#Make empty files to be filled with barcoded sequences
for barcode in barcodes:
    with open(inFile.strip("fq") + barcode + ".fq", 'w') as f:
        pass

for read in reads:
    firstTen = read[1][0:10]
    read[1] = read[1][10:]
    read[3] = read[3][10:]
    with open(inFile.strip("fq") + firstTen + ".fq", 'a') as f:
        f.writelines(read)