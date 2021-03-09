"""
    program for computational genomics
    pull out data for apparent snps  - i.e. where the are more than one base called among the reads
    results can be used for making genotype calls

    usage for  python getgenotypes.py
        3 items required on command line:
            1)name of reference fasta file
            2)name of sam file of mapped reads
            3)first base # of reference fasta file with respect to original reference

        returns a table of results, for each variable SNP in the vcf file, each row of output contains:
        (1) the base position
        (2) the reference base
        (3) a list of read bases
        (4) the qual scores,
        (5) a list of all error probabilities,
        (6) list of error probabilities for reads that match the reference,
        (7) a list of error probabilities for reads that do not match the reference


        example run
        python getgenotypes.py hs37ds_22_40000001_40010000.fa  chr22readsample.sam  40000001
"""
import sys
import math
import argparse

def phredprob(Qchar):
    """
        calculate the probability of an error using the character code from the read data
        assumes that character coding is as follows:
            phredscore = (ASCII value of the character) minus 33
            use the ord() function to get ASCII value
    """
    phredscore = ord(Qchar)-33
    return pow(10,-float(phredscore)/10)


def getrefseq(fastafname):
    """
        return a string containing the first sequence in a fasta file
    """
    f = open(fastafname,"r")
    f.readline()
    lines = []
    for l in f:
        if l=='' or l[0] == '\n' or l[0]=='>':
            break
        else:
            lines.append(l.strip())
    return ''.join(lines)

def getreadinfo(samfname):
    """
        reads a sam file
        pull info on the base position, sequence and quality from a set of reads taken from a sam file
        excludes those reads that have a CIGAR score that is not equal to the length of the sequence
        i.e. avoid reads spanning indels
    """
    f = open(samfname)
    firstbases=[]
    seqs = []
    quals = []
    for line in f:
        if line[0] != '@':
            v = line.split()
            seq = v[9]
            seqlen = len(seq)
            noindelCIGARstr = str(seqlen) + "M"
            if v[5] == noindelCIGARstr:  ## check to see if CIGAR string indicates no indels
                firstbases.append(int(v[3]))
                seqs.append(seq)
                quals.append(v[10])
    return firstbases,seqs,quals


def matchreads(refseq,refbase1num,firstbases,seqs,quals):
    """
        make a list, r
            one item for each base in refseq in order
            each item is a list and contains, in order
                the base number
                the reference base
                a list of bases that occurred in reads
                a corresponding list of quality values that occured in reads
                the alt base
    """
    ## by python numbering the first base in refseq is at position 0
    ## need to renumber of firstbases[] values, so the base positions line up
    numbases = len(refseq)
    r = [[i,refseq[i],[],[]] for i in range(numbases)]  # make a structure to hold everything we need for each base
    numreads = len(firstbases)
    for j in range(numreads):
        k = firstbases[j]
        for ci,c in enumerate(seqs[j]):
            renum1 = (k+ci) - refbase1num
            if 0 <= renum1 < numbases:
                r[renum1][2].append(c)
                r[renum1][3].append(quals[j][ci])
    return r

def getMLs(rr, matchreferrors, alterrors, openFile):
    """Each row: base position, ref allele, alt allele,  likelihoods for each genotype (0,1,2), estimated genotype"""
    referenceProbabilities = matchreferrors
    alternateErrors = alterrors
    k = len(alternateErrors) + len(referenceProbabilities)
    probabilities = []
    for genotype in range(3):
        firstProduct = 1
        secondProduct = 1
        for refProb in referenceProbabilities:
            firstProduct *= ((2-genotype)*refProb + genotype*(1-refProb))
        for altProb in alternateErrors:
            secondProduct *= ((2-genotype)*(1-altProb)+genotype*altProb)
        probabilities.append((1/2**k)*firstProduct*secondProduct)
    m = probabilities[0]
    lowest = 0
    for i in range(1,3):
        if probabilities[i] >= m:
            m = probabilities[i]
            lowest = i
    print(" (g0):" + str(probabilities[0]),"(g1):" + str(probabilities[1]),"(g2):" + str(probabilities[2]),str(lowest), file=openFile)


def main(args):


    firstrefbase = args.firstrbase
    samf = args.samf
    ofn = args.outfile

    refseq = getrefseq(args.ref)
    print("reference file loaded:",args.ref)

    firstbases,seqs,quals = getreadinfo(samf)
    print("sam file loaded:",samf," # reads:",len(seqs))

    readinfo = matchreads(refseq,firstrefbase,firstbases,seqs,quals)
    print("recorded bases and qual scores for each position")

    varrr = []
    for rr in readinfo:
        # check to see if there are multiple reads,  and that not all reads match the reference
        isvariable = (len(rr[2]) - rr[2].count('N') > 1) and (rr[2].count(rr[1]) + rr[2].count('N') != len(rr[2]))
        if isvariable:
            varrr.append(rr)
            for v in rr[2]:
                if v != rr[1] and v != 'N':
                    alt = v
                    rr.append(v)
                    break
    f = open(ofn,'w')
    g = open("MLs_" + ofn, 'w')
    print ("positions with variable reads:",file=f )
    print (" for each position print: (1) the base position, (2) the reference base, (3) the alt base,",file=f )
    print ("   (4) the read bases, (5) a list of qual scores, (6) a list of all error probabilities,",file=f )
    print ("   (7) list of error probabilities for reads that match the reference,",file=f )
    print ("    and (8) a list of error probabilities for reads that do not match the reference",file=f )
    for rr in varrr:
        errorprobs = []
        for e in rr[3]:
            errorprobs.append(phredprob(e))
        print(firstrefbase+rr[0],rr[1],rr[4])
        print (firstrefbase+rr[0],rr[1],rr[4],rr[2],rr[3],errorprobs, end='',file=f )

        

        matchreferrors = []
        alterrors = []
        i = 0
        for baseval in rr[2]:
            if baseval == rr[1]:
                matchreferrors.append(errorprobs[i])
            else:
                if baseval != 'N':
                    alterrors.append(errorprobs[i])
            i = i+ 1
        print (matchreferrors,alterrors,file=f   )

        
        print(firstrefbase+rr[0], rr[1], rr[4], end='', file=g)
        getMLs(rr, matchreferrors, alterrors, g)

    f.close()
    g.close()


def createparser():
    """
     """


    parser = argparse.ArgumentParser()
    parser.add_argument("-r",dest="ref",required=True,type=str,help="reference fasta file")
    parser.add_argument("-s",dest="samf",required=True,help=" sam file")
    parser.add_argument("-b",dest="firstrbase",default=1,type=int,help=" first base of reference, with respect to the"
        "full reference that the SAM file numbering is based on (default 1)")
    parser.add_argument("-o",dest="outfile",default="getgenotype_outfile.txt",help="outfile name")
    return parser

if __name__ == '__main__':
    """
        returns a table of results, for each variable SNP in the vcf file, each row of output contains:
        (1) the base position
        (2) the reference base
        (3) a list of read bases
        (4) the qual scores,
        (5) a list of all error probabilities,
        (6) list of error probabilities for reads that match the reference,
        (7) a list of error probabilities for reads that do not match the reference


        example run
        python getgenotypes.py -r hs37ds_22_40000001_40010000.fa  -s chr22readsample.sam  -b 40000001 -o myresults.out
    """
    parser = createparser()
    args, unknown = parser.parse_known_args()
    main(args)



