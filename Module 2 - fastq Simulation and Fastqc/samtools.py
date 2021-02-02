import os

#First, use samtools faidx to get a 100,000 bp region from each chromosome
print("Using samtools faidx to get region 16741801-16841800 from Venter chromosoms 21...")
os.system("samtools faidx chr21_CM000482.1.fasta CM000482.1:16741801-16841800 -o chr21_CM000482.1.16741801.fasta")
print("Created file chr21_CM000482.1.16741801.fasta and a .fai file.")
os.system("samtools faidx chr21_CM000511.1.fasta CM000511.1:16741801-16841800 -o chr21_CM000511.1.16741801.fasta")
print("Created file chr21_CM000511.1.16741801.fasta and a .fai file.")

