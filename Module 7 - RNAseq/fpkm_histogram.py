import math, matplotlib.pyplot as plt
SCALING = 15 #higher number means less bins
FPKM = [math.log10(float(f.split('\t')[-4])+1) for f in open("genes.fpkm_tracking", 'r').readlines()[1:]] #get FPKMs from file, avoiding first header line, taking log10(FPKM+1)
plt.hist(FPKM, int(len(FPKM)/SCALING)), plt.xlabel("log10(FPKM+1)"), plt.ylabel("# of Genes"), plt.title("Gene FPKM Values in " + str(int(len(FPKM)/SCALING)) + " Bins") #making a hist
plt.show()