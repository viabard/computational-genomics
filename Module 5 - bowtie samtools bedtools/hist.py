"""
Makes histograms for all of the bedtools output files (.bt extension). Save histogram as png.
"""
import matplotlib.pyplot as plt, os

fileNames = [i for i in os.listdir() if i.endswith(".bt")] #getting file names that end with bt, for bedtools

for f in fileNames: #for each file with .bt extension
    histx = [] #hist x values (number of reads per base position)
    histy = [] #hist y values (number of base positions)
    with open(f, 'r') as fileThing: #open the file
        for line in fileThing:
            if line.startswith("genome"): #each line that starts with genome
                histx.append(int(line.split('\t')[1])) #add number of reads per base position
                histy.append(int(line.split('\t')[2])) #add number of base positions for that number of reads
        plt.bar(histx, histy) 
        plt.title(f + " Hist.")
        plt.xlabel("number of reads per base position")
        plt.ylabel("number of base positions")
        plt.xticks(histx)
        plt.subplots_adjust(left = 0.143, right = 0.952) #adjusts the plot for better saving
        plt.savefig(f.strip(".bt") + ".png") #save the histogram as the filename but as png
        plt.close()