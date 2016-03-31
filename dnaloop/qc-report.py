import subprocess                                     
import os
import pandas as pd
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

#Grab Sample Names
def get_subdirectories(dir):
    return [name for name in os.listdir(dir)
            if os.path.isdir(os.path.join(dir, name))]

#Read in the summary statistics            
statslist = []
sampleNames = get_subdirectories("samples/")
for s in sampleNames:
	sfile = open("samples/" + s + "/read_stats.txt", 'r')
	stats = []
	header = sfile.readline()
	for l in sfile.readlines():
		stats.append(int(l.split(": ")[1].strip()))	
	statslist.append(stats)
readstats = pd.DataFrame(statslist)
readstats = readstats.transpose()
readstats.columns = sampleNames

#Parse .bedpe files
for s in sampleNames:
	sfilename = s + ".loop_counts.bedpe"
	dist_counts = subprocess.check_output("awk '$1 == $4 {print $5-$2 \"\t\" $8}' " + sfilename, shell=True)
	print dist_counts

# Now make the plots
with PdfPages('qc-report.pdf') as pdf:

    plt.figure(figsize=(8, 11))
    readstats.ix[0].plot(kind='bar')
    plt.title('Total Reads Per Sample')
    
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()

    plt.rc('text', usetex=True)
    plt.figure(figsize=(3, 3))
    x = np.arange(0, 5, 0.1)
    plt.plot(x, np.sin(x), 'b-')
    plt.title('Page Two')

    pdf.savefig()
    plt.close()




