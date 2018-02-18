#!/bin/env python
import sys, os, time

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

def getCounts(filename):
	forward = [] #[[A,C,G,T,N],[A,C,G,T,N],...,[A,C,G,T,N]]
	reverse = []
	
	current = forward
	
	f = open(filename, 'r')
	for line in f:
		#clear whitespaces
		line = line.strip()
		
		#skip comments
		if line.startswith('#'):
			continue
		
		#did we arrive at the reverse read?
		if line == "":
			current = reverse
			continue
		
		#get the data
		line = line.split()
		current.append(list( map(lambda x: int(x), line[:5])) )

	return [forward,reverse]

if (len(sys.argv)) < 2:
	print ("Usage: %s path_to_input_file.txt [path_to_output.pdf]" % (sys.argv[0]))


else:
	#get parameters
	input_file = sys.argv[1]
	
	output_file = ""
	if len(sys.argv) == 3:
		output_file = sys.argv[2]
	else: #we take the input filename prefix
		output_file = os.path.splitext(sys.argv[1])[0] + ".pdf"

	#get the frequencies
	data = getCounts(input_file)
	
	#do we have a reverse lane?
	number_of_subplots = 2 if len(data[1]) != 0 else 1
	
	#plot generatation
	plot_titles = ["Nucleotide Frequency Distribution - Forward Read", "Nucleotide Frequency Distribution - Reverse Read"]
	fig1 = matplotlib.pyplot.figure(figsize=(16.0, 10.0)) if number_of_subplots == 2 else matplotlib.pyplot.figure(figsize=(16.0, 5.0))
	
	for i,v in enumerate(range(number_of_subplots)):
		v = v+1
		ax1 = plt.subplot(number_of_subplots,1,v)
		
		countsA = list(map(lambda x: x[0], data[i]))
		countsC = list(map(lambda x: x[1], data[i]))
		countsG = list(map(lambda x: x[2], data[i]))
		countsT = list(map(lambda x: x[3], data[i]))
		countsN = list(map(lambda x: x[4], data[i]))

		x_ticks = range(1,len(data[i])+1)
		
		plt.plot(x_ticks,countsA, "-h", color='g'    , label="A", alpha=.7)
		plt.plot(x_ticks,countsC, "-o", color='y'    , label="C", alpha=.7)
		plt.plot(x_ticks,countsG, "-p", color='blue' , label="G", alpha=.7)
		plt.plot(x_ticks,countsT, "-s", color='r'    , label="T", alpha=.7)
		plt.plot(x_ticks,countsN, "-x", color='gray' , label="N", alpha=.7)
		
		plt.plot(x_ticks,countsA, "h", color='g'    , alpha=1)
		plt.plot(x_ticks,countsC, "o", color='y'    , alpha=1)
		plt.plot(x_ticks,countsG, "p", color='blue' , alpha=1)
		plt.plot(x_ticks,countsT, "s", color='r'    , alpha=1)
		plt.plot(x_ticks,countsN, "x", color='gray' , alpha=1)
		
		#subplot properties
		plt.legend()
		plt.xlabel("Nucleotide Position")
		plt.ylabel("Nucleotide Frequency")
		plt.title(plot_titles[i], y=1.10)
		ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.0e'))
		ax1.yaxis.grid(True, color='gray')
		
		x_ticks_labels = range(1,len(data[i])+1)
		ax1.set_xlim([1, len(data[i])])
		plt.xticks(x_ticks_labels, rotation=70)

		#create second xaxis for most frequent nucleotides
		nticks = []
		ncolors = []
		for x in range(len(data[i])):
			temp = zip(['A','C','G','T','N'] , [countsA[x],countsC[x],countsG[x],countsT[x],countsN[x]], ['green','yellow','blue','red','gray'])
			mmax = max(temp, key=lambda x: x[1] )
			nticks.append(mmax[0])
			ncolors.append(mmax[2])

		ax2 = plt.twiny()
		plt2 = ax2.plot(x_ticks,[0 for x in range(len(data[i]))])
		locs,labels = plt.xticks()
		plt.xticks(range(len(data[i])), nticks, size=15)
		ax2.set_xlim([0, len(data[i])-1])
		
		topx = ax2.get_xticklabels()
		for i,a in enumerate(topx):
			a.set_color(ncolors[i])
		
	#prettify overall plot
	plt.tight_layout()
	
	#save figure
	plt.savefig(output_file, format='pdf')
