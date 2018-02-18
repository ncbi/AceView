#!/usr/bin/env python
###############################################################################################
# poisson.py -- try to graphically fit the distribution of de-uno introns

# meta-x python-shell  gives an interactive python shell inside emacs
###############################################################################################
import getopt
import os
import os.path
import sys
import csv

import matplotlib as mplt
mplt.use('QTAgg')

import string
import time
import numpy as np

import matplotlib.pyplot as plt


try:
	opts, args = getopt.gnu_getopt(sys.argv[1:], '?ho:i:t:',['input=','nreads=','-title=','help'])
except getopt.GetoptError, err:
	print str(err)
	print "try --help"
	sys.exit(2)

if len(args) > 0:
	print "unknown argument "+args[0]
	print "try --help"
	sys.exit(2)

input_file = sys.stdin
title="histo"

for o, a in opts:
	if o  == "-?" or o == "-h" or o == "--help":
		usage()
        elif o == "-o":
                out_prefix = a
	elif o == "-i":
                input_file = a
	elif o == "-t":
                title = a


tt=np.genfromtxt(input_file,delimiter='\t',autostrip=True,skip_header=1)

x=tt[:,0]
y=tt[:,1]
z = np.polyfit(y, x, 1)
p = np.poly1d(z)

print z ;
plt.plot(y,x,'.', color='green' ,label="Delta")
plt.plot(p(x),x,color='red', label="z") 
# plt.plot(tt[:,0],tt[:,1])
# plt.title(title)
plt.show()
exit (0)

## code from Jan

fig = plt.figure()
ax = fig.add_subplot(111)

plt.plot(y,x, color='green' ,ls=".",label="Delta")
plt.plot(range(len(s)),s, color='red', label="Sigma")

#ax.axhline(y=np.mean(x), color='green', alpha=.3, ls='--', label="Mean Delta")
#ax.axhline(y=np.mean(s), color='red', alpha=.3, ls='--', label="Mean Sigma")

#plt.plot(b,a, ls='--',label="fit Sigma c++")
##plt.plot(d,c, label="LS")
#plt.plot(xp,p(xp), label="fit Delta")
plt.plot(xp2,p2(xp2), ls='--',label="fit Sigma python")

plt.plot(xp2,p2deriv1(xp2), ls='--', label="python first derivative")
#plt.plot(np.linspace(0,len(x),len(first_c)),first_c, ls='--',label="c++ first derivative")

plt.plot(xp2,p2deriv2(xp2), ls='--', label="python second derivative")
#plt.plot(np.linspace(0,len(x),len(second_c)),second_c, ls='--',label="c++ second derivative")

ax.axhline(y=0, color='gray', ls='-')

plt.legend()

for x,i in zip(xp2,zip([p2deriv1(x) for x in xp2],[p2deriv2(x) for x in xp2])):
	if i[0]>0 and i[1]>0:
		ax.axvline(x=x, color='gray', ls='-')
		break

#ax.plot(range(len(d)), d)
#mean = np.mean(d)
#mean = [mean]*len(d)


#n, bins, patches = ax.hist(d, 10000, normed=0, facecolor='green', alpha=0.75)
#ax.plot(range(len(d)), mean)

	
plt.show()


print "done."	


