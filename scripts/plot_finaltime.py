#!/usr/bin/python2.7

'''
Plot the time until the right nr of cells pass the boundary
Makes a scatterplot
assumes concatenated end times
'''

import sys,math,os,subprocess,random
#from PIL import Image
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from matplotlib.lines import Line2D
import numpy as np

#manager = plt.get_current_fig_manager()
#print manager
#sys.exit(1)
#########################
###                   ###
### ---   BEGIN   --- ###
###                   ###
#########################

colours=["firebrick","royalblue", "darkgoldenrod", "green", "salmon", "lightskyblue","orchid"]
filename=""
fig, (ax0) = plt.subplots(ncols=1)
if len(sys.argv) <3:
  print "This is the program 'plot_finaltime.py'"
  print "Usage: ./plot_finaltime.py <picture file> <sim set <filename(s)>> "
  sys.exit(1)
else:
  figname=sys.argv[1]

filecounter=1

labelset=[]
fig, ax=plt.subplots()
for filename in sys.argv[2:]:

  if len(filename)<4:
    labelset.append(filename)
  else:
    timepoints=[]
    xpos=[]
    count=0
  
    ##read data from file: store raw cell positions
    with open(filename,"r") as fin:

      #read the end points
      for line in fin:
        line=line.split()
        timepoints.append(int(line[0]))
        xpos.append(float(filecounter)+random.random()*0.6-0.3)

    ##start plotting##

  
    #ax0.plot(timepoints,MSD)
    ax.scatter(xpos,timepoints, c=colours[filecounter],alpha=0.5)
    #ax0.set_yscale('log')

    filecounter+=1

ax.set_xlim([0,filecounter+1])
#custom xlabels
my_labels=ax.get_xticks().tolist()
counter=0
newlab=[]
for el in my_labels:
  if el.is_integer() and counter<len(labelset) and el>0.001:
    newlab.append(labelset[counter])
    counter+=1
  else:
    newlab.append('')

ax.set_xticklabels(newlab)
#print my_labels
#print labelset
#print newlab

plt.xlabel('simulation')
plt.ylabel ('migration time')

fig.savefig(figname, bbox_inches='tight')
