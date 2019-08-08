#!/usr/bin/python2.7

'''
Plot the time until the right nr of cells pass the boundary
Makes a scatterplot
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
  print "Usage: ./plot_finaltime.py <picture file> <filename(s)> "
  sys.exit(1)
else:
  figname=sys.argv[1]

filecounter=1

fig, ax=plt.subplots()
for filename in sys.argv[2:]:

  timepoints=[]
  xpos=[]
  count=0
  
  ##read data from file: store raw cell positions
  with open(filename,"r") as fin:

    #read first line first to get the first time point (could probably more clever, but hey)
    for line in fin:
      line=line.split()
      timepoints.append(int(line[0]))
      xpos.append(float(filecounter)+random.random()*0.6-0.3)

  ##start plotting##

  
  #ax0.plot(timepoints,MSD)
  ax.scatter(xpos,timepoints, c=colours[filecounter],alpha=0.5)
  #ax0.set_yscale('log')

  filecounter+=1


fig.savefig(figname, bbox_inches='tight')
