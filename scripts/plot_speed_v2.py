#!/usr/bin/python2.7

'''
Plot the instant speed fluctuations
Plot some selected cell tracks, centered at 0,0
https://www.itl.nist.gov/div898/handbook/eda/section3/eda35c.htm
'''

import sys,math,os,subprocess,random
#from PIL import Image
import matplotlib as mpl
# mpl.use('Agg')
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
# fig, (ax0, ax1) = plt.subplots(nrows=2,sharex=True, sharey=True)
fig, (ax0, ax1) = plt.subplots(nrows=2,sharex=True)
#fig, ax0 = plt.subplots()

if len(sys.argv) <3:
  print "This is the program 'plot_speed_v2.py'"
  print "Usage: ./plot_speed_v2.py <output figure name> <filename1> <filename2>"
  sys.exit(1)
else:
  figname=sys.argv[1]
  # tracknr=int(sys.argv[2])

filecounter=0
for filename in sys.argv[2:]:

  print "reading file ",filename
  timepoints=[]
  speed=[]
  avspeed=[]
  stdevspeed=[]
  xpos=[]
  ypos=[]
  xpos.append([])
  ypos.append([])
  autocorr=[]
  
  count=0
  
  ##read data from file: store raw cell positions
  with open(filename,"r") as fin:

    #read first line first to get the first time point (could probably more clever, but hey)
    for line in fin:
      line=line.split()
      timepoints.append(int(line[0]))
      xpos[-1].append(float(line[3]))
      ypos[-1].append(float(line[4]))
      count+=1
      break

    #read rest of file
    for line in fin:
      line=line.split()
      if int(line[0])!=timepoints[-1]:
        xpos.append([-1.]*len(xpos[0]))
        ypos.append([-1.]*len(xpos[0]))
        count=0
        timepoints.append(int(line[0]))
      if (timepoints[-1]==timepoints[0]):
        xpos[-1].append(float(line[3]))
        ypos[-1].append(float(line[4]))
      else:
        xpos[-1][int(line[1])-1]=float(line[3])
        ypos[-1][int(line[1])-1]=float(line[4])
      count+=1

  ## DONE PRINTING, NOW CALCULATING
  timeint=float(timepoints[-1]-timepoints[-2])
  burn = int(1000/timeint)
  howmany= int(10000/timeint)
  xpos=xpos[burn:burn+howmany] #trimming data to exclude initial conditions and data overload
  ypos=ypos[burn:burn+howmany]
  
  maxint=len(xpos)
  nrcells=len(xpos[0])
  print maxint, nrcells
  count2=0
  totav=0.
  cellstdev=[0.0]*nrcells
  
  for i in range(maxint-1):
      deltax = [ a-b for a,b in zip( xpos[i+1] , xpos[i] ) ]    #Delta x
      deltay = [ a-b for a,b in zip( ypos[i+1] , ypos[i] ) ]    #Delta y
      speed.append([ math.sqrt(x**2. + y**2.)/timeint for x,y in zip(deltax,deltay) ]) #average speed = displacement/time interval
      avspeed.append( np.mean(speed[-1]) )
      # stdevspeed.append(np.std(speed))
 
  # ax0.plot(timepoints[burn:burn+howmany-1],avspeed,label=filename.split('/')[-1])
  # ax0.plot( [timepoints[burn],timepoints[burn+howmany-1]],[avspeed[-1],avspeed[-1]],label=filename.split('/')[-1])
  
  # ax0.hist(avspeed)
  howmany_bins = 50
  bin_edges = np.linspace(0. , 0.20, howmany_bins)
  hist, bin_edges = np.histogram( avspeed, bins=bin_edges, density=True) #average speed hist
  ax0.plot(bin_edges[:-1],hist)
  # ax0.fill_between(timepoints[burn:burn+howmany-1],[ x-y for x,y in zip(avspeed,stdevspeed)],[ x+y for x,y in zip(avspeed,stdevspeed)], alpha=0.5)
  howmany_bins = 26
  bin_edges = np.linspace(0. , 0.26, howmany_bins)
  hist2, bin_edges2 = np.histogram( [s[5] for s in speed] , bins=bin_edges, density=True) #indiv cells hist
  ax1.plot(bin_edges2[:-1],hist2)
  # ax1.plot(timepoints[burn:burn+howmany-1], [s[5] for s in speed],lw=0.5)

# ax0.set_xlabel('time (MCS)')
ax0.set_ylabel('Average cell')
# ax0.set_xlim(0., 5000.)
ax0.legend()
ax1.set_ylabel('individual cell')
ax1.set_xlabel('speed (pix/MCS)')

plt.suptitle('Instantaneous speed')
fig.savefig(figname, bbox_inches='tight')
plt.show()

'''
bla
'''