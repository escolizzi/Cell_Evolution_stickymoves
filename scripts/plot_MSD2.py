#!/usr/bin/python2.7

'''
Plot the MSD over time for one or more simulations.
Plot some selected cell tracks, centered at 0,0
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
fig, (ax0, ax1) = plt.subplots(ncols=2)
if len(sys.argv) <3:
  print "This is the program 'plot_MSD.py'"
  print "Usage: ./plot_MSD.py <nr of tracks to plot> <figure name> <filename(s)> "
  sys.exit(1)
else:
  figname=sys.argv[1]
  tracknr=int(sys.argv[2])

filecounter=0
for filename in sys.argv[3:]:

  timepoints=[]
  MSD=[]
  SDEV=[]
  xpos=[]
  ypos=[]
  xn=[]
  xpos.append([])
  ypos.append([])
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
        xpos.append([])
        ypos.append([])
        count=0
        timepoints.append(int(line[0]))
      xpos[-1].append(float(line[3]))
      ypos[-1].append(float(line[4]))
      count+=1

  ##calculate MSD and standard error
  maxint=len(xpos)
  nrcells=len(xpos[0])
  print maxint
  for i in range(maxint):
    count=0
    count2=0
    MSD.append(0.0)
    sd=0.0
    xn=[]
    while (count+i<maxint):
      for c in range(len(xpos[count+i])):  #problem when cells die...
        xn.append((xpos[count+i][c]-xpos[count][c])*(xpos[count+i][c]-xpos[count][c])+(ypos[count+i][c]-ypos[count][c])*(ypos[count+i][c]-ypos[count][c]))
        MSD[-1]+=xn[-1]
        count2+=1
      count+=1
    MSD[-1]=MSD[-1]/float(count2)
    for el in xn:
      sd+=(el-MSD[-1])*(el-MSD[-1])
    SDEV.append(math.sqrt(sd/float(count2)))

  #center cell track start at 0,0  

  for i in range(1,len(xpos)):
    for j in range(len(xpos[0])): 
      xpos[i][j]-=xpos[0][j]
      ypos[i][j]-=ypos[0][j] 

  for j in range(len(xpos[0])): 
    xpos[0][j]-=xpos[0][j]
    ypos[0][j]-=ypos[0][j] 

  ##start plotting##

  
  #ax0.plot(timepoints,MSD)
  ax0.errorbar(timepoints,MSD, yerr=SDEV, fmt='o', c=colours[filecounter])
  #ax0.set_yscale('log')

  #plot cell tracks
  #invert position matrices for plotting (each row is 1 cell instead of 1 timepoint)
  xx=zip(*xpos)
  yy=zip(*ypos)
  toplot=range(0,len(xx))
  random.shuffle(toplot)
  for el in toplot[:tracknr]:
    print el

  count=0
  count2=0
  for xrow,yrow in zip(xx,yy):
    if count in toplot[:tracknr]:
      ax1.plot(xrow,yrow, c=colours[filecounter])
      count2+=1
    count+=1

  filecounter+=1

ax0.set_xlabel('time (MCS)')
ax0.set_ylabel('MSD (pix^2)')
ax0.set_title('Mean squared displacement')

ax1.set_xlabel('x position')
ax1.set_ylabel('y position')
ax1.set_title('Re-centered cell tracks')

fig.savefig(figname, bbox_inches='tight')
#plt.show()
