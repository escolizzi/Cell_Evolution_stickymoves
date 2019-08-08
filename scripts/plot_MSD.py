#!/usr/bin/python2.7

'''
Plot the MSD over time for a simulation
Plot some selected cell tracks, centered at 0,0
'''

import sys,math,os,subprocess,random
#from PIL import Image
import matplotlib as mpl
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

colours=["firebrick", "royalblue", "forestgreen", "darkgoldenrod","salmon"]

filename=""

if len(sys.argv)!=3:
  print "This is the program 'plot_MSD.py'"
  print "Usage: ./plot_MSD.py <filename> <nr of tracks to plot>"
  sys.exit(1)
else:
  filename=sys.argv[1]
  tracknr=int(sys.argv[2])


timepoints=[]
MSD=[]
SDEV=[]
xpos=[]
ypos=[]
xn=[]
xpos.append([])
ypos.append([])
#open filename once to see what is the first time step... for now it is typically zero but in the future?
with open(filename,"r") as fin:

  for line in fin:
    line=line.split()
    inittime=int(line[0])
    xpos[-1].append(float(line[3]))
    ypos[-1].append(float(line[4]))
    break

  for line in fin:
    line=line.split()
    if int(line[0])!=inittime:
      break 
    xpos[-1].append(float(line[3]))
    ypos[-1].append(float(line[4]))

timepoints.append(inittime)
count=0
with open(filename,"r") as fin:
  for line in fin:
    line=line.split()
    if int(line[0])!=timepoints[-1]:
      if len(MSD):  #calculate the average and standard deviation
        MSD[-1]=MSD[-1]/float(count)
        sd=0.0
        for el in xn:
          sd+=(el-MSD[-1])*(el-MSD[-1])
        SDEV.append(math.sqrt(sd/float(count)))
      #reset arrays
      MSD.append(0.0)
      xpos.append([])
      ypos.append([])
      xn=[]
      timepoints.append(int(line[0]))
      count=0
    if int(line[0])!=inittime:
      #store all current displacements
      #print xpos
      xpos[-1].append(float(line[3])-xpos[0][count])
      ypos[-1].append(float(line[4])-ypos[0][count])

      xn.append((float(line[3])-xpos[0][count])*(float(line[3])-xpos[0][count]) + (float(line[4])-ypos[0][count])*(float(line[4])-ypos[0][count])) 
      MSD[-1]+=xn[-1] #add to calculation
      
      count+=1

MSD[-1]=MSD[-1]/float(count)

#MSD=[math.sqrt(x) for x in MSD]
sd=0.0
for el in xn:
  sd+=(el-MSD[-1])*(el-MSD[-1])
SDEV.append(math.sqrt(sd/float(count)))

for ind in range(len(xpos[0])):
  xpos[0][ind]=0.0
  ypos[0][ind]=0.0 

xx=zip(*xpos)
yy=zip(*ypos)
toplot=range(0,len(xx))
random.shuffle(toplot)
for el in toplot[:tracknr]:
  print el

print len(MSD), len(SDEV)
fig, (ax0, ax1) = plt.subplots(ncols=2)
#ax0.plot(timepoints[1:],MSD)
ax0.errorbar(timepoints[1:],MSD, yerr=SDEV, fmt='-o')
#ax0.set_yscale('log')
count=0
count2=0
for xrow,yrow in zip(xx,yy):
  if count in toplot[:tracknr]:
    ax1.plot(xrow,yrow)
    count2+=1
  count+=1
plt.show()
