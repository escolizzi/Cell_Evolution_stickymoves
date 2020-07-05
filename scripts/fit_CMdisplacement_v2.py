#!/usr/bin/python2.7

'''
Plot the displacement of the total cell mass over a fixed period for one or more simulations.
simulations can be grouped by parameter set by numbering a subset of the files:
./plot_CMdisplacement.py 5000 1 file1 file2 file2 2 file3 file4 3 file5 file6 ....
'''

import sys,math,os,subprocess,random
#from PIL import Image
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from numpy.polynomial.polynomial import polyfit
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
fig, ax = plt.subplots()
if len(sys.argv) <3:
  print "This is the program 'plot_CMdisplacement.py'"
  print "Usage: ./plot_CMdisplacement.py <figure name> <starttime> <endtime> <(short) simset <filename(s)>> "
  sys.exit(1)
else:
  figname=sys.argv[1]
  starttime=int(sys.argv[2])
  endtime=int(sys.argv[3])
  
xpos=[]
ypos=[]
time=[]
dist=[]

l_m =[]
l_popsize=[]

filecount=0

for filename in sys.argv[4:]:

  if len(filename)<4:
    print "new set. Not doing anything with that right now"
  else:
    filecount+=1
    xpos=[]
    ypos=[]
    time=[]
    dist=[]
    
    count=0
    popsize=0

    ##read data from file: store raw cell positions
    with open(filename,"r") as fin:
      print filename
      #read first line first to get the first time point (could probably more clever, but hey)
      for line in fin:
        line=line.split()
        timepoint=int(line[0])
        if (timepoint==starttime):
          xpos.append(0.)
          ypos.append(0.)
          xpos[-1]+=float(line[3])
          ypos[-1]+=float(line[4])
          time.append(timepoint)
          count+=1
          break

      #read rest of file, to find positions at beginning and end of period
      for line in fin:
        line=line.split()
        if int(line[0])!=timepoint:
          xpos[-1]=xpos[-1]/float(count)
          ypos[-1]=ypos[-1]/float(count)
          dist.append(math.sqrt(xpos[-1]**2+ypos[-1]**2))
          xpos.append(0.)
          ypos.append(0.)
          timepoint=int(line[0])
          time.append(timepoint)
          popsize=count
          count=0
          if(timepoint>=endtime): break
        
        xpos[-1]+=float(line[3]) #append x and y pos
        ypos[-1]+=float(line[4])
        count+=1                 #counts number of cells

    xpos[-1]=xpos[-1]/float(popsize)
    ypos[-1]=ypos[-1]/float(popsize)
    dist.append(math.sqrt(xpos[-1]**2+ypos[-1]**2))    
    
    print len(time), len(dist)
    print "Sampling every 10 points, because it's too many"
    time=time[:-1:10] # -1 is there because last point is empty
    dist=dist[:-1:10]
    
    t=np.array(time) 
    d=np.array(dist)
    b,m=polyfit(t, d, 1)
    
    print filecount, b, m
    l_m.append(-m)
    
    popsize=int(filename.split('/')[-1].split('_')[3][:-1])
    l_popsize.append(popsize)
    
    t=np.array(time)
    #ax.plot(time, dist, color=colours[(filecount-1)%len(colours)], marker='.', markersize=8) 
    #ax.plot(t, b+m*t, color=colours[(filecount-1)%len(colours)],linestyle='solid', linewidth=2)
 


ax.scatter(l_popsize,l_m)

l_popsize_sorted = sorted(list(set(l_popsize)))
l_m_sorted = len(l_popsize_sorted)*[0.]
l_howmany_m = len(l_popsize_sorted)*[0]

for j,p1 in enumerate(l_popsize_sorted):
  for i,(p2,m) in enumerate(zip(l_popsize,l_m)):
    if p1==p2:
      l_m_sorted[j]+=m
      l_howmany_m[j] +=1
  
l_m_sorted = [ x/float(y) for x,y in zip(l_m_sorted,l_howmany_m) ]  

ax.plot(l_popsize_sorted,l_m_sorted)

plt.xlabel('pop size')
plt.ylabel('speed')
fig.savefig(figname, bbox_inches='tight')
#plt.show()



