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
###     function      ###
###                   ###
#########################

def blob_plot(ax,data,list_pos, extra=False):
    lextra=[]
    for d,p in zip(data,list_pos):
      
      #this makes the histogram that constitutes the blob plot
      m=min(d)
      M=max(d)
      nbins=20
      x = np.linspace(m,M,nbins) # support for histogram
      his,bins = np.histogram(d, bins=nbins)
      
      #scale the histogram to a proper size (you may have to fiddle with this), then place it at position 'pos'
      scale_blob=2000.
      shift_his_plus_pos =  [ p + h*scale_blob/float(len(d))  for h in his]
      shift_his_minus_pos = [ p - h*scale_blob/float(len(d))  for h in his]
      
      facecolor,alpha=color_alpha_blobs # assign color and transparency
      #this is the matplotlib function that does the trick
      ax.fill_betweenx(x,shift_his_minus_pos, shift_his_plus_pos, linewidth=0.0, facecolor= facecolor, alpha=alpha)
      #calculates median or mean, if you want
      if extra=='median': lextra.append( np.median(d) )
      elif extra=='mean': lextra.append( np.mean(d) )
    #and plots it
    if extra != False:
      color = 'orangered'
      ax.plot(list_pos,lextra, color=color_mean_or_median,linestyle='-',marker='D',markersize=5, lw=1.5)
      

#########################
###                   ###
### ---   BEGIN   --- ###
###                   ###
#########################
color_alpha_blobs=('midnightblue',0.7)	
colours=["firebrick","royalblue", "darkgoldenrod", "green", "salmon", "lightskyblue","orchid"]
filename=""
fig, ax = plt.subplots()
#fig, ax0 = plt.subplots()

if len(sys.argv) <3:
  print "This is the program 'plot_migrationangles.py'"
  print "Usage: ./plot_migrationangles.py <figure name> <refvecx> <refvecy> <interval> <filename(s)> "
  sys.exit(1)
else:
  figname=sys.argv[1]
  xx=float(sys.argv[2])
  yy=float(sys.argv[3])
  interval=int(sys.argv[4])
filecounter=0
for filename in sys.argv[5:]:

  print "reading file ",filename
  timepoints=[]
  xpos=[]
  ypos=[]
  angles=[]
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

  ##calculate angles of migration
  maxint=len(xpos)
  nrcells=len(xpos[0])
  print maxint
  for i in range(1,maxint):
    count=0
    count2=0
    angles.append([])
    while (count+i<maxint):
      
      for c in range(nrcells):  #problem when cells die...
        if xpos[count+i][c]>-1:
          veclength= math.sqrt( (xpos[count+i][c]-xpos[count][c])*(xpos[count+i][c]-xpos[count][c]) + (ypos[count+i][c]-ypos[count][c])*(ypos[count+i][c]-ypos[count][c]) )
          if (veclength>0.0001):
            angles[-1].append( 180*math.acos(( (xpos[count+i][c]-xpos[count][c])*xx+(ypos[count+i][c]-ypos[count][c])*yy ) / veclength)/math.pi )
          else:
            print count, i, c
        #else:
         # print "cell ",c,"is no more"
      count+=i
    #print i, angles[-1]

  #
  data=[el for ind,el in enumerate(angles) if ind%interval==0 ]
  positions=[el for ind,el in enumerate(timepoints) if ind%interval==0 ]   
  blob_plot(ax, data, list_pos=positions)	#makes blob plots, if extra is not spcified, makes no median/mean



  filecounter+=1

plt.title('Migration angle distribution \n reference vector: ('+str(xx)+", "+str(yy)+")",fontsize=18)
plt.xlabel('period of measurement (MCS)', fontsize=14)
plt.ylabel('angle to ref. vector', fontsize=14)
#ax.set_xlim(-1., len(angles)) #to display all of them
#ax.set_ylim(-2., 12.) 


fig.savefig(figname, bbox_inches='tight')
#plt.show()
