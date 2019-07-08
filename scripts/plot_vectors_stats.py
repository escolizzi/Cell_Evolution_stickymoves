#!/usr/bin/python2.7

'''
minimal plotting tool for pCellEvol pred prey
the file is organised like this:
for all pop, every so many time steps (typically 1000), each individual is dumped, information is:
Time tau birthday key lock tau_contact J_contact tau_contact J_contact tau_contact J_contact ...
'''

import sys,math,os,subprocess
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


filename=""

maxTime=-1.
if len(sys.argv)>1:
  pos=1
  while pos < len(sys.argv) :
    if sys.argv[pos]=='-filename':
      #print "Changing filename to",
      pos+=1
      filename=sys.argv[pos]
      print filename
    else:
      print "This is the program 'plot_vectors_stats.py'"
      print "Here are options, filename is necessary:"
      print "-filename [filename]"
      sys.exit(1)
    pos +=1

ldata=[]
ltime=[]
l_avrgdata_time=[]
lcounter=np.zeros((3,3),dtype=int)
lJmatrix=np.zeros((3,3),dtype=int)

l_avrgdata_time2=[]
lcounter2=np.zeros((3,3),dtype=int) #version 2 keeps count of each single value,so that proportions can be made for blob or scatter plot 
                                    # *** actually this might be useless, as blobplot funtion makes internal histogram

lJmatrix2 = [[[] for j in xrange(3)] for i in xrange(3)]

#print lJmatrix2
#sys.exit(1)


lgamma_time=[]
lgamma=3*[0]  # this is G01, G02, G12
lgammapop=3*[0]
lpop=[0,0,0]
lpop_time=[]


#open filename once to see what is the first time step... for now it is typically zero but in the future?
with open(filename,"r") as fin:
  for line in fin:
    line=line.split()
    inittime=int(line[0])
    break
ltime.append(inittime)

ldata=[]

with open(filename,"r") as fin:
      for line in fin:
        line=line.split()
        ldata.append(line)
ldata=zip(*ldata)

lx=ldata[-2]
ly=ldata[-1]
print lx
print ly

ltheta=[ math.atan2(float(y),float(x))/math.pi*180. for x,y in zip(lx,ly) ]

#print ltime

#print len(ltime)
#print len(lpop_time[0])
#print len(l_avrgdata_time)

#print [ J_av[1,2] for J_av in l_avrgdata_time ]

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.hist(ltheta, normed=True, bins=360)


plt.show()
sys.exit(1)
