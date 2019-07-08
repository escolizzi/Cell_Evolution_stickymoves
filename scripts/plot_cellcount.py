#!/usr/bin/python2.7

'''
minimal plotting tool for pCellEvol pred prey
'''

import sys,math,os,subprocess
#from PIL import Image
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from matplotlib.lines import Line2D
import numpy as np

filename="data_cellcount.txt"

maxTime=9999999999.
if len(sys.argv)>1:
  pos=1
  while pos < len(sys.argv) :
    if sys.argv[pos]=='-maxtime':
      pos+=1
      maxTime=float(sys.argv[pos])
    elif sys.argv[pos]=='-filename':
      print "Changing filename to",
      pos+=1
      filename=sys.argv[pos]
      print filename
    else:
      print "This is the program 'plot_cellcount.py'"
      print "It plots files formatted as time countpreys countpredators"
      print "Here are possible options:"
      print "-maxtime [INT]"
      print "-filename [filename]"
      sys.exit(1)
    pos +=1

ldata=[]
with open(filename,"r") as fin:
      for line in fin:
        line=line.split()
        ldata.append(line)

ldata=zip(*ldata)
ltime=[int(x) for x in ldata[0]]
lprey=[int(x) for x in ldata[1]]
lpred=[int(x) for x in ldata[2]]

#print ltime
#print lprey
#print lpred
#sys.exit(1)
        
######################################
###### ----                ---- ######
###### ---     PLOTTING     --- ######
###### -----              ----- ######
######################################

bla=50
#plt.plot(lprey[:-bla],lpred[:-bla])
#plt.plot(lprey[-bla:],lpred[-bla:])
plt.plot(ltime,lprey,label='Prey')
plt.plot(ltime,lpred,label='Predator')
plt.xlabel('Time (MCS)', fontsize=18)
plt.ylabel('Cellcount', fontsize=18)
plt.legend()
plt.ylim(ymin=0)
plt.show()  










