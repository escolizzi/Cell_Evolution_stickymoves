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
### --- BLOB PLOT --- ###
###                   ###
#########################
      
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
      print "This is the program 'plot_field.py'"
      print "Here are options, filename is necessary:"
      print "-filename [filename]"
      sys.exit(1)
   pos +=1

#open filename once to see what is the first time step... for now it is typically zero but in the future?
#with open(filename,"r") as fin:
  #for line in fin:
    #line=line.split()
    #inittime=int(line[0])
    #break
#ltime.append(inittime)

with open(filename,"r") as fin:
  line=fin.readline()
  line=line.split() #returns split line
  #print line
  time=int(line[1])
  sizex=int(line[3]) -2 #because boundaries
  sizey=int(line[4]) -2
  cpm=[[0 for _ in range(sizey)] for _ in range(sizex)]
  while True: 
    line = fin.readline()
    if line is None: break
    if 'CA' in line: 
      break
  #Now we can read in the CA data
  i=0
  j=0
  while True: 
    line = fin.readline()
    if "IntPlane" in line: break
    line=[int(sigma) for sigma in line.split()]
    for sigma in line:
      if j==sizey:
        j=0
        i+=1
      cpm[i][j]=sigma
      j+=1
  print line
  print "Done reading"

#print cpm  
cpm = zip(*cpm)
cpm = [ bla for bla in cpm[::-1]]
fig, ax = plt.subplots(1, 1)
c = ax.pcolor(cpm)
ax.set_aspect('equal')
plt.show()
sys.exit(1)
