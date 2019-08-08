#!/usr/bin/python2.7

'''
Script to merge multiple datafiles for concurrent analysis (general stats like MSD, for which sigma is not used)

'''

import sys,math,os,subprocess,random
#from PIL import Image
import numpy as np

outputfile=open(sys.argv[1],"w")

lastpoint=99999999
timepoint=0

#find the timepoint that all files still have
for filename in sys.argv[2:]:

  with open(filename, "r") as fin:
    thisfile=list(fin)
    firstline=thisfile[0].split()
    lastline=thisfile[-1].split()
    if int(lastline[0])<lastpoint:
      lastpoint=int(lastline[0])
    if int(firstline[0])>timepoint:
      timepoint=int(firstline[0])

print "first and last time point are ",timepoint, " ", lastpoint

#concatenate the files by timepoint
nexttime=timepoint
while timepoint<=lastpoint:

  for filename in sys.argv[2:]:
    #print filename
    with open(filename, "r") as fin:
      for line in fin:
        larr=line.split()
        if int(larr[0])==timepoint:
          outputfile.write(line)
        elif int(larr[0])>timepoint:
          nexttime=int(larr[0])
          print "found nexttime in file ",filename,"it's ",nexttime
          break	
  
  timepoint=nexttime
  print timepoint

outputfile.close()
