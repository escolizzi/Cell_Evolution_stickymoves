#!/usr/bin/python2.7

'''
Makes a table with everybody vs eveybody heat map of distances
'''

import sys,math,os,subprocess
#from PIL import Image
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from matplotlib.lines import Line2D
import numpy as np

import KL_adhesion_module as KL

      
#########################
###                   ###
### ---   BEGIN   --- ###
###                   ###
#########################


filename=""
switchfile=""

if len(sys.argv)==1:
  print "Requires at least a file"
  sys.exit(1)
max_hammdist=48
hist_hamm_dist=(len(range(max_hammdist)))*[0]

lJc1_m=[]
lJc1_c2=[]

for filename in sys.argv[1:]:
  ldata=[]
  with open(filename,"r") as fin:
      for line in fin:
        line=line.split()
        
        #print line
        
        key=line[10]  #don't care about this right now
        lock=line[11]
        
        adh=[]
        if line[29]=='0' and len(line)>31: 
          adh = [ int(x) for x in line[30::2]]
          #print adh
        ldata.append( [[int(x) for x in key], [int(x) for x in lock] , adh])
        #lkey.append(key)
        #llock.append(lock)
        
        
      print("Done reading file, now processing")

      ldata_for_hist=[]
      ldata_for_actual=[]
      lhamm_dist=[]
      ldata_2D=[]

      for i,thiskl1 in enumerate(ldata):
        adh=thiskl1[2]
        if len(adh)>0:
          ldata_for_actual.append( adh[0] - np.mean(adh[1:])/2. )
          lJc1_m.append(adh[0])
          lJc1_c2.append(np.mean(adh[1:])/2.)
          #for  bla in adh[1:]:
          #  lJc1_m.append(adh[0])
          #  lJc1_c2.append(bla/2.)
        
        for j,thiskl2 in enumerate(ldata):
          if i==j: continue
          k1 = thiskl1[0]
          l1 = thiskl1[1]
          
          k2 = thiskl2[0]
          l2 = thiskl2[1]
          
          Jc1_m  = KL.JWithMedium(k1,KL.lookuptable_Jmedium)
          Jc1_c2 = KL.JWithOtherTau(( k1,l1 ),( k2,l2 ))
          ldata_for_hist.append( Jc1_m - Jc1_c2/2. )
          
          
          
          if j>i:
            dist=sum([1 if x!=y else 0 for x,y in zip(k1,k2)]) + sum([1 if x!=y else 0 for x,y in zip(l1,l2)])
            #dist=sum([1 if x!=y else 0 for x,y in zip(k1,k2)])
            
            #dist=sum([0 if x!=y else 1 for x,y in zip(k1,l2)])
            #if dist>20: #and dist < 23: 
              #print k1,l1
              #print k2,l2
              #print Jc1_m , Jc1_c2 , Jc1_m - Jc1_c2/2.
              #sys.exit(1)
            lhamm_dist.append(dist)
          
      #print ldata_for_hist
      
      bin_edges=[-6.0,-4.5, -3.,  -1.5,  0.  , 1.5 , 3. ,  4.5 , 6. ,  7.5 , 9. , 10.5, 12., 13.5]
      #hist, bin_edges = np.histogram(ldata_for_hist, bins = bin_edges,density=True)
      #print bin_edges
      #plt.plot(bin_edges[:-1],hist,color='black')

      #hist, bin_edges = np.histogram(ldata_for_actual, bins=bin_edges, density=True)
      #plt.plot(bin_edges[:-1],hist,color='red')
      
      bin_edges=range(max_hammdist)
      hist, bin_edges = np.histogram(lhamm_dist, bins=bin_edges, density=True)
      #plt.plot(bin_edges[:-1],hist,color='blue')
      hist_hamm_dist = [ x+y for x,y in zip(hist,hist_hamm_dist)]

#plt.plot(bin_edges[:-1],[x/float(len(sys.argv[1:])) for x in hist_hamm_dist],color='green',lw=3)


x_bins = np.linspace(8,20,14)
y_bins = np.linspace(0,30,30)

plt.hist2d(lJc1_m,lJc1_c2,bins=[x_bins,y_bins])

#plt.hist(ldata_for_hist, bins=30)
#plt.hist(ldata_for_actual,bins=30)
plt.show()

sys.exit(1)


