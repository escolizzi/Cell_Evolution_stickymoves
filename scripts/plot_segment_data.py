#!/usr/bin/python2.7

'''
Plot the length of straight segments of the cell paths against the angle of the segment with the gradient direction

'''

import sys,math,os,subprocess,random
#from PIL import Image
import matplotlib as mpl
# mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from scipy import stats
from matplotlib.lines import Line2D
import numpy as np

#manager = plt.get_current_fig_manager()
#print manager
#sys.exit(1)

def density_estimation(m1, m2):
    xmin=0
    xmax=np.pi
    ymin=0
    ymax=35
    X, Y = np.mgrid[xmin:xmax:40j, ymin:ymax:35j]                                                     
    positions = np.vstack([X.ravel(), Y.ravel()])                                                       
    values = np.vstack([m1, m2])                                                                        
    kernel = stats.gaussian_kde(values)                                                                 
    Z = np.reshape(kernel(positions).T, X.shape)
    return X, Y, Z


#########################
###                   ###
### ---   BEGIN   --- ###
###                   ###
#########################

colours=["firebrick","royalblue", "darkgoldenrod", "green", "salmon", "lightskyblue","orchid"]
filename=""
# fig, (ax,ax1) = plt.subplots(nrows=2,sharex=True)
# fig.subplots_adjust(hspace=0)

fig, (ax1) = plt.subplots()

#fig, ax0 = plt.subplots()

if len(sys.argv) <5:
  print "This is the program 'plot_segment_data.py'"
  print "Usage: ./plot_segment_data.py <output figure name> <peak_row> <peak_col> <filename> <filename> <filename> ..."
  sys.exit(1)
else:
  figname=sys.argv[1]
  peak_row=int(sys.argv[2])
  peak_col=int(sys.argv[3])

threshold=3.0 	
# print "# WARNING: " 
# print "Warning, distance uses only column, for linear gradient"
# print "# WARNING: "
segsize=[]
segangle=[]
smo_r=[]
smo_c=[]
count=0
for filename in sys.argv[4:]:
  rowpos=[]  #this is an empty list, evaulates to False
  colpos=[]  

  xn=[]
  
  cpos_thistime=[]
  rpos_thistime=[]
  
  colpos.append([]) #this is a list with one element, which is empty - evaluates to True
  rowpos.append([])
  
  count=0

  print filename

  ##################################################
  ##read data from file: store raw cell positions ##
  ##################################################
  
  #read time in the first line
  output = subprocess.Popen(['head', '-1', filename], stdout=subprocess.PIPE).communicate()[0]
  time0 = int(output.split(' ')[0])
  timepoint=time0
  #read file
  with open(filename,"r") as fin:
    for line in fin:
      line=line.split()
      if int(line[0])!=timepoint:
        colpos.append(cpos_thistime)
        rowpos.append(rpos_thistime)
        cpos_thistime=[]
        rpos_thistime=[]
        timepoint=int(line[0])
      # print line[7],line[8]
      rowpos[-1].append(float(line[3]))
      colpos[-1].append(float(line[4]))
  
  #now create straight segments of various lengths
  seg_stor_r=[]
  seg_stor_c=[]	#this is agnostic about the cell a segment belongs to
  seg_pos_r=[]
  seg_pos_c=[]
  print len(rowpos), len(rowpos[0])
  cellrow=[list(i) for i in zip(*rowpos)]
  cellcol=[list(i) for i in zip(*colpos)]
  #go through stored data and look for straight segments
  for i,(cr,cc) in enumerate(zip (cellrow,cellcol)):
    if i>50:
      break
    #print "cell is ",i
    startpoint=0
    if count==0 and i==0:        
      smo_r.append(cr[startpoint])
      smo_c.append(cc[startpoint])
    while startpoint<len(cr)-2:
      error=0 #distance of intermediate points to line segment
      interval=1
      #print "start ", startpoint
      while ( error<threshold and (startpoint+interval)<len(cr)-1 ):
        interval+=1
        distr=cr[interval+startpoint]-cr[startpoint]
        distc=cc[interval+startpoint]-cc[startpoint]
        #calculate the distance of intermediate points from the line
        for (ir,ic) in zip(cr[startpoint+1:startpoint+interval],cc[startpoint+1:startpoint+interval]):
          if np.hypot(distr,distc)>0.:
            nuerror=abs(distr*ic-distc*ir+cr[startpoint]*cc[interval+startpoint]-cc[startpoint]*cr[interval+startpoint])/np.hypot(distr,distc)
          else:
            nuerror=0.
          if nuerror>error:
            error=nuerror

      interval-=1
      #print interval, startpoint, startpoint+interval, len(cr)
      seg_stor_r.append(cr[interval+startpoint]-cr[startpoint])
      seg_stor_c.append(cc[interval+startpoint]-cc[startpoint])
      seg_pos_r.append(cr[startpoint])
      seg_pos_c.append(cc[startpoint])
      if count==0 and i==0:
        smo_r.append(cr[interval+startpoint])
        smo_c.append(cc[interval+startpoint])
      startpoint=startpoint+interval

  #if count==0:
    #ax1.plot(cellcol[0],cellrow[0], color='forestgreen', marker='o', linestyle='solid',linewidth=1, markersize=6)
   # c = np.linspace(0, 10, len(cellcol[0]))
   # ax1.set_aspect('equal')
   # ax1.scatter(cellcol[0],cellrow[0], c=c, cmap=cm.viridis, s=6)
   # ax1.plot(smo_c,smo_r, color='darkgrey', marker='o', linestyle='dashed',linewidth=1, markersize=4)
   # count+=1

  #done with storing the segments, time for plotting
  for (segr, segc, posr,posc) in zip(seg_stor_r, seg_stor_c,seg_pos_r, seg_pos_c):
    segsize.append(np.hypot(segr,segc))
    gradrow=peak_row-posr
    gradcol=peak_col-posc
    segangle.append( math.acos( (gradrow*segr + gradcol*segc) / (segsize[-1]*np.hypot(gradrow,gradcol)) ) )


print len(segsize), len(segangle)
lhist=[[] for i in range(20)]
stor=[0]*20
for (ang,siz) in zip(segangle,segsize):
  binny=int(ang/np.pi*20)
  lhist[binny].append(siz)
  stor[binny]+=1
hhist=[]
for (el,num) in zip(lhist,stor):
  hhist.append(np.median(el))
  #print hhist[-1]	, num
#for bla in np.linspace(0.,np.pi,20):
#  lhist.append( np.median( [x for x,y in zip(segsize,segangle) if (y>=bla and y<bla+np.pi/20. ) ] ) )


ax1.set_xticks([0, np.pi/4., np.pi/2., 3*np.pi/4., np.pi])
ax1.set_xticklabels(['$0$', r'$\frac{\pi}{4}$', r'$\frac{\pi}{2}$', r'$\frac{3\pi}{4}$', r'$\pi$'])
ax1.set_ylim(0,35)
ax1.scatter(segangle, segsize, s=5,c='dimgray')
X, Y, Z = density_estimation(segangle, segsize)
ax1.imshow(np.rot90(Z), cmap=plt.cm.gist_earth_r,                                                    
          extent=[0, np.pi, 0, 35])  #plt.cm.gist_earth_r
ax1.contour(X, Y, Z, colors='black') 
ax1.set_aspect(1.0/ax1.get_data_ratio())            
#ax1.set_aspect('auto')                                                              
#ax1.plot(m1, m2, 'k.', markersize=2)    
#ax1.scatter(np.linspace(0.,np.pi,20),hhist,s=18,c='orange')
#ax1.boxplot(lhist,positions=np.linspace(0.,np.pi,20))

# plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, hspace = 0, wspace = 0)
plt.margins(0,0)
fig.savefig(figname)
plt.show()
  
# print np.mean(rowgrad_cm),np.mean(colgrad_cm)
sys.exit(1)
  
  
  
