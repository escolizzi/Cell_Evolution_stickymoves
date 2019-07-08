#!/usr/bin/python2.7

'''
minimal plotting tool for pCellEvol pred prey
the file is organised like this:
for all pop, every so many time steps (typically 1000), each individual is dumped, information is:
Time tau birthday key lock tau_contact J_contact tau_contact J_contact tau_contact J_contact ...
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
### --- BLOB PLOT --- ###
###                   ###
#########################
def blob_plot(ax,data,list_pos, colour,datatype,extra=False):
    lextra=[]
    for d,p in zip(data,list_pos):
      if not d: 
        continue
      #this makes the histogram that constitutes the blob plot
      m=min(d) # the largest value, upper edge of the histogram
      M=max(d) # the lowest value, lower...
      
      #print p
      #his=np.bincount(d)
      #print his
      
      if datatype=='int':
        his, bins = np.histogram(d, bins=[ 5*x for x in xrange(2+ M/5 ) ])
        #if len(d)>4:
          #print bins
          #print his
          #x = (bins[:-1] + bins[1:]) / 2.  #centers 
          #print x
          #sys.exit(1)
        
        #x = (bins[:-1] + bins[1:]) / 2.  #centers 
        x = bins[:-1] # seems to work better than centers
        
        #his=np.trim_zeros(his)
        #his=np.append(his,[0])
        #his=np.append([0],his)
        #x=np.linspace(m-1,M+1,M-m+3 )
        
      elif datatype=='float':
        nbins=50
        x = np.linspace(m,M,nbins) # support for histogram of floats
        
        his,bins = np.histogram(d, bins=nbins)
      
      maxwidht=10
      #maxwidht=20000
      #scale_blob=1.
      max_his=max(his)
      #print max_his
      if max_his>0.:
        scale_blob=maxwidht/float(max_his)
        #scale_blob=10
      else: 
        scale_blob=0.
        print "Warning, max_his is 0." 
        
      shift_his_plus_pos =  [ p + h*scale_blob  for h in his]
      shift_his_minus_pos = [ p - h*scale_blob  for h in his]
            
      color_alpha_blobs=colour[:3],colour[-1]
      facecolor,alpha=color_alpha_blobs # assign color and transparency
      #this is the matplotlib function that does the trick
      
      #print "len x:",len(x), "len shift_his_minus_pos", len(shift_his_minus_pos)
      ax.fill_betweenx(x,shift_his_minus_pos, shift_his_plus_pos, linewidth=0., facecolor= facecolor, alpha=alpha)
      #calculates median or mean, if you want
      if extra=='median': lextra.append( np.median(d) )
      elif extra=='mean': lextra.append( np.mean(d) )
    #and plots it
    if extra != False:
      color = 'orangered'
      #ax.plot(list_pos,lextra, color=facecolor,linestyle='-',marker='D',markersize=5, lw=1.5)
      ax.plot(list_pos,lextra, color=colour,linestyle='None',marker='D',markersize=5, lw=1.5)
      
#########################
###                   ###
### ---   BEGIN   --- ###
###                   ###
#########################


filename=""
lfilename=[]

maxTime=-1.
if len(sys.argv)>1:
  pos=1
  while pos < len(sys.argv) :
    if sys.argv[pos]=='-time':
      pos+=1
      maxTime=int(sys.argv[pos])
    elif sys.argv[pos]=='-filename':
      #print "Changing filename to",
      pos+=1
      filename=sys.argv[pos]
      lfilename.append(filename)
      print filename
    else:
      print "This is the program 'plot_data_cellcount.py'"
      print "It plots the birthday of pred and preys alive at a time point"
      print "Here are options, filename is necessary:"
      print "-filename [filename]"
      print "-time [INT]"
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
lt_time=[]
lt_now=[]

#open filename once to see what is the first time step... for now it is typically zero but in the future?
with open(filename,"r") as fin:
  for line in fin:
    line=line.split()
    inittime=int(line[0])
    break
ltime.append(inittime)

prevtime = inittime

for filename in lfilename:
  print("Opening",filename)
  with open(filename,"r") as fin:
      for line in fin:
        line=line.split()
        time=int(line[0])
        tau=int(line[1])
        birthdate=int(line[2])
        #key=line[3]  #don't care about this right now
        #lock=line[4]
        #contacts=line[5:]
        vx= float(line[-2])
        vy= float(line[-1])
        
        if time != prevtime: 
          lt_time.append(lt_now)
          lt_now=[]
          lpop=[0 for _ in lpop]
          ltime.append(time)
          prevtime=time
          
        lpop[tau]+=1
        
        theta = math.atan2(vy,vx)/math.pi*180.
        lt_now.append(theta)
        
print("Done reading files") 
lpop_time=zip(*lpop_time)
ltime=ltime[:-1]  #removing last data point because it is not appended in the loop above
#This shouldn\t be a big deal, given that I get ~30 poins per time step per file
#because it takes a lot of time to make the blob plot when you put together a bunch of these 
##I will just sample the files up to... 1000?
##lt_time_sampled=[]
#for bla in lt_time:
  #print len(bla)
  ##lt_time_sampled.append(random.sample(bla,1000))

##lt_time=lt_time_sampled
#sys.exit(1)

#print ltime
#print lt_time
#print len(ltime)
#print len(lpop_time[0])
#print len(l_avrgdata_time)

#print [ J_av[1,2] for J_av in l_avrgdata_time ]


fig = plt.figure()
ax1 = fig.add_subplot(111)
#ax1.plot( ltime,lpop_time[1], label="#prey" )
#ax1.plot( ltime,lpop_time[2], label="#pred" )
#ax1.legend()

#ax2 = fig.add_subplot(212)
#THIS BELOW WORKS FINE, BUT NOW WE USE GAMMAs, NOTICE THAT I AM USIGN AVERAGE J VALUES

blobplot_every=10
#blob_plot(ax2,[x[1][0] for x in l_avrgdata_time2][::blobplot_every] ,ltime[::blobplot_every], (a[0].get_color(),0.5),'int')
#blob_plot(ax2,[x[1][1] for x in l_avrgdata_time2][::blobplot_every] ,ltime[::blobplot_every], (b[0].get_color(),0.5),'int')
blob_plot(ax1,lt_time[::blobplot_every] ,ltime[::blobplot_every], (0.9,0.1,0.1,0.5),'float')
#blob_plot(ax2,[x[2][0] for x in l_avrgdata_time2][::blobplot_every] ,ltime[::blobplot_every], (d[0].get_color(),0.5),'int')
#blob_plot(ax2,[x[2][1] for x in l_avrgdata_time2][::blobplot_every] ,ltime[::blobplot_every], (e[0].get_color(),0.5),'int') #overlaps perfectly with [1,2], as it should
#blob_plot(ax2,[x[2][2] for x in l_avrgdata_time2][::blobplot_every] ,ltime[::blobplot_every], (f[0].get_color(),0.5),'int')
#print [x[2][2] for x in l_avrgdata_time2]

#blob_plot(ax2,[x[2][2] for x in l_avrgdata_time2][::10] ,ltime[::10], ('yellow',0.8),'int')

#plt.show()
#print ltime
#sys.exit(1)

#PLOT GAMMAS  
# ~G10~ = <J10> - <J11>/2
# ~G20~ = J20 - J22/2
# ~G12~ = J12 - (J11+J22)/2
          
#ax2.plot( ltime, [ J_av[1,0]-J_av[1,1]/2. for J_av in l_avrgdata_time ], label="G10" )
#ax2.plot( ltime, [ J_av[2,0]-J_av[2,2]/2. for J_av in l_avrgdata_time ], label="G20" )
#ax2.plot( ltime, [ J_av[1,2]- (J_av[1,1]+J_av[1,1])/2. for J_av in l_avrgdata_time ], label="G12" )


ax1.legend()
ax1.set_ylim(-180,180)
title=filename.split('/')[-3]+'_'+filename.split('/')[-2]

fig.suptitle(title)

# Option 2
# TkAgg backend
manager = plt.get_current_fig_manager()
manager.resize(*manager.window.maxsize())
#plt.savefig(title+'.pdf')

plt.show()
sys.exit(1)
