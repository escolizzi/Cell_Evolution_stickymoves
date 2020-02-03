#!/usr/bin/python2.7

'''
Plot position in time, takes all cells in a time step and calculates average position
Plots the average of the average position, divided by filename
Also depending on file name plots single cellfile name dependence is VERY brittle.
Notice that 
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

def AverageFrom_lg_lists(lg):
  av_lg=[]
  min_lg=[]
  max_lg=[]
  stdev_lg=[]
  
  len_lg = min( [ len(x) for x in lg ] )
  for i in range(len_lg):
    av_lg.append( np.mean( [x[i] for x in lg] ) )
    stdev_lg.append( np.std( [x[i] for x in lg] ) )
    min_lg.append( min([x[i] for x in lg]) )
    max_lg.append( max([x[i] for x in lg]) )
  ltime = [ 10*x for x in range(len_lg) ]
  
  return ltime,av_lg,stdev_lg,min_lg,max_lg


#########################
###                   ###
### ---   BEGIN   --- ###
###                   ###
#########################

colours=[(0.7,0.13,0.),"royalblue", "darkgoldenrod", (0.,0.5,0.2), "salmon", "lightskyblue","orchid"]
filename=""
fig, (ax0, ax1) = plt.subplots(ncols=2,sharey=True)
plt.subplots_adjust(wspace=0.05)
#fig, ax0 = plt.subplots()
# av_xpos=[[],[],[]]
# av_ypos=[[],[],[]]
# lcount=[[],[],[]]

if len(sys.argv) <3:
  print "This is the program 'plot_diplacement_in_time_v2'"
  print "Usage: ./plot_displ.....py <figure name> <peak_row> <peak_col> <filename1> ..."
  sys.exit(1)
else:
  figname=sys.argv[1]
  peakx=int(sys.argv[2])
  peaky=int(sys.argv[3])
  #tracknr=int(sys.argv[2])

# peakx = 500
# peaky = 0
print "Calculating distance from peak at (row,col) = ", peakx,peaky

lsingle=[]
lgmin4=[]
lg0=[]
lg4=[]
lg6=[]
l1c=[]

filecounter=0
for bla,filename in enumerate(sys.argv[4:]):

    print "reading file ",filename
    # timepoints=[]
    # speed=[]
    # avspeed=[]
    # stdevspeed=[]
    # xpos=[]
    # ypos=[]
    # xpos.append([])
    # ypos.append([])
    # autocorr=[]

    # count=0

    av_xpos=[]
    av_ypos=[]
    min_xpos=[]
    min_ypos=[]
    ltime=[]
    lxpos=[]
    lypos=[]
    c1_pos=[]
    ldist=[]
    av_dist=[]
    min_dist=[]
    
    # read time in the first line
    output = subprocess.Popen(['head', '-1', filename], stdout=subprocess.PIPE).communicate()[0]
    time = int(output.split(' ')[0])
    ltime = [time]        
    with open(filename,"r") as fin:
        print 'opening file:', filename
        #read first line first to get the first time point (could probably more clever, but hey)
        for line in fin:
            line=line.split()
      
            time = int(line[0])
            xpos = float(line[3])
            ypos = float(line[4])
            dist = np.hypot( xpos-peakx,ypos-peaky )
            # print xpos,ypos
            #if time of this line different from time of previous line
            # we do aggreagate statistics
            if time != ltime[-1]:
                
                av_xpos.append(  np.mean(  lxpos )  )
                min_xpos.append( min(lxpos) )
                lxpos=[]
                
                av_ypos.append(  np.mean(  lypos )  )
                min_ypos.append( min(lypos) )
                lypos=[] 
                
                av_dist.append( np.mean(ldist) )
                min_dist.append( min(ldist))
                ldist=[]
                
                #print ltime[-1], min_xpos[-1],av_xpos[-1]
                ltime.append(time)
              
            lxpos.append(xpos)
            lypos.append(ypos)
            ldist.append(dist)
          
    #printthis=0
    print 'yo'
    if 'vlines_50c' in filename: 
         if not lg6:
           printthis=1
         else: 
           printthis=0
         lg6.append(av_dist) 
         color=(200/255.,10/255.,10/255.)
         
    elif 'data_vlines_1c' in filename:
         printthis=1
         print "Hello 1c"
         l1c.append(av_dist)
         color=(150/255.,184/255.,179/255.) #light green somewhat
    else: 
         index=4
         printthis = 0 # zeeeeeroooooo
     
    
    
    # if 'g6_1c_f1000' in filename: 
    #     index=0
    #     printthis=1
    # elif 'g6_50c_f1000' in filename:
    #     index=1
    #     printthis=1
    # else:
    #     index=-1
    #     printthis=1
    #index=bla%len(colours)
    #printthis=1
    
    # if 1:
    if printthis:
        #disregard last time point because algotirthm above exsits by then
        
        # ax0.plot( ltime[:-1],av_ypos, c= colours[index], lw=0.5  )
        # ax1.plot( ltime[:-1],min_ypos, c= colours[index], lw=0.25  )
        # 
        
        
        #ax0.plot( ltime[:-1],av_dist, c= colours[index], lw=0.5  ) # this is now printed at the end as average
        ax1.plot( ltime[:-1],min_dist, lw=0.25, color=color )
    
    
    #if not printthis:
    #    ax1.plot( ltime[:-1],min_ypos, c= colours[index], lw=0.25  )
    #    if not c1_pos:
    #        c1_pos = av_ypos
    #        howmany1c=1
    #        ltime_1c = ltime[:]
    #    else:
    #        c1_pos = [ x+y for x,y in zip(c1_pos, av_ypos) ]   #zp works only over the smallest list, which is exactly what we want
    #        howmany1c+=1
    #        if len(ltime)<ltime_1c:
    #            ltime_1c = ltime[:]
            
if lgmin4:
    #PlotMeanAndStd(lgmin4)
    ltime, av_lgmin4, std_lgmin4,minbla,maxbla=AverageFrom_lg_lists(lgmin4)
    color=(0,151/255.,1.) # blue-ish
    ax0.plot( ltime, av_lgmin4 ,color=color, label='-4',zorder=6)
    ax0.fill_between(ltime, [ x-y for x,y in zip(av_lgmin4,std_lgmin4)],[ x+y for x,y in zip(av_lgmin4,std_lgmin4)], color=color, alpha=0.5 ,zorder=1)
    
if lg0:
    ltime, av_lg0,std_lg0,minbla,maxbla=AverageFrom_lg_lists(lg0)
    color=(240/255.,218/255.,104/255.) #royal yellow
    ax0.plot( ltime, av_lg0, color=color , label='0',zorder=7)
    ax0.fill_between(ltime, [ x-y for x,y in zip(av_lg0,std_lg0)],[ x+y for x,y in zip(av_lg0,std_lg0)] , color=color, alpha=0.5,zorder=2)

if lg4:
    ltime, av_lg4,std_lg4,minbla,maxbla=AverageFrom_lg_lists(lg4)
    color=(252/255.,141/255.,89/255.,) #orange
    ax0.plot( ltime, av_lg4, color=color , label='4',zorder=8)
    ax0.fill_between(ltime, [ x-y for x,y in zip(av_lg4,std_lg4)],[ x+y for x,y in zip(av_lg4,std_lg4)] , color=color, alpha=0.5,zorder=3)

if lg6:
    ltime, av_lg6,std_lg6,minbla,maxbla=AverageFrom_lg_lists(lg6)
    color=(200/255.,10/255.,10/255.) #red
    ax0.plot( ltime, av_lg6, color=color , label='6',zorder=9)
    ax0.fill_between(ltime, [ x-y for x,y in zip(av_lg6,std_lg6)],[ x+y for x,y in zip(av_lg6,std_lg6)] , color=color, alpha=0.5,zorder=4)
if l1c:
    ltime, av_l1c,std_l1c,minbla,maxbla=AverageFrom_lg_lists(l1c)
    #color=(90/255.,180/255.,172/255.) #light seagreen
    #color=(199/255.,234/255.,229/255.) #too light seagreen
    color=(150/255.,184/255.,179/255.) #light green somewhat
    ax0.plot( ltime, av_l1c , color=color, label = 'one\ncell',zorder=5)
    ax0.fill_between(ltime, [ x-y for x,y in zip(av_l1c,std_l1c)],[ x+y for x,y in zip(av_l1c,std_l1c)] , color=color, alpha=0.5,zorder=0)
    
    
  # ax0.plot( ltime,av_xpos, c=colour, lw=0.2  )
  # ax1.plot( ltime,min_xpos, c=colour, lw=0.2  )

# print [len(x) for x in av_xpos]    
# for i in range(len(av_xpos)):
#     print i
#     ax0.plot( timepoints[:len(av_xpos[i])],[ x/float(lcount[i]) for x in av_xpos[i] ], c=colours[i], lw=1.  )
#     ax1.plot( timepoints[:len(av_ypos[i])],[ y/float(lcount[i]) for y in av_ypos[i] ], c=colours[i], lw=1.  )


ax0.set_title('Center of mass', fontsize=16)
ax0.legend(title='$\gamma_{c,m}$')
start, end = ax0.get_xlim() # this and next changes axis location
ax0.xaxis.set_ticks(np.arange(0, 150000, 50000))

ax0.set_xlim([0,75000])
ax0.set_ylim([0,None])
ax0.set_xlabel('time (MCS)', fontsize=14)
ax0.set_ylabel('distance from gradient peak', fontsize=14)

ax1.set_title('Cell closest to\ngradient peak', fontsize=16)
ax1.set_xlim([0,75000])
ax1.set_ylim([0,None])
ax1.set_xlabel('time (MCS)', fontsize=14)
fig.savefig(figname, bbox_inches='tight')
plt.show()
