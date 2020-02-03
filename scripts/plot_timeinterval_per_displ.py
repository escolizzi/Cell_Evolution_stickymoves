#!/usr/bin/python2.7

'''
Measure angle of displacement between two time points.

'''

import sys,math,os,subprocess,random
#from PIL import Image
import matplotlib as mpl
# mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from matplotlib.lines import Line2D
import numpy as np

#manager = plt.get_current_fig_manager()
#print manager
#sys.exit(1)

def plot_mean_fillstd(bin_edges,listoflists,color='black',label='label'):
    ll=zip(*listoflists)
    mean_ll=[ np.mean(x) for x in ll ]
    std_ll=[ np.std(x) for x in ll ]
    ax1.plot(bin_edges[:-1], mean_ll,label=label,color=color)
    ax1.fill_between(bin_edges1[:-1], [x-y for x,y in zip(mean_ll,std_ll)],[x+y for x,y in zip(mean_ll,std_ll)],color=color,alpha=0.5)

#########################
###                   ###
### ---   BEGIN   --- ###
###                   ###
#########################

colours=["firebrick","royalblue", "darkgoldenrod", "green", "salmon", "lightskyblue","orchid","black"]
filename=""
# fig, (ax,ax1) = plt.subplots(nrows=2,sharex=True)
# fig.subplots_adjust(hspace=0)

# fig, (ax,ax1) = plt.subplots(nrows=2)

fig, ax = plt.subplots()

if len(sys.argv) <5:
  print "This is the program 'plot_timeinterval_per_displ.py'"
  print "Usage: ./plot_timeinterval_....py <output figure name> <peak_row> <peak_col> <filename> <filename> <filename> ..."
  sys.exit(1)
else:
  figname=sys.argv[1]
  peak_row=int(sys.argv[2])
  peak_col=int(sys.argv[3])

# print "# WARNING: " 
# print "Warning, distance uses only column, for linear gradient"
# print "# WARNING: "
macoldist=0
filecounter=0
l_tot_hist=[]

hist_list_1from50=[]
hist_list_1from50_g6=[]
hist_list_1from50_gmin4=[]
hist_list_1from50=[]
hist_list_1c=[]
hist_list_50c=[]
not_recognised=False
count=0
for filename in sys.argv[4:]:
    rowpos=[]  #this is an empty list, evaulates to False
    colpos=[]  
    rowgrad=[]
    colgrad=[]

    xn=[]

    cpos_thistime=[]
    rpos_thistime=[]
    cg_thistime=[]
    rg_thistime=[]

    colpos.append([]) #this is a list with one element, which is empty - evaluates to True
    rowpos.append([])
    rowgrad.append([])
    colgrad.append([])
    
    timepoints=[]
    
    print "filename: ",filename
    # avangle=[0.]*200*int(binmult)   #1 bin is binsize celldiam 
    # binamount=[0]*200*int(binmult)  #to keep track of how many in each bin
    # angles=[[] for x in range(200*int(binmult))]

    ##################################################
    ##read data from file: store raw cell positions##
    ##################################################

    #read time in the first line
    output = subprocess.Popen(['head', '-1', filename], stdout=subprocess.PIPE).communicate()[0]
    time0 = int(output.split(' ')[0])
    timepoint=time0
    timepoints.append(timepoint)
    
    #read file
    with open(filename,"r") as fin:
        for line in fin:
          line=line.split()
          if int(line[0])!=timepoint:
            colpos.append(cpos_thistime)
            rowpos.append(rpos_thistime)
            colgrad.append(cg_thistime)
            rowgrad.append(rg_thistime)
            cpos_thistime=[]
            rpos_thistime=[]
            cg_thistime=[]
            rg_thistime=[]
            timepoint=int(line[0])
            timepoints.append(timepoint)
          # print line[7],line[8]
          rowpos[-1].append(float(line[3]))
          colpos[-1].append(float(line[4]))
          rowgrad[-1].append(float(line[7]))
          colgrad[-1].append(float(line[8]))
         
    print "Done reading", filename
    
    skip=100
    # howlong=-200 #4000
    howlong=1000 #4000
    #we skip beginning and end
    rowpos=rowpos[skip:skip+howlong]
    colpos=colpos[skip:skip+howlong]
    rowgrad=rowgrad[skip:skip+howlong]
    colgrad=colgrad[skip:skip+howlong]
    
    timepoints=timepoints[:(-2*skip+howlong) ] #this has to be like this because timepoints is also used as x axis
    timeresolution=timepoints[1]-timepoints[0]
    
    len_rpos = len(rowpos)
    nrcells = len(rowpos[0])
    ldist_time = []
    
    rowpos = map(list, zip(*rowpos)) 
    colpos = map(list, zip(*colpos)) 
    for ic in range(nrcells):
        for it1 in range(len_rpos):
            for it2 in range(len_rpos-it1):
                time_interval = it2 # +it1 - it1
                row_displ = rowpos[ic][it2+it1] - rowpos[ic][it1]
                col_displ = colpos[ic][it2+it1] - colpos[ic][it1]
                displ = np.hypot(row_displ,col_displ)
                #displ is on x axis, time interval on y axis
                ldist_time.append( (displ,time_interval*timeresolution) )
    # ldist_time.sort(key=lambda x: x[0])
    # maxdist= ldist_time[-1][0] 
    maxdist=max( [x[0] for x in ldist_time] )
    ldiststep=[]
    lavrg_timeoftravel=[]
    step=1
    for i in range(0,int(maxdist)+1,step):
        ldiststep.append( i+0.5 )
        lavrg_timeoftravel.append( np.mean( [ x[1] for x in ldist_time if (x[0]>=i and x[0]<i+step) ] ) )
            
    # ax.scatter([x[0] for x in ldist_time],[x[1] for x in ldist_time],s=0.01,c=colours[count],alpha=0.5)
    # count +=1
    ax.plot(ldiststep,lavrg_timeoftravel,c=colours[count], label=filename.split('/')[-1])
    count +=1
ax.set_xlabel('distance (pixels)')
ax.set_ylabel('time of travel (MCS)')
ax.legend()
plt.savefig(figname)
plt.show()
print "bye"
sys.exit(1)            

        
        
    
                
                
for _ in range(1):
    #calculate displacement vector for two consecutive time steps
    # then calculate the angle between them
    time_of_dipl=500 # to this you should multiply 10 to get time
    for (rp1,rp2,rp3,cp1,cp2,cp3) in zip(rowpos[:-time_of_dipl],rowpos[time_of_dipl:],rowpos[2*time_of_dipl:],colpos[:-time_of_dipl],colpos[time_of_dipl:],colpos[2*time_of_dipl:]):
        rowdispl_1=[ r2-r1 for (r1,r2) in zip(rp1,rp2) ] 
        coldispl_1=[ c2-c1 for (c1,c2) in zip(cp1,cp2) ] 
        rowdispl_2=[ r3-r2 for (r2,r3) in zip(rp2,rp3) ] 
        coldispl_2=[ c3-c2 for (c2,c3) in zip(cp2,cp3) ] 
        displ1=[np.hypot(x,y) for x,y in zip(rowdispl_1,coldispl_1)]
        displ2=[np.hypot(x,y) for x,y in zip(rowdispl_2,coldispl_2)]
        dotprod=[(rd1*rd2+cd1*cd2)/(d1*d2) for (rd1,rd2,cd1,cd2,d1,d2) in zip(rowdispl_1,rowdispl_2,coldispl_1,coldispl_2,displ1,displ2)  if (d1>0 and d2>0)]
        for bla in dotprod:
            try:
                langle.append(math.acos(bla))
            except:
                pass
        
    #print rowpos
    #print popav_row_direction
    # print follow_your_heart
    howmany_bins = 20
    max = np.pi
    bin_edges = np.linspace(0. , max  , howmany_bins)
    hist1, bin_edges1 = np.histogram( langle, bins=bin_edges, density=True) #indiv cells hist
    # ax.plot(bin_edges1[:-1], hist1,label=filename.split('/')[-1] )#, color=colours[1])
    
    # hist1, bin_edges1 = np.histogram( pop_accuracy, bins=bin_edges, density=True) #pop measure histogram
    
    if '1cellfrom50' in filename: 
        if '1cellfrom50_g6' in filename:
            colorindex=1
            hist_list_1from50_g6.append(hist1)
            
        elif '1cellfrom50_g-4' in filename:
            colorindex=3
            hist_list_1from50_gmin4.append(hist1)
            
        else:
            colorindex=-2
            hist_list_1from50.append(hist1)
    elif 'bc_g6_1c_f500' or 'bc_1c_g6_f500' in filename:
        colorindex=2
        hist_list_1c.append(hist1)
    elif 'bc_g6_50c_f500' in filename:
        colorindex=3
        hist_list_50c.append(hist1)
    else:
        print "Filename not recognised (printed in black)"
        colorindex = -1    
        not_recognised=True
    
    if not_recognised:
        # ax.plot(bin_edges1[:-1], hist2,label=filename.split('/')[-1], color=colours[colorindex])
        ax1.plot(bin_edges1[:-1], hist1,label=filename.split('/')[-1])#,color=colours[colorindex])
        not_recognised=False
    # 
if max == np.pi:
    ax.set_xticks([0, np.pi/4., np.pi/2., 3*np.pi/4., np.pi])
    ax.set_xticklabels(['$0$', r'$\frac{\pi}{4}$', r'$\frac{\pi}{2}$', r'$\frac{3\pi}{4}$', r'$\pi$'])
    ax1.set_xticks([0, np.pi/4., np.pi/2., 3*np.pi/4., np.pi])
    ax1.set_xticklabels(['$0$', r'$\frac{\pi}{4}$', r'$\frac{\pi}{2}$', r'$\frac{3\pi}{4}$', r'$\pi$'])
# ax.set_ylim([0,1.])
ax.set_title('Indiv measure')
ax.legend()
# 
if hist_list_1from50:
    plot_mean_fillstd(bin_edges,hist_list_1from50,color=colours[1],label='1from50')
if hist_list_1from50_g6:
    plot_mean_fillstd(bin_edges,hist_list_1from50_g6,color=colours[1],label='1from50_g6')
if hist_list_1from50_gmin4:
    plot_mean_fillstd(bin_edges,hist_list_1from50_gmin4,color=colours[3],label='1from50_g-4')
if hist_list_1c:
    plot_mean_fillstd(bin_edges,hist_list_1c,color=colours[2],label='1c')
if hist_list_50c:
    plot_mean_fillstd(bin_edges,hist_list_50c,color=colours[0],label='50c')

ax1.set_ylim([0,1])
ax.set_title('Pop measure, time of displ = '+str(time_of_dipl*10)+' MCS')
ax1.legend()
fig.savefig(figname)
plt.show()
sys.exit(1)

'''      
sys.exit(1)

sys.exit(1)
'''  
  
  