#!/usr/bin/python2.7

'''
Measure how much gradient measure now is different from displacement later.

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

fig, (ax,ax1) = plt.subplots(nrows=2)

#fig, ax0 = plt.subplots()

if len(sys.argv) <5:
  print "This is the program 'plot_error_....py'"
  print "Usage: ./plot_error_....py <output figure name> <peak_row> <peak_col> <filename> <filename> <filename> ..."
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
hist_list_1c=[]
hist_list_50c=[]
not_recognised=False

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

    count=0

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
          # print line[7],line[8]
          rowpos[-1].append(float(line[3]))
          colpos[-1].append(float(line[4]))
          rowgrad[-1].append(float(line[7]))
          colgrad[-1].append(float(line[8]))
         
    print "Done reading", filename
    
    skip=100
    howlong=4000
    #we skip beginning and end
    rowpos=rowpos[skip:skip+howlong]
    colpos=colpos[skip:skip+howlong]
    rowgrad=rowgrad[skip:skip+howlong]
    colgrad=colgrad[skip:skip+howlong]
    
    #Results lists
    popav_row_direction =[]
    popav_col_direction =[]
    
    #center of mass of the blob
    peakdir_row =[]
    peakdir_col =[]
    
    row_displ=[]
    col_displ=[]
    
    #calculate displacement vector
    time_of_dipl=100 # to this you should multiply 10 to get time
    for (rp1,rp2,cp1,cp2) in zip(rowpos[:-time_of_dipl],rowpos[time_of_dipl:],colpos[:-time_of_dipl],colpos[time_of_dipl:]):
        rowdispl_now=[ r2-r1 for (r1,r2) in zip(rp1,rp2) ] 
        coldispl_now=[ c2-c1 for (c1,c2) in zip(cp1,cp2) ] 
        row_displ.append( rowdispl_now )
        col_displ.append( coldispl_now )
    #calculate average local gradient measure over time_of_dipl
    t_rg=zip(*rowgrad)
    t_cg=zip(*colgrad)
    time_rowgrad=[]
    time_colgrad=[]
    for rg_icell in t_rg:
        time_rowgrad.append([ np.mean( rg_icell[ i:i+ time_of_dipl]) for i in range( len(rg_icell) - time_of_dipl ) ])
    for cg_icell in t_cg:
        time_colgrad.append([ np.mean( cg_icell[ i:i+ time_of_dipl]) for i in range( len(cg_icell) - time_of_dipl ) ])
    time_rowgrad=map(list,zip(*time_rowgrad))
    time_colgrad=map(list,zip(*time_colgrad))
    print len(row_displ)
    print len(col_displ)
    print len(time_rowgrad)
    print len(time_colgrad)
    # sys.exit(1)
    #quantify how much cells have gone where they felt like going according to measured gradient: dot product
    follow_your_heart=[]
    for rd1,cd1,rg1,cg1 in zip(row_displ,col_displ,time_rowgrad,time_colgrad):
        # print rd1,cd1,rg1,cg1
        for rd,cd,rg,cg in zip(rd1,cd1,rg1,cg1): 
            # print rd,cd,rg,cg
            dist_displ= math.sqrt( rd**2. + cd**2. )
            dist_grad = math.sqrt( rg**2. + cg**2. )
            if dist_displ==0 or dist_grad==0: continue
            this_guy=(rd*rg+cd*cg)/(dist_displ*dist_grad)
            try: 
                follow_your_heart.append( math.acos(this_guy) )
            except:
                pass
            # print this_guy
            # follow_your_heart.append( this_guy ) 
    
    #print rowpos
    #print popav_row_direction
    # print follow_your_heart
    howmany_bins = 20
    max = np.pi
    bin_edges = np.linspace(0. , max  , howmany_bins)
    hist1, bin_edges1 = np.histogram( follow_your_heart, bins=bin_edges, density=True) #indiv cells hist
    # ax.plot(bin_edges1[:-1], hist1,label=filename.split('/')[-1] )#, color=colours[1])
    
    # hist1, bin_edges1 = np.histogram( pop_accuracy, bins=bin_edges, density=True) #pop measure histogram
    
    if '1cellfrom50' in filename: 
        colorindex=1
        hist_list_1from50.append(hist1)
    elif 'bc_g6_1c_f500' in filename:
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
if hist_list_1c:
    plot_mean_fillstd(bin_edges,hist_list_1c,color=colours[2],label='1c')
if hist_list_50c:
    plot_mean_fillstd(bin_edges,hist_list_50c,color=colours[0],label='50c')

# ax1.set_ylim([0,0.7])
ax.set_title('Pop measure, time of displ = '+str(time_of_dipl*10)+' MCS')
ax1.legend()
fig.savefig(figname)
plt.show()
sys.exit(1)

'''      
sys.exit(1)

sys.exit(1)
'''  
  
  