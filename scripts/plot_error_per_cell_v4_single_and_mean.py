#!/usr/bin/python2.7

'''
Plot the error a cell makes in detecting gradient at different time scales, both average and mean square error?
Compare to direction???

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
    
    skip=1000
    howlong=2000
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
    
    # i indicises time
    for i,(rp,rg,cp,cg) in enumerate( zip(rowpos,rowgrad,colpos,colgrad) ):
        #Direction to the peak of the gradient from center of mass of blob
        #  (needed at printing time to calculate accuracy)\
        
        cm_row = np.mean(rp) 
        cm_col = np.mean(cp)
                
        dist_from_peak = math.sqrt( (peak_row - cm_row)**2. + (peak_col - cm_col)**2. )
        peakdir_row.append( (peak_row - cm_row)/dist_from_peak )
        peakdir_col.append( (peak_col - cm_col)/dist_from_peak )
        #find average chem vector of the population at this time point 
        av_rg = np.mean(rg)
        av_cg = np.mean(cg)
        
        popav_row_direction.append( av_rg )
        popav_col_direction.append( av_cg )
    # ACCURACY
    pop_accuracy = []
    for (pr,pc,ar,ac)  in zip( peakdir_row,peakdir_col, popav_row_direction,popav_col_direction ) :
        x = (pr*ar+pc*ac)/(math.sqrt( ar**2. + ac**2. )) 
        # x = (pr*ar+pc*ac)/(math.sqrt( ar**2. + ac**2. )*math.sqrt(pr**2. + pc**2.)) 
        #if x>1.: x=1.
        #elif x<-1.: x=-1.
        # elif x < -0.9:
        #     print "weird x=",x
        try:
            pop_accuracy.append(math.acos(x))
        except:
            print "This is not working, x =",x
        # pop_accuracy.append(x)
    #individual cell accuracy, we just flatten all the data we have
    print
    indiv_cell_accuracy=[]
    for (rp,rg,cp,cg) in zip(rowpos,rowgrad,colpos,colgrad): 
        for (rp1,rg1,cp1,cg1) in zip(rp,rg,cp,cg):
            peak_dist_row1 = peak_row - rp1
            peak_dist_col1 = peak_col - cp1
            peak_dist1 = math.sqrt( peak_dist_row1**2. + peak_dist_col1**2.)
            
            x = (rg1*peak_dist_row1+cg1*peak_dist_col1)/( peak_dist1 *math.sqrt(rg1**2. + cg1**2.)) 
            try:
                indiv_cell_accuracy.append(math.acos(x))
            except:
                print "This is not working, x =",x
    #DONE WITH DATA ANALYSIS, NOW PRINTING
    
    #print rowpos
    #print popav_row_direction
    
    howmany_bins = 20
    bin_edges = np.linspace(0. , np.pi, howmany_bins)
    hist2, bin_edges2 = np.histogram( indiv_cell_accuracy, bins=bin_edges, density=True) #indiv cells hist
    hist1, bin_edges1 = np.histogram( pop_accuracy, bins=bin_edges, density=True) #pop measure histogram
    
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
        ax.plot(bin_edges1[:-1], hist2,label=filename.split('/')[-1], color=colours[colorindex])
        ax1.plot(bin_edges1[:-1], hist1,label=filename.split('/')[-1],color=colours[colorindex])
        not_recognised=False
    
ax.set_xticks([0, np.pi/4., np.pi/2., 3*np.pi/4., np.pi])
ax.set_xticklabels(['$0$', r'$\frac{\pi}{4}$', r'$\frac{\pi}{2}$', r'$\frac{3\pi}{4}$', r'$\pi$'])
# ax.set_ylim([0,1.])
ax.set_title('Indiv measure')
ax.legend()

if hist_list_1from50:
    plot_mean_fillstd(bin_edges,hist_list_1from50,color=colours[1],label='1from50')
if hist_list_1c:
    plot_mean_fillstd(bin_edges,hist_list_1c,color=colours[2],label='1c')
if hist_list_50c:
    plot_mean_fillstd(bin_edges,hist_list_50c,color=colours[0],label='50c')
    
ax1.set_xticks([0, np.pi/4., np.pi/2., 3*np.pi/4., np.pi])
ax1.set_xticklabels(['$0$', r'$\frac{\pi}{4}$', r'$\frac{\pi}{2}$', r'$\frac{3\pi}{4}$', r'$\pi$'])
ax1.set_ylim([0,0.7])
ax1.set_title('Pop measure')
ax1.legend()
fig.savefig(figname)
plt.show()
sys.exit(1)

'''      
sys.exit(1)

sys.exit(1)
'''  
  
  