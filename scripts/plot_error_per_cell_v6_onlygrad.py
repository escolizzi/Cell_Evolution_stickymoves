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



#########################
###                   ###
### ---   BEGIN   --- ###
###                   ###
#########################

colours=["firebrick","royalblue", "darkgoldenrod", "green", "salmon", "lightskyblue","orchid"]
filename=""
# fig, (ax,ax1) = plt.subplots(nrows=2,sharex=True)
# fig.subplots_adjust(hspace=0)

# fig, (ax,ax1) = plt.subplots(nrows=2)

fig, ax = plt.subplots()

if len(sys.argv) <5:
  print "This is the program 'plot_error_....py'"
  print "Usage: ./plot_error_....py <output figure name> <cell targetarea> <peak_row> <peak_col> <filename> <filename> <filename> ..."
  sys.exit(1)
else:
  figname=sys.argv[1]
  cell_tarsize=float(sys.argv[2])
  peak_row=int(sys.argv[3])
  peak_col=int(sys.argv[4])

# print "# WARNING: " 
# print "Warning, distance uses only column, for linear gradient"
# print "# WARNING: "
macoldist=0
filecounter=0
l_tot_hist=[]
for filename in sys.argv[5:]:
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

    print filename
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
    
    skip=500
    howlong=2000
    #we skip beginning and end
    rowpos=rowpos[skip:skip+howlong]
    colpos=colpos[skip:skip+howlong]
    rowgrad=rowgrad[skip:skip+howlong]
    colgrad=colgrad[skip:skip+howlong]
    
    l_cell_accuracy=[]
    
    # i indicises time
    for i,(rp,rg,cp,cg) in enumerate( zip(rowpos,rowgrad,colpos,colgrad) ):
        #Direction to the peak of the gradient from center of mass of blob
        #  (needed at printing time to calculate accuracy)
        for r1,c1,rg1,cg1 in zip(rp,cp,rg,cg):
            peakdir_row = (peak_row - r1) 
            peakdir_col = (peak_col - c1)
            dist_from_peak = np.hypot( peakdir_row, peakdir_col ) 
            
            gr_length = np.hypot( rg1,cg1 )
            dp = (peakdir_row*rg1 + peakdir_col*cg1)/float(dist_from_peak*gr_length)
            
            l_cell_accuracy.append( math.acos(dp) )
        
    #DONE WITH DATA ANALYSIS, NOW PRINTING
    
    howmany_bins = 30
    bin_edges = np.linspace(0. , np.pi, howmany_bins)
    hist1, bin_edges1 = np.histogram( l_cell_accuracy, bins=bin_edges, density=True)
    ax.plot(bin_edges1[:-1], hist1,label=filename.split('/')[-1])
    
ax.set_xticks([0, np.pi/4., np.pi/2., 3*np.pi/4., np.pi])
ax.set_xticklabels(['$0$', r'$\frac{\pi}{4}$', r'$\frac{\pi}{2}$', r'$\frac{3\pi}{4}$', r'$\pi$'])
ax.set_ylim([0,0.5])
ax.legend()
ax.set_title('Accuracy of gradient measure')
# ax1.legend()
plt.margins(0,0)
fig.savefig(figname)
plt.show()
sys.exit(1)

'''      
    sys.exit(1)
    
    
    #histogram of correctness of chemvec, for one cell over time
    # hist1, bin_edges1 = np.histogram([[correctness[0] for correctness in l_correct_grad]], bins =howmany_bins, density=True) 
    # ax.plot(bin_edges1[:-1], hist1,label="One cell, all times")
    # 

    # bla=sorted([x for subs in l_correct_grad for x in subs])
    # # ax.plot( range(len(bla)),bla )
    # ax.plot( range(len(bla)-1), [ y2-y1 for y1,y2 in zip(bla[:-1],bla[1:]) ] )
    # 
    #histogram of all measures of all cells at one time point
    # bin_edges = np.linspace(-2*np.pi , 2.*np.pi, howmany_bins)
    bin_edges = np.linspace(0. , np.pi, howmany_bins)
    #bin_edges = np.linspace(-1, 1, howmany_bins)
    hist3, bin_edges3 = np.histogram([x for subs in l_correct_grad for x in subs] ,bins =bin_edges, density=True) 
    ax.plot(bin_edges3[:-1], hist3,label="All cell, all times")


    acos_dotprod_neigh = [] 
    for pr,pc,avr,avc in zip(rowpos[skip:-skip],colpos[skip:-skip],l_av_row_chemvec,l_av_col_chemvec):
      for i_p_r,i_p_c,i_av_r,i_av_c in zip(pr,pc,avr,avc):
          peakdist=np.hypot(peak_row-i_p_r, peak_col-i_p_c)
          delta_rp=(peak_row-i_p_r)/peakdist
          delta_cp=(peak_col-i_p_c)/peakdist
          dot_prod_neigh = (delta_rp*i_av_r+delta_cp*i_av_c)/np.hypot(i_av_r,  i_av_c)
          i_acos_dotprod_neigh = math.acos(dot_prod_neigh)
          acos_dotprod_neigh.append(i_acos_dotprod_neigh)
          
    # acos_dotprod_neigh = [ math.acos( ( ( (peak_row-i_p_r)/(np.hypot(peak_row-i_p_r, peak_col-i_p_c) ) )*i_av_r + ((peak_col-i_p_c)/(np.hypot(peak_row-i_p_r, peak_col-i_p_c))*i_av_c)/np.hypot(i_av_r,  i_av_c) ) for pr,pc,avr,avc in zip(rowpos[skip:-skip],colpos[skip:-skip],l_av_row_chemvec,l_av_col_chemvec) for i_p_r,i_p_c,i_av_r,i_av_c in zip(pr,pc,avr,avc) ]
    hist4, bin_edges4 = np.histogram(acos_dotprod_neigh ,bins =bin_edges, density=True) 
    ax.plot(bin_edges4[:-1], hist4,label="All average neigh, all times")
    # 

    # #histogram of average angle each time step
    # hist1, bin_edges1 = np.histogram([np.mean(x) for x in l_correct_grad] ,bins =howmany_bins, density=True) 
    # ax.plot(bin_edges1[:-1], hist1,label="Average angle, per time step")
    # 

    # hist2, bin_edges2 = np.histogram( [x for y in l_av_row_chemvec for x in y] , bins =250, density=True) 
    # hist4, bin_edges4 = np.histogram( [x for y in l_av_col_chemvec for x in y] , bins =250, density=True) 
    # ax1.plot(bin_edges2[:-1], hist2,label='<chem vec row>')
    # ax1.plot(bin_edges4[:-1], hist4,label='<chem vec col>')

    hist5, bin_edges5 = np.histogram( [np.mean(x) for x in rowgrad[skip:-skip]],bins=howmany_bins, density=True) 
    hist6, bin_edges6 = np.histogram( [np.mean(x) for x in colgrad[skip:-skip]],bins=howmany_bins, density=True) 
    ax1.plot(bin_edges5[:-1], hist5,label='Avrg. row/time')
    ax1.plot(bin_edges6[:-1], hist6,label='Avrg. col/time')

    pop_measure_row = [np.mean(x)+np.mean(y) for x,y in zip(rowgrad[skip:-skip],l_sum_row_chemvec)]
    pop_measure_col = [np.mean(x)+np.mean(y) for x,y in zip(colgrad[skip:-skip],l_sum_col_chemvec)]
    hist7, bin_edges7 = np.histogram(pop_measure_row ,bins=howmany_bins, density=True) 
    hist8, bin_edges8 = np.histogram(pop_measure_col ,bins=howmany_bins, density=True) 
    ax1.plot(bin_edges7[:-1], hist7,label='full, row')
    ax1.plot(bin_edges8[:-1], hist8,label='full, col')

    #And now we can get the accuracy of the measure:
    acos_dotprod_pop = [] 
    for pr,pc,mr,mc in zip(rowpos[skip:-skip],colpos[skip:-skip],pop_measure_row,pop_measure_col):
      av_pos_row=np.mean(pr)
      av_pos_col=np.mean(pc)
      peakdist=np.hypot(peak_row - av_pos_row, peak_col - av_pos_col)
      delta_rp=(peak_row - av_pos_row)/peakdist
      delta_cp=(peak_col - av_pos_col)/peakdist
      dotprod_pop = (delta_rp*mr+delta_cp*mc)/np.hypot(mr, mc)
      i_acos_dotprod_pop = math.acos(dotprod_pop)
      acos_dotprod_pop.append(i_acos_dotprod_pop)
    hist9, bin_edges9 = np.histogram(acos_dotprod_pop ,bins =bin_edges, density=True) 
    ax.plot(bin_edges9[:-1], hist9,label="Full Pop")

    print "Avrg. pop row,col",np.mean(pop_measure_col), np.mean(pop_measure_col)

    # bin_edges = np.linspace(0. , 2.*np.pi, howmany_bins)
    # hist5, bin_edges5 = np.histogram( [x for y in l_av_accuracy for x in y] , bins =bin_edges, density=True) 
    # print "howmany_bins=",howmany_bins,"len bins is ", len(bin_edges5)
    # if len(l_tot_hist) !=0:
    #     l_tot_hist = [x+y for x,y in zip (l_tot_hist,hist5)]
    # else: 
    #     l_tot_hist = hist5
    # ax1.plot(bin_edges5[:-1], hist5,label='<accuracy angle>')

    ax1.plot([0,0],[0,4])

    # for i in range(len(l_correct_grad)/2,len(l_correct_grad), len(l_correct_grad)/10):
    #     hist2, bin_edges2 = np.histogram(l_correct_grad[i] ,bins =50,density=True) #histogram of all measures of all cells at one time point
    #     ax1.plot(bin_edges2[:-1], hist2)
    #     # break
    # 
    # fig, ax = plt.subplots()

    # THERE IS SOMETHING ODD HERE ... IF YOU PRINT THE FOLOWING IT SHOULD BE ONE LIST OF FLOATS


    # ax.plot([x for x in range(len(l_correct_grad))], [np.mean(x) for x in l_correct_grad])
    # print [x[0] for x in l_correct_grad] 

    # ax.plot([x for x in range(len(l_correct_grad))], [x[0] for x in l_correct_grad])
    # 
    # ax.plot([0, len(l_correct_grad) ], [0,0])
    # ax.plot([0, len(l_correct_grad) ], 2*[np.mean([np.mean(x) for x in l_correct_grad])])
    # 

# fig.savefig(figname)

# ax.set_xticks([0, np.pi/4., np.pi/2., 3*np.pi/4., np.pi])
# ax.set_xticklabels(['$0$', r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$'])
ax.set_ylim([0,1.0])
ax.legend()

# print bin_edges5[:-1]
# print l_tot_hist
# ax1.plot(bin_edges5[:-1], [x/float(len(sys.argv[5:])) for x in l_tot_hist],label='<accuracy angle>')
# ax1.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2., 2*np.pi])
# ax1.set_xticklabels(['$0$', r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$'])
ax1.legend()
fig.savefig(figname)
plt.show()
  
# print np.mean(rowgrad_cm),np.mean(colgrad_cm)
sys.exit(1)
'''  
  
  