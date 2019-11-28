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

fig, (ax,ax1) = plt.subplots(nrows=2)

#fig, ax0 = plt.subplots()

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
      
  # Center of mass
  # row_cm = [np.mean(rows) for rows in rowpos]
  # col_cm = [np.mean(cols) for cols in colpos]
  # rowgrad_cm = [np.mean(rows) for rows in rowgrad]
  # colgrad_cm = [np.mean(cols) for cols in colgrad]
  
  l_correct_grad = []
  l_delta_rp=[]
  l_delta_cp=[]
  l_av_row_chemvec=[]
  l_av_col_chemvec=[]
  l_av_accuracy=[]
  for i,(rp,rg,cp,cg) in enumerate(zip (rowpos,rowgrad,colpos,colgrad)):
      
      # vector from cm to 250,0 is the perfect measurment of the gradient
      # the dot product of meansured gradient with the perfect direction is the quality of the measure
      
      delta_rp=[peak_row-x for x in rp]
      delta_cp=[peak_col-x for x in cp]    #because minus zero
      
      dist_from_peak = np.hypot(delta_rp,delta_cp) #this is a vector
      # dist_from_peak = delta_cp #this is a vector
      # print dist_from_peak
      # sys.exit(1)
      # my_check = np.mean( np.hypot(rowgrad,colgrad) ) #this was just to check that hypot of grad vector is 1 - it is!
      # print my_check
      
      delta_rp = np.divide(delta_rp,dist_from_peak) #normalising the "right direction vector"
      delta_cp = np.divide(delta_cp,dist_from_peak)
      
      #dot product between the right direction vector and the measured gradient
      #correct_grad = [ irg*idelta_rp + icg*idelta_cp  for irg,idelta_rp,icg,idelta_cp in zip(rg,delta_rp,cg,delta_cp)  ]
      #this is the absolute angle of the cells chem vector 
      # correct_grad = [ math.atan2( icg,irg)  for irg,idelta_rp,icg,idelta_cp in zip(rg,delta_rp,cg,delta_cp)  ]
      correct_grad = [  np.mod( math.atan2( icg,irg) - math.atan2( idelta_cp,idelta_rp)  , 2*np.pi) for irg,idelta_rp,icg,idelta_cp in zip(rg,delta_rp,cg,delta_cp)  ]
      l_correct_grad.append(correct_grad)
      
      l_delta_rp.append(delta_rp)
      l_delta_cp.append(delta_cp)
      
      #Distance for cells in contact with this one
      contact_distance = 12.
      period_lookingback=0 # zero is we look at only now
      l_av_row_chemvec_thistime=len(rowpos[i])*[0.]
      l_av_col_chemvec_thistime=len(rowpos[i])*[0.]
      l_av_accuracy_thistime=len(rowpos[i])*[0.]
      av_row_chemvec=0.
      av_col_chemvec=0.
      howmany_contact=0
      
      if True:
      # for period in range(period_lookingback+1): 
          # print "period :", period
          period=0
          if i-(period_lookingback+1)<0: continue
          for j,(r1,c1,rg1,cg1) in enumerate(zip(rowpos[i-period],colpos[i-period],rowgrad[i-period],colgrad[i-period])):
              av_accuracy=0
              howmany_contact=0
              delta_rp = peak_row - r1
              delta_cp = peak_col - c1
              # dist_from_peak = np.hypot(delta_rp,delta_cp)
              # delta_rp /= dist_from_peak    # now these are normalised directions
              # delta_cp /= dist_from_peak
              accuracy = (math.atan2(delta_cp, delta_rp) - math.atan2(cg1,rg1))%(2*np.pi) 
              av_accuracy += accuracy
              # print accuracy
              av_row_chemvec+=rg1
              av_col_chemvec+=cg1
              howmany_contact+=1
              
              
              for r2,c2,rg2,cg2 in zip(rowpos[i-period],colpos[i-period],rowgrad[i-period],colgrad[i-period]):
                  if r1==r2 and c1==c2: continue
                  if np.hypot(r1-r2,c1-c2) > contact_distance: continue
                  #those that remain are at a distance of less than contact_distance
                  #so we average their vectors, also with self
                  av_row_chemvec+=rg2
                  av_col_chemvec+=cg2
                  howmany_contact+=1
                  
                  delta_rp = peak_row - r2
                  delta_cp = peak_col - c2
                  accuracy = (math.atan2(delta_cp, delta_rp) - math.atan2(cg1,rg1))%(2*np.pi)
                  av_accuracy += accuracy
                
              av_row_chemvec/=float(howmany_contact)
              av_col_chemvec/=float(howmany_contact)
              l_av_row_chemvec_thistime[j] += av_row_chemvec
              l_av_col_chemvec_thistime[j] += av_col_chemvec
              
              av_accuracy/=float(howmany_contact)
              l_av_accuracy_thistime[j] += av_accuracy
      l_av_row_chemvec.append([x/float(period_lookingback+1) for x in l_av_row_chemvec_thistime])
      l_av_col_chemvec.append([x/float(period_lookingback+1) for x in l_av_col_chemvec_thistime])
      
      l_av_accuracy.append([x/float(period_lookingback+1) for x in l_av_accuracy_thistime])
  
  #sys.exit(1)
  # print l_correct_grad
  # print "Avrg. quality of measure =", np.mean([np.mean(x) for x in l_correct_grad])
  # sys.exit(1)
  
  howmany_bins = 100
  #histogram of correctness of chemvec, for one cell over time
  hist1, bin_edges1 = np.histogram([[correctness[0] for correctness in l_correct_grad]], bins =howmany_bins, density=True) 
  #histogram of all measures of all cells at one time point
  hist3, bin_edges3 = np.histogram([x for subs in l_correct_grad for x in subs] ,bins =howmany_bins, density=True) 
  ax.plot(bin_edges1[:-1], hist1,label="One cell, all times")
  ax.plot(bin_edges3[:-1], hist3,label="All cell, all times")
  
  # hist2, bin_edges2 = np.histogram( [x for y in l_av_row_chemvec for x in y] , bins =250, density=True) 
  # hist4, bin_edges4 = np.histogram( [x for y in l_av_col_chemvec for x in y] , bins =250, density=True) 
  # ax1.plot(bin_edges2[:-1], hist2,label='<chem vec row>')
  # ax1.plot(bin_edges4[:-1], hist4,label='<chem vec col>')
  
  bin_edges = np.linspace(0. , 2.*np.pi, howmany_bins)
  hist5, bin_edges5 = np.histogram( [x for y in l_av_accuracy for x in y] , bins =bin_edges, density=True) 
  print "howmany_bins=",howmany_bins,"len bins is ", len(bin_edges5)
  if len(l_tot_hist) !=0:
      l_tot_hist = [x+y for x,y in zip (l_tot_hist,hist5)]
  else: 
      l_tot_hist = hist5
  # ax1.plot(bin_edges5[:-1], hist5,label='<accuracy angle>')
  
  # ax1.plot([0,0],[0,2])
  
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
ax.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2., 2*np.pi])
ax.set_xticklabels(['$0$', r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$'])
ax.set_ylim([0,0.2])
ax.legend()

print bin_edges5[:-1]
print l_tot_hist
ax1.plot(bin_edges5[:-1], [x/float(len(sys.argv[5:])) for x in l_tot_hist],label='<accuracy angle>')
ax1.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2., 2*np.pi])
ax1.set_xticklabels(['$0$', r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$'])
ax1.legend()
fig.savefig(figname)
plt.show()
  
# print np.mean(rowgrad_cm),np.mean(colgrad_cm)
sys.exit(1)
  
  
  