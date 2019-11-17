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
# fig, (ax0,ax1) = plt.subplots(nrows=2,sharex=True)
# fig.subplots_adjust(hspace=0)

#fig, ax0 = plt.subplots()

if len(sys.argv) <3:
  print "This is the program 'plot_flowfield_fromCM.py'"
  print "Usage: ./plot_flowfield_fromCM.py <figure name> <cell targetarea> <binsize [pix]> <filename> <filename> <filename> ..."
  sys.exit(1)
else:
  figname=sys.argv[1]
  cell_tarsize=float(sys.argv[2])
  binsize=int(sys.argv[3])
  
macoldist=0
filecounter=0


field_size = 2*500


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
  row_cm = [np.mean(rows) for rows in rowpos]
  col_cm = [np.mean(cols) for cols in colpos]
  rowgrad_cm = [np.mean(rows) for rows in rowgrad]
  colgrad_cm = [np.mean(cols) for cols in colgrad]
  
  l_correct_grad = []
  l_delta_rp=[]
  l_delta_cp=[]
  for rp,rg,cp,cg in zip (rowpos,rowgrad,colpos,colgrad):
      # vector from cm to 250,0 is the perfect measurment of the gradient
      # the dot product of meansured gradient with the perfect direction is the quality of the measure
      delta_rp=[x-250 for x in rp]
      # print rp,delta_rp
      delta_cp=cp[:]    #because minus zero
      dist_from_peak = np.hypot(delta_rp,delta_cp)
      # my_check = np.mean( np.hypot(rowgrad,colgrad) ) #this was just to check that hypot of grad vector is 1 - it is!
      # print my_check
      
      delta_rp = [x/dist_from_peak for x in delta_rp] #normalising the "right direction vector"
      delta_cp = [x/dist_from_peak for x in delta_cp]
      
      print np.hypot(delta_cp,delta_rp)
      #dot product between the right direction vector and the measured gradient
      correct_grad = [ irg*idelta_rp + icg*idelta_cp  for irg,idelta_rp,icg,idelta_cp in zip(rg,delta_rp,cg,delta_cp)  ]
      # print len(correct_grad)
      l_correct_grad.append(correct_grad)
      
      l_delta_rp.append(delta_rp)
      l_delta_cp.append(delta_cp)
  
  print l_correct_grad
  print "Avrg. quality of measure =", np.mean([np.mean(x) for x in l_correct_grad])
  # sys.exit(1)
  
  hist1, bin_edges1 = np.histogram([[correctness[0] for correctness in l_correct_grad]])
  # hist2, bin_edges2 = np.histogram(colgrad_cm)
  
  fig, ax = plt.subplots()
  
  # THERE IS SOMETHING ODD HERE ... IF YOU PRINT THE FOLOWING IT SHOULD BE ONE LIST OF FLOATS
  
  
  # ax.plot([x for x in range(len(l_correct_grad))], [np.mean(x) for x in l_correct_grad])
  # print [x[0] for x in l_correct_grad] 
  
  ax.plot([x for x in range(len(l_correct_grad))], [x[0] for x in l_correct_grad])
  
  ax.plot([0, len(l_correct_grad) ], [0,0])
  ax.plot([0, len(l_correct_grad) ], 2*[np.mean([np.mean(x) for x in l_correct_grad])])
  
  # ax.plot(bin_edges1[:-1], hist1)
  
  # ax.plot(bin_edges2[:-1], hist2)
  # fig.savefig(figname)
  plt.show()
  
  print np.mean(rowgrad_cm),np.mean(colgrad_cm)
  sys.exit(1)
  
  
  