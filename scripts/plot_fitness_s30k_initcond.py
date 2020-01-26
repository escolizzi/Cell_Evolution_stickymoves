#!/usr/bin/python2.7

'''
Plot fitness as box plot for different initial configurations, 
and different ratios of adhereing and repelling cells

'''

import sys,math,os,subprocess,random
#from PIL import Image
import matplotlib as mpl
import KL_adhesion_module as klam
# mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from matplotlib.lines import Line2D
from matplotlib import gridspec 
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

fig = plt.figure() 
gs = gridspec.GridSpec(1, 2, width_ratios=[20, 1]) 
ax = plt.subplot(gs[0])
ax1 = plt.subplot(gs[1],sharey=ax)

# fig, (ax,ax1) = plt.subplots(ncols=2, sharey=True)

#fig, ax0 = plt.subplots()

if len(sys.argv) <5:
  print "This is the program 'plot_fitness_30k....py'"
  print "Usage: ./plot_fitness_....py <output figure name> <cell targetarea> <peak_row> <peak_col> <fitness file1> <fitness file2> ..."
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
l_tot_boxplot=[]
l_boxplot_pos=[]
l_boxplot_ticklabels=[]
i=0
for filename in sys.argv[5:]:
  rowpos1=[]  #this is row pos for one of the two types
  colpos1=[]  
  rowpos2=[]  #this is for the other
  colpos2=[]  
  
  print filename
  # avangle=[0.]*200*int(binmult)   #1 bin is binsize celldiam 
  # binamount=[0]*200*int(binmult)  #to keep track of how many in each bin
  # angles=[[] for x in range(200*int(binmult))]
  
  #########################################################
  ## -- read data from file: store raw cell positions -- ##
  #########################################################
  
  output = subprocess.Popen(['tail', '-1', filename], stdout=subprocess.PIPE).communicate()[0]
  outsplit = output.split(' ')
  key1 = outsplit[4]
  lock1 = outsplit[5][:-1] #to remove the '\n' character
  timepoint = int(output.split(' ')[0])
  
  print "key1", key1
  print "lock1", lock1
  # sys.exit(1)
  #read file
  with open(filename,"r") as fin:
    first_line=True
    for line in fin:
      
      # Because first line is sigma==0, for some reason
      if first_line: 
          first_line=False
          continue
      
      line=line.split()
      if int(line[0])!=timepoint:
        break
      # print line[4],line[5]
      # continue
      # print line[7],line[8]
      if line[4]==key1 and line[5]==lock1:
          # print "hello1"
          rowpos1.append(float(line[2]))
          colpos1.append(float(line[3]))
      else: 
          # print key1,lock1
          key2=line[4]
          lock2= line[5]
          rowpos2.append(float(line[2]))
          colpos2.append(float(line[3]))
      
  # Center of mass
  # row_cm = [np.mean(rows) for rows in rowpos]
  # col_cm = [np.mean(cols) for cols in colpos]
  # rowgrad_cm = [np.mean(rows) for rows in rowgrad]
  # colgrad_cm = [np.mean(cols) for cols in colgrad]
  
  lfitness=[]
  dist_from_peak1=[ np.hypot( peak_row-x , peak_col-y)  for x,y in zip(rowpos1,colpos1)]
  dist_from_peak2=[ np.hypot( peak_row-x , peak_col-y)  for x,y in zip(rowpos2,colpos2)]
  key1=[int(x) for x in key1]
  lock1=[int(x) for x in lock1]
  key2=[int(x) for x in key2]
  lock2=[int(x) for x in lock2]
  
  gamma1 = klam.JWithMedium(key1,klam.lookuptable_Jmedium) - 0.5*klam.JWithOtherTau(( key1,lock1 ),( key1,lock1 ))
  gamma2 = klam.JWithMedium(key2,klam.lookuptable_Jmedium) - 0.5*klam.JWithOtherTau(( key2,lock2 ),( key2,lock2 ))
  
  l_tot_boxplot.append(dist_from_peak1)
  l_boxplot_ticklabels.append(filename+'_g'+str(gamma1))
  l_boxplot_pos.append(i-0.5)
  
  l_tot_boxplot.append(dist_from_peak2)
  l_boxplot_ticklabels.append(filename+'_g'+str(gamma2))
  l_boxplot_pos.append(i+0.5)
  i+=3
  # print gamma1
  # print gamma2
  # 
  # bla = filename.split('_')[-2]
  # position = int(bla[0])/float( int(bla[0])+int(bla[-1]) )
  # print position
  # howmanyfiles = len(sys.argv[5:])
  # # ax.boxplot([dist_from_peak1,dist_from_peak2], positions=[position+0.05*howmanyfiles/3,position-0.5*howmanyfiles/3,])
  # print "i =",i
  # 
  # 
  # ax.boxplot([dist_from_peak1,dist_from_peak2], positions = [i-0.1,i+0.1])
  # labels = [item.get_text() for item in ax.get_xticklabels()]
  # labels[0]=filename
  # labels[1]=filename
  # 
  # print "labels now",labels
  # # xtickslabels=[]
  # # for x in sys.argv[5:]:
  # #     xtickslabels.append(x)
  # #     xtickslabels.append(x)
  # # ax.set_xticklabels(xtickslabels,rotation=45)
  # i+=1
ax.boxplot(l_tot_boxplot,positions=l_boxplot_pos)
ax.set_xlim([-2,i+2])
# labels = [item.get_text() for item in ax.get_xticklabels()]
# if you wish to explicitly set tick labels
# ax.set_xticklabels( len(labels)*['www'])
# print "labels final",labels

# ax1 should be a plot of the fitness function
h_dist=50. #the_line in cell_evolution

X=range(0,500,1)
Y=[1. / ( 1. + pow( dist/h_dist , 2.) ) for dist in X]
#Fitness = 1. / ( 1. + pow( dist/h_dist , 2.) );
ax1.plot(Y,X) #notice that it is reversed
ax.set_xticklabels(l_boxplot_ticklabels,rotation=45)
plt.show()
sys.exit(1)

  