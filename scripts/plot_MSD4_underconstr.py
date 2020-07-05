#!/usr/bin/python2.7

'''
Plot the MSD over time for one or more simulations.
Plot some selected cell tracks, centered at 0,0
'''

import sys,math,os,subprocess,random
#from PIL import Image
import matplotlib as mpl
# mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from matplotlib.lines import Line2D
from scipy import optimize
from scipy import interpolate
from scipy.signal import savgol_filter
from scipy.ndimage.filters import uniform_filter1d
import numpy as np

#manager = plt.get_current_fig_manager()
#print manager
#sys.exit(1)
#########################
###                   ###
### ---   BEGIN   --- ###
###                   ###
#########################

def piecewise_linear(x, x0, y0, k1, k2):
  return np.piecewise(x, [x < x0], [lambda x:k1*x + y0-k1*x0, lambda x:k2*x + y0-k2*x0])

colours=[(0/255.,151/255.,255/255.),(240/255.,218/255.,104/255.), "darkgoldenrod", "green", "salmon", "lightskyblue","orchid"]
filename=""
fig, (ax0, ax1,ax2) = plt.subplots(ncols=3)
#fig, ax0 = plt.subplots()

if len(sys.argv) <3:
  print "This is the program 'plot_MSD2.py'"
  print "Usage: ./plot_MSD.py <figure name> <nr of tracks to plot> <filename(s)> "
  sys.exit(1)
else:
  figname=sys.argv[1]
  tracknr=int(sys.argv[2])

filecounter=0
for filename in sys.argv[3:]:

  print "reading file ",filename
  timepoints=[]
  MSD=[]
  SDEV=[]
  xpos=[]
  ypos=[]
  xn=[]
  xpos.append([])
  ypos.append([])
  count=0
  
  ##read data from file: store raw cell positions
  with open(filename,"r") as fin:

    #read first line first to get the first time point (could probably more clever, but hey)
    for line in fin:
      line=line.split()
      timepoints.append(int(line[0]))
      xpos[-1].append(float(line[3]))
      ypos[-1].append(float(line[4]))
      count+=1
      break

    #read rest of file
    for line in fin:
      line=line.split()
      if int(line[0])!=timepoints[-1]:
        xpos.append([-1.]*len(xpos[0]))
        ypos.append([-1.]*len(xpos[0]))
        count=0
        timepoints.append(int(line[0]))
      if (timepoints[-1]==timepoints[0]):
        xpos[-1].append(float(line[3]))
        ypos[-1].append(float(line[4]))
      else:
        xpos[-1][int(line[1])-1]=float(line[3])
        ypos[-1][int(line[1])-1]=float(line[4])
      count+=1
  
  #I should eliminate the initial data, because it is initial condition
  toremove=len(timepoints)/9
  # step=5
  xpos = xpos[toremove:]
  # xpos = xpos[::step]
  ypos = ypos[toremove:]
  # ypos = ypos[::step]
  timepoints=timepoints[: -toremove ] #this has to be like this because timepoints is also used as x axis
  # timepoints=timepoints[::step]
  
  print "Done reading, now calculating MSD"
  ##calculate MSD and standard error
  # We take non overlapping data so that we minimise correlations
  
  #data points are collected every so many time steps:
  timeresolution=timepoints[1]-timepoints[0]
  print "timeresolution", timeresolution 
  #we use this to convert between time points and indexes in xpos and ypos
  len_xpos=len(xpos)
  nrcells=len(xpos[0])
  count_toremove_fromtimepoints=0
  #we exclude the first time point because it is zero which would give a timeinterval and step size of zero
  # we axclude later time points, where we do not have much data (much is at least 5...not so much)
  max_interval=timepoints[-1]/int(timeresolution)
  min_interval=50/timeresolution    #50 is how often we change persistence
                                    # THIS HAS TO BE AT LEAST 1
  all_sd=[]
  for timeinterval in timepoints[min_interval: max_interval]:
      # print "timeinterval",timeinterval
      xn=[] #for standard deviation
      MSD.append(0.0)
      sd=0.0
      count2=0  #just a counter in case some cells have died
      stepsize=timeinterval/timeresolution
      #we use i to iterate over xpos indexes
      for i in range( 0,len_xpos, stepsize/2 ):
          for c in range(nrcells):
              try:
                  a=xpos[i+stepsize][c]
              except:
                  #print "Exceeding stepsize: ", i+stepsize
                  #print "for cell c:",c
                  break
              if xpos[i+stepsize][c]>-1:
                  sq_displ = (xpos[i+stepsize][c]-xpos[i][c])**2. + (ypos[i+stepsize][c]-ypos[i][c])**2.
                  xn.append(sq_displ)   #we append this here and below we use it for stdev
                  MSD[-1]+=xn[-1]
                  count2+=1
              else:
                  print "cell ",c,"is no more"
      
      MSD[-1]=MSD[-1]/float(count2)
      all_sd.append(MSD[-1])
      
      for el in xn:
        sd+=(el-MSD[-1])**2.
      SDEV.append( math.sqrt(sd/float(count2)) )   # if divided by sqrt(count2) -> std error of the mean, else stdev
  
  ax0.plot(timepoints[min_interval:max_interval],MSD,'o',c=colours[filecounter], markersize=0.5)
  ax0.errorbar(timepoints[min_interval:max_interval],MSD, yerr=SDEV, fmt='', c=colours[filecounter],errorevery=500)
  # print all_sd[::500],timepoints[min_interval:max_interval:500]
  # ax0.boxplot(all_sd[filecounter*250::500], positions=timepoints[min_interval+filecounter*250:max_interval:500], widths=500,showfliers=False)
  print "Error bars standard error of the mean = standard dev / sqrt(how many data points) "
  print "howmany data points = ", len(xpos[0])
  log_t=timepoints[min_interval:max_interval]
  log_MSD=MSD[:]
  if log_t[0]==0:
      log_t=log_t[1:]
      log_MSD=log_MSD[1:]
  
  log_t = [np.log(x) for x in log_t]
  log_MSD = [np.log(x) for x in log_MSD]
  
  #PLOTS LOG LOG DATA
  #print "Uncomment this if you want scatter plot of loglog data"
  
  ax1.plot(log_t,log_MSD, 'o',c=colours[filecounter], markersize=0.5) 
  # ~ ax1.plot(log_t,log_MSD, 'o',c=colours[filecounter], markersize=0.5) 
  error_bars_log_plus = [ x+y/z for x,y,z in zip(log_MSD,SDEV,MSD)]
  error_bars_log_minus = [ x-y/z for x,y,z in zip(log_MSD,SDEV,MSD)]
  ax1.fill_between(log_t,error_bars_log_minus,error_bars_log_plus, color=colours[filecounter], alpha=0.5) 
  
  ax2.plot(log_t,log_MSD, 'o',c=colours[filecounter], markersize=0.5) 
  # ~ ax1.plot(log_t,error_bars_log_minus, c=colours[filecounter], alpha=0.5) 
  
  # ~ ax1.plot(timepoints[min_interval:max_interval],MSD,'o',c=colours[filecounter], markersize=0.5)
  # ~ ax1.set_xscale('log')
  # ~ ax1.set_yscale('log')
  # ~ ax1.plot(log_t,[8 for _ in log_t], 'o',c=colours[filecounter], markersize=0.5) 
  # ~ ax1.plot(log_t,[10 for _ in log_t], 'o',c=colours[filecounter], markersize=0.5) 
  # ~ ax1.plot([8 for _ in log_MSD], log_MSD,c=colours[filecounter], markersize=0.5) 
  # ~ ax1.plot([10 for _ in log_MSD], log_MSD,c=colours[filecounter], markersize=0.5)
  # ~ print "hello"
  # ~ print "hello"
  # ~ print "hello"
  #ax1.plot(log_t[:len(log_t)/100],log_MSD[:len(log_t)/100], 'o',c=colours[filecounter], markersize=0.5) 
  
  #ax1.plot(log_t[:-1],[(x-y)/(tx-ty) for x,y,tx,ty in zip(log_MSD[1:],log_MSD[:-1],log_t[1:],log_t[:-1])], 'o',c=colours[filecounter], markersize=0.5) 
  #running average
  
  #mmm = uniform_filter1d( [(x-y)/(tx-ty) for x,y,tx,ty in zip(log_MSD[1:],log_MSD[:-1],log_t[1:],log_t[:-1])], size=2)
  #ax1.plot(log_t[:len(log_t)/100],mmm[:len(log_t)/100])
  
  #mmm = uniform_filter1d( [(x-y)/(tx-ty) for x,y,tx,ty in zip(log_MSD[1:],log_MSD[:-1],log_t[1:],log_t[:-1])], size=20)
  #ax1.plot(log_t[:-1],mmm)
  
  #mmm = uniform_filter1d( [(x-y)/(tx-ty) for x,y,tx,ty in zip(log_MSD[1:],log_MSD[:-1],log_t[1:],log_t[:-1])], size=500)
  #ax1.plot(log_t[:-1],mmm)
  
  #Idea print running average with logarithmically increasing size: 1->2, 2->4, 4 -> 8 ,...
  #print log_t
  l_run_avr=[]
  maxx=11
  step=56
  stepsize=step/float(maxx)
  xaxis=np.linspace(0, 11, step )
  for i in xaxis:
    #print "i=",i
    bla=0
    count =0
    for t,x in zip(log_t,log_MSD):
      if t<i and t>=i-stepsize:
        bla+=x
        count += 1
      else:
        #print i,bla,count
        l_run_avr.append(bla/float(count) if count>0 else 0.)
        break
  print l_run_avr
  #sys.exit(1)
  #print "length xaxis, lrunavr", len(xaxis), len(l_run_avr)
  #ax1.scatter(xaxis[:-1] , [ x-y for x,y in zip(l_run_avr[1:],l_run_avr[:-1]) ])
  #ax1.scatter(xaxis , l_run_avr)
  #plt.show()
  #sys.exit(1)
  
  #Plot normal, but in log log scale
  # ax1.plot(timepoints[min_interval:max_interval],MSD,'o',c=colours[filecounter], markersize=0.5)
  # ax0.plot( [timepoints[min_interval], timepoints[max_interval-1]],[MSD[0], MSD[0] -timepoints[min_interval] + timepoints[max_interval-1]],'-')
  # ax1.plot( [timepoints[min_interval], timepoints[max_interval-1]],[MSD[0], MSD[0] -timepoints[min_interval] + timepoints[max_interval-1]],'-')
  # ax1.set_xscale('log')
  # ax1.set_yscale('log')
  # continue
  
  #this interpolates two lines - doesn't work very well
  if False:
      #to get two lines interpolation
      p , e = optimize.curve_fit(piecewise_linear, log_t, log_MSD)
      #xd = np.linspace(0, 15, 100)    
      ax1.plot(log_t, piecewise_linear(log_t, *p), '--',c=colours[filecounter])
      #ax1.plot(log_t, poly1d_fn(log_t), '--',c=colours[filecounter])
  
  
  #Linear interpolation of log log transformed data, 
  # to get exponent
  elif False:
      m_and_b = np.polyfit(log_t, log_MSD, 1)
      poly1d_fn = np.poly1d(m_and_b) # CAREFUL HERE: poly1d_fn is now a function which takes in x and returns an estimate for y
      print "filename: ", filename, "m,b:", m_and_b
      ax1.plot(log_t, poly1d_fn(log_t), '--',c=colours[filecounter])
  #interpolates with cuve 80
  elif True:
      # ltime=timepoints[min_interval:max_interval]
      # tck1=interpolate.splrep(ltime, MSD, s=0.0)
      # ynew1 = interpolate.splev(ltime, tck1, der=0)
      yhat = savgol_filter(log_MSD, 5, 3) # window size, polynomial order
      # tck2=interpolate.splrep(log_t, yhat, s=1.0)
      tck2=interpolate.splrep(log_t, yhat, s=0.)
      # xnew = np.arange(4,10,0.05)
      xnew = np.arange(log_t[0],log_t[-1],0.05)
      ynew2 = interpolate.splev(xnew, tck2, der=0)
      # ax1.plot(xnew,ynew2)
      
      yhat = savgol_filter(ynew2, 21, 2) # window size 51, polynomial order 3
      print "Hello"
      tck2=interpolate.splrep(xnew, yhat, s=0.1)
      # xnew = np.arange(4,10,0.1)
      ynew2 = interpolate.splev(xnew, tck2, der=0)
      ax2.plot(xnew,ynew2)
      
      xnew = np.arange(log_t[0],log_t[-1],0.05)
      ynew2 = interpolate.splev(xnew, tck2, der=1)
      
      
      # yhat2 = savgol_filter(ynew2, 1001, 2) # window size 51, polynomial order 3
      # 
      # ax0.plot(ltime,yhat)
      # ax0.plot(ltime,ynew2)
      print tck2
      ax2.plot(xnew,ynew2)
      # timeresolution=float(timepoints[1]-timepoints[0])
      # hand_made_derivative = [(y1-y0)/timeresolution for y1,y0 in zip(MSD[1:],MSD[:-1] )]
      # ax1.plot(ltime[1:] , hand_made_derivative,label='hello' )
      # 
      
      #I should maybe interpolate the original function
      # ax1.plot(log_t[5:-5],[np.log(x) for x in yhat[5:-5]])
      # ax1.plot(timepoints[min_interval:max_interval],MSD)
      # ax1.plot(timepoints[min_interval:max_interval],yhat)
      # tck=interpolate.splrep(ltime, yhat, s=1.0)
      # # xnew=np.arange(ltime[0],ltime[-1],100)
      # ynew = interpolate.splev(ltime, tck, der=0)
      # yder = interpolate.splev(ltime, tck, der=1)
      # ax0.plot(xnew,ynew, '-',c=colours[filecounter])
      # 
      # tck = interpolate.splrep(log_t, log_MSD, s=2.0)
      # xnew=np.arange(4,11,0.001)
      # ynew = interpolate.splev(xnew, tck, der=0)
      # # yder = interpolate.splev(xnew, tck, der=1)
      # # print xnew,ynew
      # ax1.plot(log_t,[np.log(x) for x in ynew], '-',c=colours[filecounter])
      # ax1.plot(log_t,[np.log(x) for x in yder], '-',c=colours[filecounter])
      # yder = interpolate.splev(xnew, tck, der=1)
  
  #ax0.set_yscale('log')


  # #plot cell tracks
  # #invert position matrices for plotting (each row is 1 cell instead of 1 timepoint)
  # xx=zip(*xpos)
  # yy=zip(*ypos)
  # toplot=range(0,len(xx))
  # random.shuffle(toplot)
  # #for el in toplot[:tracknr]:
  # #  print el
  # 
  # count=0
  # count2=0
  # for xrow,yrow in zip(xx,yy):
  #   if count in toplot[:tracknr]:
  #     ax1.plot(xrow,yrow, lw=0.5 ,c=colours[filecounter])
  #     count2+=1
  #   count+=1

  filecounter+=1
  if filecounter>=len(colours): 
      filecounter %= len(colours)
      print "Watning, more files than colours"

#No interpolation - just a 45 degree line to see what grows faster than 1
# ax1.plot( [log_t[0], log_t[-1] ] , [log_MSD[0]-log_t[0]+log_t[0], log_MSD[0]-log_t[0]+log_t[-1] ] )
# ax1.plot( [log_t[0], log_t[-1] ] , [log_MSD[0]-log_t[0]+log_t[0], log_MSD[0]-log_t[0]+2.*log_t[-1] ] )
ax2.plot( [log_t[0], log_t[-1] ] , [1,1 ], 'k--' )
# ax1.plot( [log_t[0], log_t[-1] ] , [2,2 ] )

ax0.set_xlabel('time interval (MCS)')
ax0.set_ylabel('MSD (pix^2)')
ax0.set_title('Mean squared displacement')
ax0.set_aspect('equal')
ax1.set_xlabel('log(time interval)')
ax1.set_ylabel('log-log transformed (MSD/$\Delta t$)')
ax2.set_ylabel('Derivative of log-log transformed (MSD/$\Delta t$)')
ax1.set_title('Diffusive exponent')
ax1.legend()
ax1.set_aspect('equal')
ax1.set_xticks([np.log(100), np.log(1000), np.log(10000)])
ax1.set_xticklabels(['100', '1000', '10000'])
ax1.set_ylim([0.5,None])
ax2.set_ylim([0.5,None])
#ax1.set_yticks([0,1,2])
#ax1.set_yticklabels(['0', '1', '2'])

ax2.get_shared_x_axes().join(ax1, ax2)
ax2.set_xticklabels([])
ax2.get_shared_y_axes().join(ax1, ax2)
ax2.set_yticklabels([])

fig.savefig(figname, bbox_inches='tight')
plt.show()


