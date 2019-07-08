#!/usr/bin/python2.7

'''
minimal plotting tool for pCellEvol pred prey
the file is organised like this:
for all pop, every so many time steps (typically 1000), each individual is dumped, information is:
Time tau birthday key lock tau_contact J_contact tau_contact J_contact tau_contact J_contact ...
'''

import sys,math,os,subprocess
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
### ---   BEGIN   --- ###
###                   ###
#########################


#lx=ldata[-2]
#ly=ldata[-1]

lx = np.random.uniform(low=-1., high=1., size=(100000,))
ly = np.random.uniform(low=-1., high=1., size=(100000,))

ltheta=[ math.atan2(float(y),float(x))/math.pi*180. for x,y in zip(lx,ly) ]

#print ltime

#print len(ltime)
#print len(lpop_time[0])
#print len(l_avrgdata_time)

#print [ J_av[1,2] for J_av in l_avrgdata_time ]

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.hist(ltheta, normed=True, bins=180)


plt.show()
sys.exit(1)
