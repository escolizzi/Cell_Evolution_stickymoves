#!/usr/bin/python2.7

'''
Makes a set of blob plots with the last so many time steps for a series of files (ordered by season, maybe?)
Usage: 
./plot_steadystate..blabla.py <output file name> <#timesteps, from the last> file1 file2 file3 ...
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


#########################
###                   ###
### ---   BEGIN   --- ###
###                   ###
#########################

colours=["firebrick","royalblue", "darkgoldenrod", "green", "salmon", "lightskyblue","orchid"]
filename=""
fig, ax = plt.subplots()
if len(sys.argv) < 3:
  print "This is the program 'plot_steadysate...blabla.py'"
  print "Usage:"
  print "./plot_steadystate..blabla.py <output file name> <#steps> <season dur1>,<season dur2>,<season dur3> <file1> <file2> <file3>  ..."
  print "./plot_steadystate..blabla.py <output file name> <#steps> -e <file1> <file2> <file3>  ..."
  print "  # steps is how many steps back you want to see, like... 5"
  print "  Comma-separate season duration, and space separate filenames - the two have to match"
  print "  Alternatively: -e option, attempts extracting season duration from file name (but may fail)"
  sys.exit(1)
else:
  figname=sys.argv[1]
  timesteps = int(sys.argv[2])
  if sys.argv[3]=='-e':
    llabel=[]
    for bla in sys.argv[4:]:
      print "filename: ", bla
      start = bla.find("_s")
      lss = len("_s")
      end = bla.find("k_")
      label = bla[ start+lss: end ]
      # ~ print "label: ", label
      llabel.append(label+'000')
  else: 
    llabel=sys.argv[3].split(',')
    if llabel[-1]=='': llabel=llabel[:-1]
  lfilename=sys.argv[4:]

if len(llabel) != len(lfilename): 
  print "Error, different number of labels and filenames?"
  print llabel
  print lfilename
  sys.exit(1)

#the figure
fig, ax = plt.subplots()
llJmed=[]
llJcel=[]
lmeanJc=[]
lmeanJm=[]

#especially if we got option -e, chances are that filenames and labels are sorted aplhanumerically, instead of numerically
# so we sort
llabel.sort(key = lambda x: int(x))
try:
  lfilename.sort( key = lambda x: int(   x[ x.find("_s") +2 : x.find("k_") ]  ) )
except:
  try:
    lfilename.sort( key = lambda x: int(   x[ x.replace('_s', 'XXX', 1).find('_s') +2 : x.find("k_") ]  ) ) #looks for the second
  except:
    print "Error in sorting lfilename"
    sys.exit(1)

print llabel
print lfilename

#some smart looping to get season and filename
for label,filename in zip(llabel,lfilename):
    #we get season from label
    season = int(label)
    #last saved time step is
    output = subprocess.Popen(['tail', '-1', filename], stdout=subprocess.PIPE).communicate()[0]
    lasttime = int(output.split(' ')[0])
    # by the way, data is saved every 10000
    # so if the season is 10000 or less we take the last 10 time steps and that's it
    
    # if season<=10000:
    #     ltime = [lasttime - i*10000 for i in range(10)]
    # 
    # #else, we should take a time point for each of the last 10 seasons
    # # which one? ... the middle one!
    # # BUT WHAT IF SEASON / 2 is not divisible by 10000? -> should get closest point to that...
    # # so instead of middle one we do trimmed_lasttime - 100000
    # else:
    # I really don't know what the hell I was doing here
    # ~ trimmed_lasttime = lasttime - lasttime%season 
    # ~ timepoint = trimmed_lasttime - 10000
    # ~ ltime = [lasttime - i*season for i in range(10)]
    
    ltime =[lasttime]
    
    # print ltime
    # Now we go get the date
    lJmed=[]
    lJcel=[]
    with open(filename,"r") as fin:
        print "Opening", filename
        # ~ print "ltime", ltime
        # for line in fin[::-1]:
        for line in reversed(fin.readlines()):
            # print line
            line=line.split(' ')
            time = int(line[0])
            # print time
            if time not in ltime:
              # ~ print time,"not in ltime: ", ltime
              if len(ltime)<=timesteps:
                # ~ print "Hello"
                ltime.append(time)
              else:
                # ~ print "we break because", len(ltime),">",timesteps
                break
                # print "Hello"
                # print line
                #save lJo, lJ1
            jvalwithnei_startpoint=31 # in normal files it is 31
            if line[jvalwithnei_startpoint]=='0': 
                    lJmed.append(int(line[jvalwithnei_startpoint+1])) #the line[30] is J val with med
                    lJcel.extend([ int(x) for x in line[jvalwithnei_startpoint+3::2] ] ) 
            else:
                    lJcel.extend([ int(x) for x in line[jvalwithnei_startpoint+1::2] ] ) 
    
    print ltime
    llJcel.append(lJcel)
    llJmed.append(lJmed)
    # At this point we can start filling up our figure
    # meanJc=blob_plot(ax,[lJcel],[season], ('royalblue',1.),'int' ,extra='mean')
    # meanJm=blob_plot(ax,[lJmed],[season], ('firebrick',1.),'int' ,extra='mean')
    
    # lmeanJc.extend(meanJc)
    # lmeanJm.extend(meanJm)
    # ax.boxplot(lJmed, positions=[math.log(season)])
    # ax.boxplot(lJcel, positions=[math.log(season)])

#sets flier properties    
#flierprops = dict(marker='.', markerfacecolor='k', markersize=1,linestyle='none') <- this works perfectly, I just don't want to show fliers any more, because nothing Changes

position_offset = 1000 # so that Jc and Jm don't overlap
labels=[int(label[:-3]) for label in llabel ]
print "labels:", labels

ax.boxplot(llJcel, positions=[int(x)-2.5*position_offset if i==llabel.index(x) else int(x)+1.25*position_offset  for i,x in enumerate(llabel)], widths=(2000), showfliers=False, medianprops=dict(color='royalblue') )
lmeanJc= [ np.mean(lJcel) for lJcel in llJcel]
# ax.plot(llabel,lmeanJc,color='royalblue')
ax.plot(llabel,lmeanJc,linestyle = '-',marker='o', markerfacecolor='None', markeredgecolor='royalblue', markersize = 1)

ax.boxplot(llJmed, positions=[int(x)-1.25*position_offset if i==llabel.index(x) else int(x)+2.5*position_offset  for i,x in enumerate(llabel)], widths=(2000), showfliers=False, medianprops=dict(color='darkgoldenrod') )
lmeanJm= [ np.mean(lJm) for lJm in llJmed]
# ax.plot(llabel,lmeanJm,color='darkgoldenrod')
ax.plot(llabel,lmeanJm,linestyle = '-',marker='o', markerfacecolor='None', markeredgecolor='darkgoldenrod', markersize = 1)

# ax.plot(llabel , [ jm-jc/2. for jc,jm in zip(lmeanJc,lmeanJm)], color='seagreen')
ax.plot(llabel , [ jm-jc/2. for jc,jm in zip(lmeanJc,lmeanJm)], linestyle = '-',marker='o', markerfacecolor='None', markeredgecolor='seagreen', markersize = 1)

ax.plot(llabel,[28. for _ in llabel], lw = 0.5, color='royalblue')
ax.plot(llabel,[14. for _ in llabel], lw = 0.5, color='darkgoldenrod')
ax.plot(llabel,[0. for _ in llabel], lw = 0.5, color='seagreen')

ax.set_xticklabels(labels)
ax.set_xticks([int(label) for label in llabel ])
minor_ticks = np.arange(-10, 46, 5)

ax.set_yticks(minor_ticks,minor=True)

ax.set_xlabel('$\\tau_s$ [$\cdot 10^3$ mcs]', fontsize=12)

ax.set_xlim([0,510000])
ax.set_ylim([-10,10])

# plt.show()    
plt.savefig(figname)

sys.exit(1)


