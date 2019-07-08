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
### --- BLOB PLOT --- ###
###                   ###
#########################
      
#########################
###                   ###
### ---   BEGIN   --- ###
###                   ###
#########################


filename=""

maxTime=-1.
if len(sys.argv)>1:
  pos=1
  while pos < len(sys.argv) :

   if sys.argv[pos]=='-filename':
      #print "Changing filename to",
      pos+=1
      filename=sys.argv[pos]
      print filename
   else:
      print "This is the program 'plot_field.py'"
      print "Here are options, filename is necessary:"
      print "-filename [filename]"
      sys.exit(1)
   pos +=1

#open filename once to see what is the first time step... for now it is typically zero but in the future?
#with open(filename,"r") as fin:
  #for line in fin:
    #line=line.split()
    #inittime=int(line[0])
    #break
#ltime.append(inittime)

with open(filename,"r") as fin:
  line=fin.readline()
  line=line.split() #returns split line
  #print line
  time=int(line[1])
  sizex=int(line[3]) -2 #because boundaries
  sizey=int(line[4]) -2
  cpm=[[0 for _ in range(sizey)] for _ in range(sizex)]
  while True: 
    line = fin.readline()
    if line is None: break
    if 'CA' in line: 
      break
  #Now we can read in the CA data
  i=0
  j=0
  while True: 
    line = fin.readline()
    if "IntPlane" in line: break
    line=[int(sigma) for sigma in line.split()]
    for sigma in line:
      if j==sizey:
        j=0
        i+=1
      cpm[i][j]=sigma
      j+=1
  print line
  print "Done reading"

#print cpm  
cpm = zip(*cpm)
cpm = [ bla for bla in cpm[::-1]]
fig, ax = plt.subplots(1, 1)
c = ax.pcolor(cpm)
ax.set_aspect('equal')
plt.show()
sys.exit(1)
'''
    
    line=line.split()
    time=int(line[0])
    tau=int(line[1])
    birthdate=int(line[2])
    #key=line[3]  #don't care about this right now
    #lock=line[4]
    contacts=line[5:-2]
    #print contacts
    fract_metab=float(line[-1])
    
    if time != ltime[-1]: 
      #calculate statistics
      lJ_av_matrix = np.divide(lJmatrix.astype(float) , lcounter.astype(float), out=np.zeros_like(lJmatrix.astype(float)), where=lcounter!=0)
      
      #lJmatrix.astype(float)/lcounter.astype(float)
      l_avrgdata_time.append(lJ_av_matrix)
      l_avrgdata_time2.append(lJmatrix2)
      #append pop sizes
      lpop_time.append(lpop)
      #append the new time
      ltime.append(time)
      lfract_metab.append( lfract_metab_thistime )
      # zero the counters
      lcounter=np.zeros((3,3),dtype=int)
      lJmatrix=np.zeros((3,3),dtype=int)
      lJmatrix2=[[[] for j in xrange(3)] for i in xrange(3)]
      lpop=[0 for _ in lpop]
      lfract_metab_thistime=[]
    
    lfract_metab_thistime.append(fract_metab)
    
    lpop[tau]+=1 
    for ctau, cJ in [ contacts[pos:pos + 2] for pos in xrange(0, len(contacts), 2) ]:
      ctau=int(ctau)
      lcounter[tau,ctau]+=1
      lJmatrix[tau,ctau]+=int(cJ)
      lJmatrix2[tau][ctau].append(int(cJ))
      if int(cJ)<=0: 
        print "Got cJ<=0 in line:", line
    # Now I want to calculate gamma on the fly, and then average, 
    # but it might not be possible if: not enough data e.g. prey touching only preys
    # also, heterogeneous neighbourhood: gamma is a "large" scale property
    # G10 = J10 - J11/2
    # G20 = J20 - J22/2
    # G12 = J12 - (J11+J22)/2
    
    # we cannot calculate gamma with medium if there's no medium contact
    #if int(contacts[0])==0:
      #for ctau, cJ in [ contacts[pos:pos + 2] for pos in xrange(0, len(contacts), 2) ]:
        #if ctau==0: J0tau=int(cJ)
        #else: ...
          
      #if tau==1:
        #bla
      #elif tau==2
      #else:
        #error
      
      #for ctau, cJ in [ contacts[pos:pos + 2] for pos in xrange(0, len(contacts), 2) ]:
        
      # G10 = J10 - J11/2
      # G20 = J20 - J22/2
      # G12 = J12 - (J11+J22)/2
      
      
    #ldata.append(line)

lpop_time=zip(*lpop_time)
ltime=ltime[:-1]  #removing last data point because it is not appended in the loop above

#print ltime

#print len(ltime)
#print len(lpop_time[0])
#print len(l_avrgdata_time)

#print [ J_av[1,2] for J_av in l_avrgdata_time ]

fig = plt.figure()
ax1 = fig.add_subplot(311)
ax1.plot( ltime,lpop_time[1], label="#prey" )
try: 
  ax1.plot( ltime,lpop_time[2], label="#pred" )
except:
  pass
ax1.legend()

ax2 = fig.add_subplot(312)
#THIS BELOW WORKS FINE, BUT NOW WE USE GAMMAs, NOTICE THAT I AM USIGN AVERAGE J VALUES
a=ax2.plot( ltime, [ J_av[1,0] for J_av in l_avrgdata_time ], label="J10" )
b=ax2.plot( ltime, [ J_av[1,1] for J_av in l_avrgdata_time ], label="J11" )
c=ax2.plot( ltime, [ J_av[1,2] for J_av in l_avrgdata_time ], label="J12" )

d=ax2.plot( ltime, [ J_av[2,0] for J_av in l_avrgdata_time ], label="J20" )
#e=ax2.plot( ltime, [ J_av[2,1] for J_av in l_avrgdata_time ], label="J21" )  #overlaps perfectly with [1,2], as it should
f=ax2.plot( ltime, [ J_av[2,2] for J_av in l_avrgdata_time ], label="J22" )


blobplot_every= 20
blob_plot(ax2,[x[1][0] for x in l_avrgdata_time2][::blobplot_every] ,ltime[::blobplot_every], (a[0].get_color(),0.5),'int', extra='median')
blob_plot(ax2,[x[1][1] for x in l_avrgdata_time2][blobplot_every/2::blobplot_every] ,ltime[blobplot_every/2::blobplot_every], (b[0].get_color(),0.5),'int', extra='median')
#blob_plot(ax2,[x[1][2] for x in l_avrgdata_time2][::blobplot_every] ,ltime[::blobplot_every], (c[0].get_color(),0.5),'int',extra='mean')
#blob_plot(ax2,[x[2][0] for x in l_avrgdata_time2][::blobplot_every] ,ltime[::blobplot_every], (d[0].get_color(),0.5),'int')
#blob_plot(ax2,[x[2][1] for x in l_avrgdata_time2][::blobplot_every] ,ltime[::blobplot_every], (e[0].get_color(),0.5),'int') #overlaps perfectly with [1,2], as it should
#blob_plot(ax2,[x[2][2] for x in l_avrgdata_time2][::blobplot_every] ,ltime[::blobplot_every], (f[0].get_color(),0.5),'int')
#print [x[2][2] for x in l_avrgdata_time2]

#blob_plot(ax2,[x[2][2] for x in l_avrgdata_time2][::10] ,ltime[::10], ('yellow',0.8),'int')

#plt.show()
#print ltime
#sys.exit(1)

#PLOT GAMMAS  
# ~G10~ = <J10> - <J11>/2
# ~G20~ = J20 - J22/2
# ~G12~ = J12 - (J11+J22)/2
          
#ax2.plot( ltime, [ J_av[1,0]-J_av[1,1]/2. for J_av in l_avrgdata_time ], label="G10" )
#ax2.plot( ltime, [ J_av[2,0]-J_av[2,2]/2. for J_av in l_avrgdata_time ], label="G20" )
#ax2.plot( ltime, [ J_av[1,2]- (J_av[1,1]+J_av[1,1])/2. for J_av in l_avrgdata_time ], label="G12" )


ax2.legend()

ax3 = fig.add_subplot(313)
a = ax3.plot(ltime,[np.mean(x) for x in lfract_metab],label="avrg fr metab")
blob_plot(ax3,lfract_metab[::blobplot_every] ,ltime[::blobplot_every], (a[0].get_color(),0.5),'float')
ax3.legend()

title=filename.split('/')[-3]+'_'+filename.split('/')[-2]

fig.suptitle(title)

# Option 2
# TkAgg backend
manager = plt.get_current_fig_manager()
manager.resize(*manager.window.maxsize())
#plt.savefig(title+'.pdf')

plt.show()
sys.exit(1)
'''
