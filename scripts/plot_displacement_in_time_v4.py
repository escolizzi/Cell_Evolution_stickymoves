#!/usr/bin/python2.7

'''
v4: the whole point of this version is to strengthen name dependence, 
so that arbitrary groups can be made, each with their name tag (optional)


Plot position in time, takes all cells in a time step and calculates average position
Plots the average of the average position, divided by filename
Also depending on file name plots single cellfile name dependence is VERY brittle.
Notice that 
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

def AverageFrom_lg_lists(lg):
  av_lg=[]
  min_lg=[]
  max_lg=[]
  stdev_lg=[]
  
  len_lg = min( [ len(x) for x in lg ] )
  for i in range(len_lg):
    av_lg.append( np.mean( [x[i] for x in lg] ) )
    stdev_lg.append( np.std( [x[i] for x in lg] ) )
    min_lg.append( min([x[i] for x in lg]) )
    max_lg.append( max([x[i] for x in lg]) )
  ltime = [ 10*x for x in range(len_lg) ]
  
  return ltime,av_lg,stdev_lg,min_lg,max_lg


#########################
###                   ###
### ---   BEGIN   --- ###
###                   ###
#########################


# ~ colours=[(0.7,0.13,0.),"royalblue", "darkgoldenrod", (0.,0.5,0.2), "salmon", "lightskyblue","orchid"]

fig, (ax0, ax1) = plt.subplots(ncols=2,sharey=True)
plt.subplots_adjust(wspace=0.05)

if len(sys.argv) <3:
  print "This is the program 'plot_diplacement_in_time_v2'"
  print "Usage: ./plot_displ.....py <figure name> <peak_row> <peak_col> <filename1> ..."
  print "For groups of files use the letter g for separating groups, and space between groups:" 
  print "./plot_displ.....py <figure name> <peak_row> <peak_col> g <file1> <file2> g <file3> <file4> ..."
  sys.exit(1)
else:
  figname=sys.argv[1]
  peakx=int(sys.argv[2])
  peaky=int(sys.argv[3])

print "Calculating distance from peak at (row,col) = ", peakx,peaky

lgroups_filename=[]
if( sum([1 if x=='g' else 0 for x in sys.argv[4:]]) == 0 ):
  lgroups_filename=[ [x] for x in sys.argv[4:] ] # each file is its own group
else:
  for x in sys.argv[4:]:
    if x=='g':
      lgroups_filename.append([])
    else:
      #check if file exists
      if(os.path.exists(x)):
        lgroups_filename[-1].append(x)
      else:
        print "File not found: ", x
        sys.exit(1)
print lgroups_filename
len_lgroups_filename = len(lgroups_filename)
print "# groups", len_lgroups_filename

# ~ choose colours now that we know how many groups we have
colours=cm.gist_earth( np.linspace(0.1,0.9,len_lgroups_filename ) )

lldist=[[] for groups in lgroups_filename] # ldist for all files, organised to match file names
# ~ sys.exit(1)
for group_number,group in enumerate(lgroups_filename):

  for bla,filename in enumerate(group):

    print "reading file ",filename
    
    av_xpos=[];av_ypos=[];min_xpos=[];min_ypos=[];
    ltime=[]
    lxpos=[]
    lypos=[]
    c1_pos=[]
    
    ldist=[]
    av_dist=[]
    min_dist=[]
    
    # read time in the first line
    output = subprocess.Popen(['head', '-1', filename], stdout=subprocess.PIPE).communicate()[0]
    time = int(output.split(' ')[0])
    ltime = [time]        
    
    with open(filename,"r") as fin:
			for line in fin:
				line=line.split()
		  
				time = int(line[0])
				xpos = float(line[3])
				ypos = float(line[4])
				dist = np.hypot( xpos-peakx,ypos-peaky )
				#if time of this line different from time of previous line
				# we do aggreagate statistics
				if time != ltime[-1]:
					#xpos and ypos are not really used ..
					#av_xpos.append(  np.mean( lxpos )  ); min_xpos.append( min(lxpos) ); lxpos=[]
					#av_ypos.append(  np.mean( lypos )  ); min_ypos.append( min(lypos) ); lypos=[] 
					
					av_dist.append( np.mean(ldist) ) #average position of every cell at this time
					min_dist.append( min(ldist)) #cell closest to peak at this point
					ldist=[] # resets dist to zero for next time step
					
					ltime.append(time)
				
				ldist.append(dist) 
				lxpos.append(xpos);lypos.append(ypos)
    
    # ~ lldist[group_number].append(av_dist)
    ax0.plot( ltime[:-1],av_dist, c= colours[group_number], lw=0.71  )
    
	# At this point we have ldist for this file	
  #ax0.plot( ltime[:-1],av_dist, c= colours[index], lw=0.5  )
  #lldist[group_number].append(ldist)
#for groups in lldist:
#  print "Hello"

ax0.set_title('Center of mass', fontsize=16)
ax0.legend(title='$\gamma_{c,m}$')
start, end = ax0.get_xlim() # this and next changes axis location

# ~ ax0.xaxis.set_ticks(np.arange(0, 200000, 25000))

# ~ start = 'mup'
# ~ end ='_'
# ~ s = lgroups_filename[-1][-1]
# ~ try:
  # ~ startpoint=s.find(start)+len(start)
  # ~ endpoint= s.find(start) + s[s.find(start):].find(end)

  # ~ print "starting from: ", s
  # ~ print 
  # ~ print "Got: ",s[startpoint: endpoint]
  # ~ print
  # ~ mup = float( s[startpoint: endpoint] )

  # ~ if mup<=1.:
    # ~ xmax=200000 
  # ~ elif mup<4.:
    # ~ xmax=100000
  # ~ else:
    # ~ xmax=50000
# ~ except:
  # ~ xmax=100000
  # ~ pass
  
xmax = 200000
ax0.set_xlim([0,xmax])
ax0.set_ylim([0,500])

ax0.set_xlabel('time (MCS)', fontsize=14)
ax0.set_ylabel('distance from gradient peak', fontsize=14)

ax1.set_title('Cell closest to\ngradient peak', fontsize=16)
ax1.set_xlim([0,xmax])
ax1.set_ylim([0,None])
ax1.set_xlabel('time (MCS)', fontsize=14)
fig.savefig(figname, bbox_inches='tight')
plt.show()
sys.exit(1)
  
