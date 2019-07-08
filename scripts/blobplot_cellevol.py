#!/usr/bin/python2.7

'''
The idea is to compare the variances in interdiv time of two poulations with the variances of their 
respective ancestors.
This should make me see higher fitness variance in pop but lower in anc when mut are biased
(By looking at the distribution of variances vs distribution of ancestors).
Also, I want to make 2D plots of active and inactive genes, both for ancestors and for populations.
This should help me see that biased mut rate explore their mut nei better.
'''

import sys,math,os,string
from subprocess import Popen, PIPE
#from PIL import Image
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec
import numpy as np

lanc=[]
lancgen=[]
lancT=[];lancQ=[];lancR=[];lancB=[];lancP=[]
lanct=[];lancq=[];lancr=[];lancp=[]

# COLOR MAP #
#colorT='blue'
#colorQ='green'
#colorRr='red'
#colorRp='cyan'

# COLOR MAP #
colorT=(56./255.,146./255.,56./255.)# green (medium blue)
colorQ=(0/255.,109./255.,240/255.)	# blue (lotsa blue)
colorRr=(166/255.,0/255.,0/255.)	# red (no blue)
colorRp=(26/255.,26/255.,30/255.)	# grey - antracitis

color1=(6./255.,6./255.,98/255.)
color2=(128/255.,0/255.,0/255.)
color3=colorT
color4=colorRp

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


#########################
###                   ###
### --- BLOB PLOT --- ###
###                   ###
#########################
def blob_plot(ax,data,list_pos, colour,datatype,extra=False):
    lextra=[]
    for d,p in zip(data,list_pos):
      if not d: 
        
        continue
        
      #this makes the histogram that constitutes the blob plot
      m=min(d) # the largest value, upper edge of the histogram
      M=max(d) # the lowest value, lower...
      if datatype=='int':
        his=np.bincount(d)
        
        his=np.trim_zeros(his)
        his=np.append(his,[0])
        his=np.append([0],his)
        x=np.linspace(m-1,M+1,M-m+3 )
        
        #print "his",his
        #print "x",x
        '''
        if M-m==0:
          #in these cases histograms are awkard or do not display at all, so we add an extra "zero" data point
          nbins=3
          x = np.linspace(m-1,M+1,nbins) # support for histogram: this is where the histogram is placed
          print "M-m=0; x",x
        elif M-m==1:
          nbins=3
          x = np.linspace(m-1,M,nbins) # support for histogram
          print "M-m=1; x",x
        else:
          nbins=M-m+1
          x = np.linspace(m,M,nbins) # support for histogram
          print "M-m>1; x",x
          
        his,bins = np.histogram(d, bins=x)
        
        #plt.plot(np.append(his,0),bins)
        #plt.show()
      '''
      elif datatype=='float':
        nbins=25
        x = np.linspace(m,M,nbins) # support for histogram of floats
        
        #nbins=25  
        #print "min,max,nbins=",m,M,nbins
      
        his,bins = np.histogram(d, bins=nbins)
      '''  
      if m==M or M-m==1 and datatype=='int': 
        print "m==M or M-m==1"
        print "his",his
        print "bins",bins
      '''     
      
      #scale the histogram to a proper size (you may have to fiddle with this), then place it at position 'pos'
      
      
      #scale_blob=0.5
      
      maxwidht=0.15
      max_his=max(his)
      #print max_his
      if max_his>0.:
        scale_blob=maxwidht/max_his
      else: 
        scale_blob=0.
        print "Warning, max_his is 0." 
        
      shift_his_plus_pos =  [ p + h*scale_blob  for h in his]
      shift_his_minus_pos = [ p - h*scale_blob  for h in his]
      # OLD ->   shift_his_plus_pos =  [ p + h*scale_blob/float(len(d))  for h in his]
      #          shift_his_minus_pos = [ p - h*scale_blob/float(len(d))  for h in his]
      
      color_alpha_blobs=colour
      facecolor,alpha=color_alpha_blobs # assign color and transparency
      #this is the matplotlib function that does the trick
      ax.fill_betweenx(x,shift_his_minus_pos, shift_his_plus_pos, linewidth=0., facecolor= facecolor, alpha=alpha)
      #calculates median or mean, if you want
      if extra=='median': lextra.append( np.median(d) )
      elif extra=='mean': lextra.append( np.mean(d) )
    #and plots it
    if extra != False:
      color = 'orangered'
      ax.plot(list_pos,lextra, color=color,linestyle='-',marker='D',markersize=5, lw=1.5)
      
 

    ##############################
#---##########   MAIN   ##########---#
    ##############################

d1={}
ltime=[]

howlong=250000
ldir=[]
lhowlong=[]
#go through files in current dir
data_file="data_cellcount_1.txt"

#### ARGUMENTS ####

if len(sys.argv)==1:
  print "Specify some directories for data"
  print "Also make sure directories exist, and that their path is either in ./ or in ../ or that you are providing a correct absolute path"
elif len(sys.argv)>1:
  #we find path names and howlong
  for x in sys.argv[1:]:
    ldir.append(x)
  print ldir
  
#### CHECK DIRECTORIES EXISTENCE ####
lpath=[]
for d in ldir:
  if "../" in d and os.path.isdir(d):
    #print 'hellooooooooooooooooo'
    path=d
    if path[-1]!='/': path+="/"	# to make a dir name end with '/'
    lpath.append(path)	# it's an absolute path and we append it as is
    continue
  
  #alternatively we look for directories either where the program is called 
  # or one level up
  path="./"+d
  if os.path.isdir(path): 
    
    if path[-1]!='/': path+="/"	# to make a dir name end with '/'
    lpath.append(path)
    continue
  
  path="../"+d
  if os.path.isdir(path): 
    
    if path[-1]!='/': path+="/"
    lpath.append(path)
    continue

print "confirmed paths:"
print lpath
print

#### FETCH FILES ####

lfile=[]
# We go get the files
for path in lpath:
  if os.path.isfile(path+data_file): 
    lfile.append(path+data_file)
  else:
    print "data file not found in ", path


#### GET DATA ####

#lin=[]
#ldata=[]
#lpop=[]
#lJmatrix2 = [[[] for j in xrange(3)] for i in xrange(3)]

#fig, ax = plt.subplots(nrows=1, ncols=5,sharey=True)
#fig = plt.figure(constrained_layout=True)
#fig = plt.figure()
gs = gridspec.GridSpec(2, 5,height_ratios=[3, 1])

ax00 = plt.subplot(gs[0,0])
ax01 = plt.subplot(gs[0,1],sharex=ax00,sharey=ax00)
#ax02 = plt.subplot(gs[0,2],sharey=ax01)   # <- just to show that you can share axis
ax02 = plt.subplot(gs[0,2],sharey=ax00,sharex=ax00)
ax03 = plt.subplot(gs[0,3],sharey=ax00,sharex=ax00)
ax04 = plt.subplot(gs[0,4],sharey=ax00,sharex=ax00)
ax1all = plt.subplot(gs[1,:])

#how many time steps do we go back?
howmanytimesteps=1
counter=0
gotalldata=False
totfiles=len(lfile)
colors = cm.get_cmap('Dark2')    # 11 discrete colors
#print colors(0.5)
alpha=1.
colors=[(colors(x/float(totfiles+1))[:3],alpha) for x in xrange(totfiles)]

for filenumber,(bla,path) in enumerate(zip(lfile,lpath)):
  lrawdata=[] #transient variables
  lrawpop=[]
  lJmatrix2 = [[[] for j in xrange(3)] for i in xrange(3)]
  
  aprocess = Popen(['tail', '-1', bla], stdout=PIPE, stderr=PIPE) #get last line for last time step
  stdout, stderr = aprocess.communicate()
  lasttime=int(stdout.split()[0])
  thistime=int(lasttime)
  
  aprocess = Popen(['tail', '-500000', bla], stdout=PIPE, stderr=PIPE) #get bunch of lines
  stdout, stderr = aprocess.communicate()
  data=stdout.split('\n')[::-1][1:] #reverse it so that we start from end, and remove the last empty line
  counter=0
  for line in data:
    line=line.split()
    currenttime= int(line[0])
    
    if currenttime!=thistime:
      thistime=currenttime
      counter+=1
      #print counter
    if counter==howmanytimesteps:
      gotalldata=True
      break
    
    #ldata.append(line)
    tau=int(line[1])
    contacts=line[5:-2]
    for ctau, cJ in [ contacts[pos:pos + 2] for pos in xrange(0, len(contacts), 2) ]:
      ctau=int(ctau)
      #lcounter[tau,ctau]+=1
      #lJmatrix[tau,ctau]+=int(cJ)
      lJmatrix2[tau][ctau].append(int(cJ))
 
  if not gotalldata:
    print
    print "Warning, did not get all the data"
    print
  
  
  if lJmatrix2[1][0]:
    blob_plot(ax00,[lJmatrix2[1][0]] ,[filenumber], colors[filenumber],'int',extra='mean');
  if lJmatrix2[1][1]:  
    blob_plot(ax01,[lJmatrix2[1][1]] ,[filenumber], colors[filenumber],'int',extra='mean');  
  if lJmatrix2[1][2]:
    blob_plot(ax02,[lJmatrix2[1][2]] ,[filenumber], colors[filenumber],'int',extra='mean');  
  if lJmatrix2[2][0]:
    blob_plot(ax03,[lJmatrix2[2][0]] ,[filenumber], colors[filenumber],'int',extra='mean');  
  if lJmatrix2[2][2]:
    blob_plot(ax04,[lJmatrix2[2][2]] ,[filenumber], colors[filenumber],'int',extra='mean');  
  


ltext=[ str(i)+': '+bla.split('/')[-3]+'_'+bla.split('/')[-2]+'\n' for i,bla in enumerate(lfile) ]
#print ltext
ax1all.text(0.1, 0., ''.join(ltext) )

#title=filename.split('/')[-3]+'_'+filename.split('/')[-2]

ax00.set_title('J10')
ax01.set_title('J11')
ax02.set_title('J12')
ax03.set_title('J20')
ax04.set_title('J22')

# TkAgg backend
manager = plt.get_current_fig_manager()
manager.resize(*manager.window.maxsize())
#plt.tight_layout()
plt.show()
sys.exit(0)

