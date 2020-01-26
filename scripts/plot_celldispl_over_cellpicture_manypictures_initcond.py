#!/usr/bin/python2.7

'''
Plot the displacement of some cells over the time period of two pictures.

'''

import sys,math,os,subprocess,random
from PIL import Image
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

def GetTimeFromPath(bla):
    bla=bla.split('.')[-2]
    s_init_time = ''
    for x in bla[::-1]:
        try: 
            int(x)
        except:
            break
        s_init_time += x
    # print s_init_time[::-1]    
    return int(s_init_time[::-1])

def CheckIfNeiIsYellow(x,y,image):
    for i in range(-1,2):
        for j in range(-1,2):
            # print "hello"
            try: 
                r,g,b,alpha = image.getpixel((x+i,y+j))
                # print "x,y=",x,y,"rgb=",r,g,b
                if r==254 and g>200 and b==0:
                    # print "got yellow nei"
                    return True
            except:
                pass    
    return False

def CheckIfNeiIsNotCell(x,y,image):
    for i in range(-1,2):
        for j in range(-1,2):
            # print "hello"
            try: 
                r,g,b,alpha = image.getpixel((x+i,y+j))
                # print "x,y=",x,y,"rgb=",r,g,b
                if not(r==200 and g==30 and b==30) and not(r==0 and g==0 and b==0) and not(r==0 and g==51 and b==102):
                    # print "got yellow nei"
                    return True
            except:
                pass    
    return False
  
#########################
###                   ###
### ---   BEGIN   --- ###
###                   ###
#########################

colours=["firebrick","royalblue", "darkgoldenrod", "green", "salmon", "lightskyblue","orchid"]
filename=""
fig, ax = plt.subplots()
if len(sys.argv) <3:
  print "This is the program 'plot_CMdisplacement_manypictures.py'"
  print "Usage: ./plot_celldispl_ blabla .py <output name> <peak row> <peak col> <path/to/data_filename> <path/to/png1> <path/to/png2> <path/to/png3> ..."
  sys.exit(1)
else:
  figname=sys.argv[1]
  peak_row=int(sys.argv[2])
  peak_col=int(sys.argv[3])
  filename=sys.argv[4]
  lpath_to_png=sys.argv[5:]
  
# Get initial and final time, to get trajectories from datafiles
picture_times=[GetTimeFromPath(path_to_png) for path_to_png in lpath_to_png]
init_time = min(picture_times)
end_time = max(picture_times)

#I should really finish drawing the contours of the images:
limage = [Image.open(path_to_png).convert('RGBA') for path_to_png in lpath_to_png]
width, height = limage[0].size
for x in range(width):
    for y in range(height):
        for image in limage:
            r,g,b,alpha = image.getpixel((x,y))
            # print r,g,b
            if (r==200 and g==30 and b==30) or (r==0 and g==51 and b==102):
                # print "a red pixel at x,y",x,y
                if CheckIfNeiIsNotCell(x,y,image): 
                    image.putpixel((x,y),(0,0,0))
        

print "Done with coloring borders"
# Extract all that is black
print "w,h = ", width,height
# sys.exit(1)
ldata = [np.array(image) for image in limage]

#blacks = np.zeros(shape=(width,height)) #it is 1 where black
# blacks = np.zeros(shape=(width,height,3)) #it is 1 where black
data_colours = np.zeros((height,width), dtype=(float,4))
for i in xrange(height):
    for j in xrange(width):
        for data in ldata:
            if data[i][j][0] == 0 and data[i][j][1]== 0 and data[i][j][2]==0: 
                # blacks[i][j]=(0.2,0.2,0.2,1.)
                data_colours[i][j]=(0.,0.,0.,1.)
            elif data[i][j][0] == 200 and data[i][j][1]== 30 and data[i][j][2]==30:
                data_colours[i][j]=(200/255.,30/255.,30/255.,1.)
            elif  data[i][j][0] == 0 and data[i][j][1]== 51 and data[i][j][2]==102:
                data_colours[i][j]=(10/255.,200/255.,250/255.,1.)
print "Done extracting black" 

# dpi=40
# figsize = width / float(dpi), height / float(dpi)
# fig, ax = plt.subplots(figsize=figsize)
# ax.imshow(blacks,extent=[0, 600, 400, 0],zorder=3,interpolation='none')
# plt.gca().invert_yaxis()
# plt.savefig(figname)
# sys.exit(1)


# lgrimage = [ image.convert('L') for image in limage ]
# # grimage2.show()
# # sys.exit(1)
# print "Greys before lightneing"
# for bla in sorted( lgrimage[0].getcolors()):
#      print bla
# # sys.exit(1)     
# 
# # make every gray slightly lighter, but not white and black
# for x in range(width):
#     for y in range(height):
#         # grey,alpha = image1.getpixel((x,y))
#         # if grey<250 and grey> 10: 
#         #     # print "old", grey
#         #     # grey = int(  500.* float(grey)/( 250.+float(grey) )  )
#         #     # grey +=20
#         #     grey = grey
#         #     # print "new", grey
#         #     image1.putpixel((x,y),grey)
#         grey = grimage2.getpixel((x,y))
#         # print "grey",grey
#         if grey < 150: 
#             grey += 120
#             grimage2.putpixel((x,y),grey)
#         elif grey<250: 
#             # print "old", grey
#             grey = int(  450.* float(grey)/( 200.+float(grey) )  )
#             # grey +=20
#             #grey = grey
#             # print "new", grey
#             grimage2.putpixel((x,y),grey)
# print "Greys after lightening"
# for bla in sorted( grimage2.getcolors()):
#      print bla
# # grimage2.show()
# # sys.exit(1)
# # 
# grimage2 = grimage2.convert('RGB')
# grimage2_arr = np.asarray(grimage2) # In order to pass it to imshow

# Now open file, from point init_time to point end_time get tracks (as x,y, coord)
ltime=[]
d1={}
with open(filename,"r") as fin:
    for line in fin:
        line=line.split()
        timepoint=int(line[0])
        if timepoint<init_time: 
            continue
        elif timepoint>end_time: 
            break
        ltime.append(timepoint)
        # could make a dict of sigmas that hold two lists with x and y pos
        sigma=int(line[1])
        xpos=float(line[3])
        ypos=float(line[4])
        if(sigma in d1):
            d1[sigma][0].append(xpos)
            d1[sigma][1].append(ypos)
        else:
            d1[sigma]=[[xpos],[ypos]]
            
print "Done reading, now overlaying"

dpi=40
figsize = width / float(dpi), height / float(dpi)
fig, ax = plt.subplots(figsize=figsize)
# img = plt.imread("new.png")
# x = range(300)

# extent = [ horizontal min, horizontal max, vertical_min,vertical_max ]
# ax.imshow(grimage2_arr, extent=[0, 600, 400, 0],zorder=1,interpolation='none')
ax.imshow(data_colours,origin='lower',extent=[0, width/2, 0, height/2],zorder=1,interpolation='none') #extent specifies the coordinates of the boundaries of the picture, 
                                                                                                # shrinking or stretching according to xlim ad ylim
                                                                                                # use width/2 , height/2 because the pngs are always
                                                                                                # twice the number of pixels as the CA

# plt.show()
# for bla in img:
#     print bla
#     sys.exit(1)
#     for x in bla:
#         print x
#         sys.exit(1)

# circles=[]
for x in range(0,600,100):
    circle = plt.Circle((peak_row, peak_col), x, color='k', ls='--', linewidth=1, fill=False) 
    ax.add_artist(circle)
# for circle in circles:
#     ax.add_artist(circle)

colour=0
for i in range(400):
    # ax.plot(d1[1][1], d1[1][0], '-', linewidth=1.5, color=colours[0],zorder=2)  #these are plotted as normalcoordinates, i.e. with origin at bottom left
    # ax.plot(d1[11][1], d1[11][0], '-', linewidth=1.5, color=colours[1],zorder=2)
    # ax.plot(d1[89][1], d1[89][0], '-', linewidth=1.5, color=colours[3],zorder=2)
    
    try:
        ax.plot(d1[i][1], d1[i][0], '-', linewidth=1.5, color=colours[colour],zorder=2)
        colour+=1
        if colour >= len(colours): break
    except: 
        pass
# ax.imshow(blacks,extent=[0, 600, 400, 0],zorder=3)

# ax.set_ylim((300,700)) #if peak at 500, 0
# ax.set_xlim((0,500))
ax.set_ylim((0,height/2)) #full pictures
ax.set_xlim((0,width/2))

# plt.savefig('1_'+figname)
plt.gca().invert_yaxis()
# plt.savefig('2_'+figname)
# SOLVED: print "Problem 1, you get a low resolution, probably interlaced picture as output..."
# SOLVED print "Problem 2, the picture is upside down -> rotate along horizontal axis"
# SOLVED print "Problem 3, it would be nice to get the grays a little lighter"
print "Problem 4, it would be cooool to add on top a small plot that tells you the average amount of gradient present (polar)"

# plt.show()
plt.savefig(figname)
sys.exit(1)



