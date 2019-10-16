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

#########################
###                   ###
### ---   BEGIN   --- ###
###                   ###
#########################

colours=["firebrick","royalblue", "darkgoldenrod", "green", "salmon", "lightskyblue","orchid"]
filename=""
fig, ax = plt.subplots()
if len(sys.argv) <3:
  print "This is the program 'plot_CMdisplacement.py'"
  print "Usage: ./plot_celldispl_ blabla .py <output name> <path/to/data_filename> <path/to/png1> <path/to/png2> "
  sys.exit(1)
else:
  figname=sys.argv[1]
  filename=sys.argv[2]
  path_to_png1=sys.argv[3]
  path_to_png2=sys.argv[4]
  
# Get initial and final time, to get trajectories from datafiles
init_time = GetTimeFromPath(path_to_png1)
end_time = GetTimeFromPath(path_to_png2)

#I should really finish drawing the contours of the images:
image1 = Image.open(path_to_png1).convert('RGBA')
image2 = Image.open(path_to_png2).convert('RGBA')
#red is (200, 30, 30), if red is next to yellow (254, 214 and higher, 0) then it should be turned to black.
width, height = image1.size
for x in range(width):
    for y in range(height):
        r,g,b,alpha = image1.getpixel((x,y))
        # print r,g,b
        if r==200 and g==30 and b==30:
            # print "a red pixel at x,y",x,y
            if CheckIfNeiIsYellow(x,y,image1): 
                image1.putpixel((x,y),(0,0,0))
        #same for image2
        r,g,b,alpha = image2.getpixel((x,y))
        if r==200 and g==30 and b==30:
            if CheckIfNeiIsYellow(x,y,image2): 
                image2.putpixel((x,y),(0,0,0))
        
print "Done with coloring borders"
# image2.show()
# sys.exit(1)
# Extract all that is black
# print image1.size
# print "w,h = ", width,height
data1 = np.array(image1)
data2 = np.array(image2)
#blacks = np.zeros(shape=(width,height)) #it is 1 where black
# blacks = np.zeros(shape=(width,height,3)) #it is 1 where black
blacks = np.zeros((height,width), dtype=(float,4))
for i in xrange(height):
    for j in xrange(width):
        if data1[i][j][0] == 0: 
            blacks[i][j]=(0.2,0.2,0.2,1.)
        if data2[i][j][0] == 0: 
            blacks[i][j]= (0.2,0.2,0.2,1.)
print "Done extracting black" 

# Load the two images and convert them to gray
# image1 = Image.open(path_to_png1).convert('LA')
# image2 = Image.open(path_to_png2).convert('LA')
grimage1 = image1.convert('L')
grimage2 = image2.convert('L')
# grimage2.show()
# sys.exit(1)
print "Greys before lightneing"
for bla in sorted( grimage2.getcolors()):
     print bla
# sys.exit(1)     
# make every gray slightly lighter, but not white and black
for x in range(width):
    for y in range(height):
        # grey,alpha = image1.getpixel((x,y))
        # if grey<250 and grey> 10: 
        #     # print "old", grey
        #     # grey = int(  500.* float(grey)/( 250.+float(grey) )  )
        #     # grey +=20
        #     grey = grey
        #     # print "new", grey
        #     image1.putpixel((x,y),grey)
        grey = grimage2.getpixel((x,y))
        # print "grey",grey
        if grey < 150: 
            grey += 120
            grimage2.putpixel((x,y),grey)
        elif grey<250: 
            # print "old", grey
            grey = int(  450.* float(grey)/( 200.+float(grey) )  )
            # grey +=20
            #grey = grey
            # print "new", grey
            grimage2.putpixel((x,y),grey)
print "Greys after lightening"
for bla in sorted( grimage2.getcolors()):
     print bla
# grimage2.show()
# sys.exit(1)
# 
#         #do it also for the other image
# 
# # Blend the two images
# new_img = Image.blend(image1, image2, 0.75)
# # new_img.show()
# print "images have been greyed and blended"
# # Print colors
# # print "Colors: "
# # for bla in sorted( new_img.getcolors()):
# #     print bla
# 
# 
# new_img.save("new.png","PNG")

# grimage2.save("new.png","PNG")

grimage2 = grimage2.convert('RGB')
grimage2_arr = np.asarray(grimage2) # In order to pass it to imshow

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

dpi=80
figsize = width / float(dpi), height / float(dpi)
fig, ax = plt.subplots(figsize=figsize)
# img = plt.imread("new.png")
# x = range(300)

ax.imshow(grimage2_arr, extent=[0, 600, 400, 0],zorder=1,interpolation='none')


# plt.show()
# for bla in img:
#     print bla
#     sys.exit(1)
#     for x in bla:
#         print x
#         sys.exit(1)

ax.plot(d1[1][1], d1[1][0], '-', linewidth=1.5, color=colours[0],zorder=2)
# ax.plot(d1[11][1], d1[11][0], '-', linewidth=1.5, color=colours[1],zorder=2)
# ax.plot(d1[89][1], d1[89][0], '-', linewidth=1.5, color=colours[3],zorder=2)
try:
    ax.plot(d1[99][1], d1[99][0], '-', linewidth=1.5, color=colours[2],zorder=2)
except: 
    pass
# ax.imshow(blacks,extent=[0, 600, 400, 0],zorder=3)
ax.imshow(blacks,extent=[0, 600, 400, 0],zorder=3,interpolation='none')
ax.set_ylim((100,300))
ax.set_xlim((0,350))
plt.gca().invert_yaxis()
# SOLVED: print "Problem 1, you get a low resolution, probably interlaced picture as output..."
# SOLVED print "Problem 2, the picture is upside down -> rotate along horizontal axis"
# SOLVED print "Problem 3, it would be nice to get the grays a little lighter"
print "Problem 4, it would be cooool to add on top a small plot that tells you the average amount of gradient present (polar)"

# plt.show()
plt.savefig(figname)
sys.exit(1)



