#!/usr/bin/python2.7

'''
take a key lock file and returns jvalues
'''

import sys,math,os,subprocess
#from PIL import Image
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from matplotlib.lines import Line2D
import numpy as np

def Distance(current,target):
  return sum([ abs(x-y) for x,y in zip(target,current) ])

def PrintSuccess(target,current,k1,l1,k2,l2):
  print
  print "Success"
  PrintResult(target,current,k1,l1,k2,l2)
  
def Print5000(target,current,k1,l1,k2,l2):
  print
  print howmanysteps,"steps passed, this is what I got"
  PrintResult(target,current,k1,l1,k2,l2)
  
def PrintResult(target,current,k1,l1,k2,l2):
  print "targ:", target
  print "curr:",current
  print "distance", Distance(current,target)
  print "with following K1,L1,K2,L2 (for KL...dat file)"
  print 
  print '3'
  for x in k1:
    print x,
  print  
  for x in l1:
    print x,
  print 
  for x in k2:
    print x,
  print  
  for x in l2:
    print x,
  print 


def ReverseKeyForMedium(hjm,lenkl=10,keypos_formedium=4):
  #reverse half jvals with medium
  # inverse formula is: Binary(hj1m - 10)|_4chars
  # so if we start from half j val 15 -> 5 -> 0101
  # of course no val higher than 10+Ten(1111)= 25
  # and min 10+Ten(0000) = 10
  if hjm<10 or hjm>25:
    print "hjm out of range [10,25]"
    sys.exit(1)
  
  key=lenkl*[0]
  hjm-=10
  i=keypos_formedium-1
  while i>=0:
    #print i,key
    key[i]=hjm%2
    hjm=hjm/2
    i-=1
  return key
    

def JWithMedium(key):
  Jval=0
  for i in range(keypos_formedium):
    Jval += int(key[i])*pow(2,keypos_formedium-i-1); #so that zeroth bit is most significant
  Jval += 10; #so that interaction with medium can't be 0
  return Jval

def JWithOtherTau(( key1,lock1 ),( key2,lock2 )):
  score=0;
  
  for i in range(len(key1)):
    score += 1 if key1[i] != lock2[i] else 0;
    score += 1 if key2[i] != lock1[i] else 0;
  
  Jval = 3 + (int)(0.5+ 40.*math.exp( -0.01 * pow(float(score),2.) ));
  return Jval
  
#filename="data_cellcount.txt"

def RandomiseKL(k1,l1,k2,l2,mutrate):

  rk1=[ 1-x if np.random.rand()<mutrate else x for i,x in enumerate(k1) ] # randomises everything
  rk2=[ 1-x if np.random.rand()<mutrate else x for i,x in enumerate(k2) ] # randomises everything
  #rk1=[ 1-x if (np.random.rand()<mutrate and i>=keypos_formedium) else x for i,x in enumerate(k1) ]
  #rk2=[ 1-x if (np.random.rand()<mutrate and i>=keypos_formedium) else x for i,x in enumerate(k2) ]
  
  rl1=[ 1-x if np.random.rand()<mutrate else x for i,x in enumerate(l1) ]
  rl2=[ 1-x if np.random.rand()<mutrate else x for i,x in enumerate(l2) ]
  return rk1,rl1,rk2,rl2

#######################################
####### BEGIN STUFF ###################
#######################################
j11 = 6 # to 41 #int(sys.argv[1])  # variable of interest

while j11<43:    
    sys.stderr.write("j11="+str(j11)+"\n")
    keypos_formedium=4
    lenkl=10

    hj1m = 10 #up to 15
    j12 = 8 # to 14  #high adh with prey - always
    hj2m = hj1m
    j21 = j12
    j22 = 38 # up to 41

    mutrate=0.025
    howmanysteps=10000

    k1=ReverseKeyForMedium(hj1m)
    k2=ReverseKeyForMedium(hj2m)

    k1=[ int(2.*np.random.rand()) if i>=keypos_formedium else x for i,x in enumerate(k1) ]
    k2=[ int(2.*np.random.rand()) if i>=keypos_formedium else x for i,x in enumerate(k2) ]

    l1=[ int(2.*np.random.rand()) for _ in range(lenkl) ]
    l2=[ int(2.*np.random.rand()) for _ in range(lenkl) ]

    #print k1,k2,l1,l2
    #exit(1)

    #print k1
    #print k2
    #print l1
    #print l2

    #Now the cool part:
    # calculate product with vector made by (hj1m,j11 ,j12 ,hj2m ,j21 ,j22 )
    # if it is very close, ok, else randomly change something, etc...

    # no
    #the max possible score is sum(hj1m^2,j11^2 ,j12^2 , etc...)
    #maxval= sum( [hj1m**2,j11**2 ,j12**2 ,hj2m **2,j21**2 ,j22**2] )

    init_targetJ = [hj1m,j11 ,j12 ,hj2m ,j21 ,j22 ]
    print "Target: ", [hj1m,j11 ,j12 ,hj2m ,j21 ,j22 ], "pos 0 -> +5, pos 2,4 -> +6, pos 5 -> +3"
    currentJ = [ JWithMedium(k1), 
              JWithOtherTau(( k1,l1 ),( k1,l1 )), 
              JWithOtherTau(( k1,l1 ),( k2,l2 )), 
              JWithMedium(k2),
              JWithOtherTau(( k1,l1 ),( k2,l2 )),
              JWithOtherTau(( k2,l2 ),( k2,l2 ))
              ]
      
    distance = Distance(init_targetJ,currentJ)

    counter=0
    counter2=0

    targetJ=init_targetJ[:]

    lresults=[]

    for j in [j12+x for x in range(7)]:
      targetJ=init_targetJ[:]
      targetJ[2]=j
      targetJ[4]=j
      
      for i in [hj1m+x for x in range(6)]:
        #target=init_target[:]
        targetJ[0]=i
        targetJ[3]=i
        #k1=ReverseKeyForMedium(hj1m)
        #k2=ReverseKeyForMedium(hj2m)
        #k1=[ int(2.*np.random.rand()) if i>=keypos_formedium else x for i,x in enumerate(k1) ]
        #k2=[ int(2.*np.random.rand()) if i>=keypos_formedium else x for i,x in enumerate(k2) ]
        
        for h in [j22+x for x in range(6)]:
          targetJ=targetJ[:-1]+[h]
        
          
          counter=0
          
          #print "Target:",targetJ
          
          distance=Distance(targetJ,currentJ) 
          
          while True:
              #if counter%100==0: print counter,
              if distance==0:
                #PrintSuccess(targetJ,currentJ,k1,l1,k2,l2)
                lresults.append( [ distance, currentJ, k1,l1,k2,l2] )
                break
                
              if counter>=howmanysteps:
                #Print5000(targetJ,currentJ,k1,l1,k2,l2)
                lresults.append( [ distance, currentJ, k1,l1,k2,l2] )
                #print "Got :", distance, currentJ
                #print
                break
              
              counter+=1
              #Randomise 
              rk1,rl1,rk2,rl2 = RandomiseKL(k1,l1,k2,l2,mutrate)
              while (rk1,rl1,rk2,rl2) == (k1,l1,k2,l2):
                rk1,rl1,rk2,rl2 = RandomiseKL(k1,l1,k2,l2,mutrate)
              #print k1,l1,k2,l2
              #rk1,rl1,rk2,rl2 = Randomise(k1,l1,k2,l2,mutrate)
              #print rk1,rl1,rk2,rl2
              #print 
              rcurrentJ = [ JWithMedium(rk1), 
                          JWithOtherTau(( rk1,rl1 ),( rk1,rl1 )), 
                          JWithOtherTau(( rk1,rl1 ),( rk2,rl2 )), 
                          JWithMedium(rk2),
                          JWithOtherTau(( rk1,rl1 ),( rk2,rl2 )),
                          JWithOtherTau(( rk2,rl2 ),( rk2,rl2 ))
                        ]
              
              rdistance = Distance(targetJ,rcurrentJ) 
              
              #if rdistance-distance<=4:
              # if rdistance<=distance or np.random.rand() < 0.05:
              #if rdistance<=distance or (rdistance-distance<=4 and np.random.rand() < 0.01):
              #print rdistance,counter 
              if rdistance<=distance:
                #print 'Got a better or equal:'
                #print currentJ
                #print rcurrentJ
                #print 'end'
                #counter+=1
                #if rdistance>=distance: 
                  #counter+=1
                #else: counter=0
                k1=rk1[:]
                k2=rk2[:]
                l1=rl1[:]
                l2=rl2[:]
                distance=rdistance
                currentJ=rcurrentJ[:]
                #print "Got new vector"
                #print "Target", target
                #print "new ve", current
                
                
    #print
    best=[[100000,'bla']]
    for result in lresults:
      if result[0]<best[0][0]:
        best=[result[:]]
      elif result[0]==best[0][0]:
        best.append( result[:] )
        

    for bla in best:
      print bla
      ## [distance, currentJ, k1,l1,k2,l2]
      #print "Distance:", result[0]
      #print result[1]
      #print
      #print '3'
      #for x in result[2]:
        #print x,
      #print
      #for x in result[3]:
        #print x,
      #print
      #for x in result[4]:
        #print x,
      #print
      #for x in result[5]:
        #print x,
      #print
      #print "--------------------"
      #print
      #print distance
      #raw_input("Press Enter for another cycle...")
    
    j11+=1
