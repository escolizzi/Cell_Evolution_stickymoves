
'''
contains standard adhesion module for cell evolution type experimentsself.
In the future it should also be able to read par files
'''

import sys,math,os,subprocess
import numpy as np

lenkl=24
lookuptable_Jmedium=[4,3,2,1,1,1]
keypos_formedium=len(lookuptable_Jmedium)

def Distance(current,target):
  return sum([ abs(x-y) for x,y in zip(target,current) ])

# in this new version we just use the lookup table    
def ReverseKeyForMedium(hjm,lenkl,lookuptable_Jmedium):
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

def JWithMedium(key,lookuptable_Jmedium):
  #Jval=0
  #for i,j in xrange(len(lookuptable_Jmedium)):
    #Jval += int(key[i])*pow(2,keypos_formedium-i-1); #so that zeroth bit is most significant
  Jval = sum( [ x*y  for x,y in  zip( key[:len(lookuptable_Jmedium)] , lookuptable_Jmedium  ) ]  )
  Jval += 8; #so that interaction with medium can't be 0
  return Jval

####   GENERALISED FORMULA   ####
def JWithOtherTau(( key1,lock1 ),( key2,lock2 )):
  score=0;
  
  for i in range(len(key1)):
    score += 1 if key1[i] != lock2[i] else 0;
    score += 1 if key2[i] != lock1[i] else 0;
  
  #Jval = 3 + (int)(0.5+ 40.*math.exp( -pow(float(score),2.) / pow( float(lenkl) , 2.) ));
  #the above formula has been changed to 
  Jval =(int)(0.5+  52.- 48./float(2*lenkl) *float(score))
  
  return Jval
  
def RandomKL(lenkl):
  rk1=[ 1 if np.random.rand()<0.5 else 0 for _ in range(lenkl) ] # randomises everything
  rk1=[ 1 if np.random.rand()<0.5 else 0 for _ in range(lenkl) ]
  return rk1,rl1
