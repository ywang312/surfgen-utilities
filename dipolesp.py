import numpy as np
import os.path
import sys
import math

dip1='propls.ci.drt1.state1.all'
dip2='propls.ci.drt1.state2.all'
dip12='trncils.FROMdrt1.state1TOdrt1.state2.all'

#Check if dipole file exist in the same directory
if os.path.isfile(dip1) == False or os.path.isfile(dip2) == False or os.path.isfile(dip12) == False:
   print('Dipole file missing!')
   sys.exit()
else:
   print('Dipole moment found, perform norm and dot product calculation...')

#filename, block to read
def readdipole(file,m,n):
    data=[]
    with open(file,'r') as inpf:
       for line in inpf:
          strarray=line.split()
          for j in range(m,n):
              data.append(float(strarray[j]))
    pts=len(data)/3
    print('Number of dipole data points read is %3d'%(pts))
    return data

#G12=1/2*(O11-O22) sign doesn't matter here
def diff(d1,d2):
    q=np.zeros(3)
    for i in range(0,3):
        q[i]=(d1[i]-d2[i])/2
    print('the dipole difference vector is (%8.6f %8.6f %8.6f)'%(q[0],q[1],q[2]))
    return q

def loopoverlines(dipole1,dipole2,dipole12):
    length=len(dipole1)/3
    if length!=len(dipole2)/3 or length!=len(dipole12)/3:
       print('Inproper number of dipole data read!')
       sys.exit()
    print('Loading %d points...'%length)
    for i in range(int(length)):
        print('Enter point %3d calculation'%(i+1))
        #Partition the data array
        d1=dipole1[3*i:3*i+3:]
        d2=dipole2[3*i:3*i+3:]
        d12=dipole12[3*i:3*i+3:]
        print(d1,d2)
        print('----------')
        q12=diff(d1,d2)
        print('----------')
        print('the transition dipole is (%8.6f %8.6f %8.6f)'%(d12[0],d12[1],d12[2]))
        print('----------')
        print('the dot product of difference dipole and transition dipole is '+str(np.dot(d12,q12)))
        normsquarediff=(np.linalg.norm(q12))**2-(np.linalg.norm(d12))**2
        print('the square norm difference of difference dipole and transition dipole is '+str(normsquarediff))
        print('====================================================')
    return None 

def grd(dipole1,dipole2,dipole12):
    print('====================================================')
    print('Enable single point dipole gradient calculation...')
    ninternal=(len(dipole1)/3-1)/2
    print('System with %d internal coordinate detected'%int(ninternal))    
    #print g,h grd vector
    grdg=[]
    grdh=[]   
    q=np.zeros(6*int(ninternal))
    for i in range(6*int(ninternal)):
        q[i]=(dipole1[i+3]-dipole2[i+3])/2
    print('Enter g gradient calculation...') 
    with open('intgrd.g','w') as outg:     
        for i in range(int(ninternal)):
            print('reads in coord%d'%(i+1))
            d1b=dipole1[3+6*i:6+6*i:]
            d2b=dipole2[3+6*i:6+6*i:]
            d12b=dipole12[3+6*i:6+6*i:]
            q12b=q[6*i:3+6*i:]
            d1f=dipole1[6+6*i:9+6*i:]
            d2f=dipole2[6+6*i:9+6*i:]
            d12f=dipole12[6+6*i:9+6*i:]
            q12f=q[3+6*i:6+6*i:]
            gb=(np.linalg.norm(q12b))**2-(np.linalg.norm(d12b))**2
            gf=(np.linalg.norm(q12f))**2-(np.linalg.norm(d12f))**2
            grdg.append((gf-gb)/0.0002) #Assuming the interval is 0.0001
            outg.write(str('%16.9f'%grdg[i]))
            outg.write("\n")
    print('g gradient written to intgrd.g')
    print('Enter h gradient calculation...') 
    with open('intgrd.h','w') as outh:      
        for i in range(int(ninternal)):
            print('reads in coord%d'%(i+1))   
            d1b=dipole1[3+6*i:6+6*i:]
            d2b=dipole2[3+6*i:6+6*i:]
            d12b=dipole12[3+6*i:6+6*i:]
            q12b=q[6*i:3+6*i:]
            hb=np.dot(d12b,q12b)
            d1f=dipole1[6+6*i:9+6*i:]
            d2f=dipole2[6+6*i:9+6*i:]
            d12f=dipole12[6+6*i:9+6*i:]
            q12f=q[3+6*i:6+6*i:]
            hf=np.dot(d12f,q12f)
            grdh.append((hf-hb)/0.0002) 
            #Make sure the sign is the same for transition dipole
            outh.write(str('%16.9f'%grdh[i]))
            outh.write("\n")
    print('h gradient written to intgrd.h')
    return grdg,grdh

def dtheta(dipole1,dipole2,dipole12):
    print('optional calculation for the derivative coupling from AtD transform...')
    ninternal=(len(dipole1)/3-1)/2
    q=np.zeros(6*int(ninternal))
    norm=0
    for i in range(6*int(ninternal)):
        q[i]=(dipole1[i+3]-dipole2[i+3])/2
    with open('intgrd.nad.AtD','w') as outc:  
        for i in range(int(ninternal)):
            print('reads in coord%d'%(i+1))
            d1b=dipole1[3+6*i:6+6*i:]
            d2b=dipole2[3+6*i:6+6*i:]
            d12b=dipole12[3+6*i:6+6*i:]
            q12b=q[6*i:3+6*i:]
            d1f=dipole1[6+6*i:9+6*i:]
            d2f=dipole2[6+6*i:9+6*i:]
            d12f=dipole12[6+6*i:9+6*i:]
            q12f=q[3+6*i:6+6*i:]
            gb=(np.linalg.norm(q12b))**2-(np.linalg.norm(d12b))**2
            gf=(np.linalg.norm(q12f))**2-(np.linalg.norm(d12f))**2
            hb=np.dot(d12b,q12b)
            hf=np.dot(d12f,q12f)
            thetab=math.atan(2*hb/gb)*0.25
            thetaf=math.atan(2*hf/gf)*0.25
            grdtheta=(thetaf-thetab)/0.0002
            outc.write(str('%16.9f'%grdtheta))
            outc.write("\n")
            norm=norm+grdtheta**2
    print('the norm of derivative coupling is %16.9f'%norm)
    print('nad gradient from AtD written to intgrd.nad.AtD')
    return None



print('====================================================')
print('==============Beginning dipolesp.py...==============')
print('====================================================')
print('Reading dipole file in the same directory...')
dipole1=readdipole(dip1,1,4)
dipole2=readdipole(dip2,1,4)
dipole12=readdipole(dip12,2,5)
loopoverlines(dipole1,dipole2,dipole12)
grd(dipole1,dipole2,dipole12)
dtheta(dipole1,dipole2,dipole12)
print('Normal termination')
