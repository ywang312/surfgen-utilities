import numpy as np
import os.path
import sys

dip1='propls.ci.drt1.state1.all'
dip2='propls.ci.drt1.state2.all'
dip12='trncils.FROMdrt1.state1TOdrt1.state2.all'

#Check if dipole file exist in the same directory

if os.path.isfile(dip1) == False or os.path.isfile(dip2) == False or os.path.isfile(dip12) == False:
   print('Dipole file missing.')
   sys.exit()
else:
   print('Dipole moment found, perform norm and dot product calculation...')

def readdipole(file,m,n):
    data=[]
    with open(file,'r') as inpf:
       for line in inpf:
          strarray=line.split()
          for j in range(m,n):
              data.append(float(strarray[j]))
    lendata=len(data)
    print('Number of dipole data points read is %3d'%(lendata/3))
    return data

def diff(d1,d2):
    print('Perform difference vector calculation...')
    q=[]
    for i in range(0,3):
        q.append(((d1[i]-d2[i])/2))
    print('the dipole difference vector is (%8.6f %8.6f %8.6f)'%(q[0],q[1],q[2]))
    print('the norm of the dipole difference vector is '+str(np.linalg.norm(q)))
    return q

def loopoverlines(dipole1,dipole2,dipole12):
    length=len(dipole1)/3
    if length!=len(dipole2)/3 or length!=len(dipole12)/3:
       print('Inproper number of dipole data read!')
       sys.exit()
    for i in range(int(length)):
        print('Enter point %3d calculation'%(i+1))
        d1=dipole1[3*i:3*i+3:1]
        d2=dipole2[3*i:3*i+3:1]
        d12=dipole12[3*i:3*i+3:1]
        print(d1,d2)
        print('----------')
        q12=diff(d1,d2)
        print('----------')
        print('the transition dipole is (%8.6f %8.6f %8.6f)'%(d12[0],d12[1],d12[2]))
        print('the norm of the transition dipole is '+str(np.linalg.norm(d12)))
        print('----------')
        print('the dot product of difference dipole and transition dipole is '+str(np.dot(d12,q12)))
        print('====================================================')
    return None 
#Write internal gradient

print('====================================================')
print('==============Beginning dipolesp.py...==============')
print('====================================================')
print('Reading dipole file in the same directory...')
dipole1=readdipole(dip1,1,4)
dipole2=readdipole(dip2,1,4)
dipole12=readdipole(dip12,2,5)
loopoverlines(dipole1,dipole2,dipole12)
print('Normal termination')
