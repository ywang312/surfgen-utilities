import numpy as np
import os.path
import sys

fdip1='propls.ci.drt1.state1.all'
fdip2='propls.ci.drt1.state2.all'
fdip12='trncils.FROMdrt1.state1TOdrt1.state2.all'

#Check if dipole file exist in the same directory
if os.path.isfile(fdip1) == False or os.path.isfile(fdip2) == False or os.path.isfile(fdip12) == False:
   print("Dipole file missing.")
   sys.exit()

def readdipole(file,m,n):
    data=[]
    with open(file,'r') as inpf:
     for line in inpf:
        strarray=line.split()
        for j in range(m,n):
            data.append(float(strarray[j]))
    return data

def diff(d1,d2):
    q=[]
    for i in range(0,90):
        q.append(float((d1[i]-d2[i])/2))
    return q

#AtD transformation
def rotateangle(diff,coup):
    theta=[]
    for i in range(0,30):
        nsqdiff=diff[3*i]**2+diff[3*i+1]**2+diff[3*i+2]**2
        nsqcoup=coup[3*i]**2+coup[3*i+1]**2+coup[3*i+2]**2
        tan4=2*(diff[3*i]*coup[3*i]+diff[3*i+1]*coup[3*i+1]+diff[3*i+2]*coup[3*i+2])/(nsqdiff-nsqcoup)
        rot=np.arctan(tan4)/4
        theta.append(float(rot))
    return theta

#Write internal gradient
def writeintgrd(theta,file):
    intgrd=[]
    with open(file,'w') as outf:
      for i in range(0,15):
          intgrd.append((theta[2*i+1]-theta[2*i])/0.002)
          outf.write(str('%16.9f'%intgrd[i]))
          outf.write("\n")

dipole1=readdipole(fdip1,1,4)
dipole2=readdipole(fdip2,1,4)
dipole12=readdipole(fdip12,2,5)
q12=diff(dipole1,dipole2)
theta=rotateangle(q12,dipole12)
writeintgrd(theta,'intgrd.dipole')
