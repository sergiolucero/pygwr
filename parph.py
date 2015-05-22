#/usr/bin/python
import sys
from numpy import random, save
from time import time
from gwrlib import softsetup, startph, phIteration, wrapitup
########################       initialization. 
quiet = 0
if len(sys.argv)<4:
    print 'USAGE: parph nx nscen nwells'
    exit()
else:
  nx,nscen,nwells = [int(sys.argv[i+1]) for i in range(3)]
###################################   setting up dirs+software
softsetup()
###################################################
#wvec = 100*random.rand(1,3); print wvec; 
wvec = [[ 11.7046, 70.8941,94.154 ]]
save('work/wvec.npy',wvec)
t1=time();k=0
startph(nx,nscen,nwells,-2,2,0.50)		# now done in parallel
#startph(nx,nscen,nwells)		# now done in parallel
k,Efk,zk = phIteration(0,nscen,0,quiet)  # k=0, rho=0
z=[zk];Ef=[Efk]
# loading existing W+qhat could be an alternative
#then update qhat,W.... and maybe rho???
#while (zk>0.0001):
while (zk>0.1):
   k,Efk,zk = phIteration(0,nscen)
   z.append(zk); Ef.append(Efk)
##################################################
#wrapitup(Ef,z,nx,t1)
