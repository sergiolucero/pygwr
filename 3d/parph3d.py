import os
import sys
from time import time
from gwrlib3d import startph3d, phIteration3d, wrapitup3d
########################       initialization. 
quiet = 0
if len(sys.argv)<6:
    print 'USAGE: parph3d nx ny nz nscen nwells'
    exit()
else:
  nx,ny,nz,nscen,nwells = [int(sys.argv[i+1]) for i in range(5)]
###################################   setting up dirs+software
os.system('make -s')  # updates propagator!
os.system('easy_install --quiet mpi4py'); # if not already available!
os.system('apt-get -qq install zip');
if not(os.path.isdir('work')): os.mkdir('work');
if not(os.path.isdir('Presolved')): os.mkdir('Presolved');
###################################################
t1=time();k=0
nvec=[nx,ny,nz]
startph3d(nvec,nscen,nwells)		# now done in parallel
k,Efk,zk = phIteration3d(0,nscen,0,quiet)  # k=0, rho=0
z=[zk];Ef=[Efk]
# loading existing W+qhat could be an alternative
#then update qhat,W.... and maybe rho???
while (zk>0.01):
   k,Efk,zk = phIteration3d(0,nscen)
   z.append(zk); Ef.append(Efk)
##################################################
wrapitup3d(Ef,z,nvec,t1)
