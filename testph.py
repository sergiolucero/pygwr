import sys
from time import time
import os
from gwrlib import startph, phIteration, wrapitup
from numpy import savez
########################       initialization. 
nxvals=[2**i for i in range(8,12)]
print nxvals
nscen=4
nwells=4
quiet = 0

for nx in nxvals:
    t1=time();k=0
    startph(nx,nscen,nwells)		# now done in parallel
    k,Efk,zk = phIteration(0,nscen,0,quiet)  # k=0, rho=0
    print 'Tiempo[%d] =%g' %(nx,time()-t1)
print 'ALL DONE!'
