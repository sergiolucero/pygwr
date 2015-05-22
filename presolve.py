import sys
from mpi4py import MPI
from gwrlib import GenHydroScen, readref
from numpy import load
import warnings
from scipy.sparse import SparseEfficiencyWarning
warnings.simplefilter('ignore', SparseEfficiencyWarning)
##################################################################
# one scenario run, getting scenario info from MPI rank
nargs = len(sys.argv)
rho,nscen,nwells,nx,LB,UB = readref()
#refcard = load('work/refcard.npy')
#rho,nscen,nwells,nx,LB,UB = readcard(refcard)
#print 'a,b,c,d=',nx,nscen,nwells,LB,UB,rho
comm = MPI.COMM_WORLD; rank = comm.Get_rank()
#print "[PreSolving scenario %d]" %(rank),
#print "[nx=%d,nscen=%d,nwells=%d]" %(nx,nscen,nwells)
#print rank
GenHydroScen(nx,nscen,nwells,rank)
