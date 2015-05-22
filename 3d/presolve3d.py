import sys
from mpi4py import MPI
from gwrlib3d import *
from numpy import load
import warnings
from scipy.sparse import SparseEfficiencyWarning
warnings.simplefilter('ignore', SparseEfficiencyWarning)
##################################################################
# one scenario run, getting scenario info from MPI rank
nargs = len(sys.argv)
rho,nscen,nwells,nx,ny,nz,LB,UB = load('work/refcard3d.npy');
comm = MPI.COMM_WORLD; rank = comm.Get_rank()
print "[PreSolving scenario %d][nx=%d,ny=%d,nz=%d,nscen=%d,nwells=%d]" %(rank,nx,ny,nz,nscen,nwells)
GenHydroScen3d([nx,ny,nz],nscen,nwells,rank)
