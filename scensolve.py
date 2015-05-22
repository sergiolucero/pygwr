import sys
from mpi4py import MPI
from gwrlib import *
from numpy import load
import warnings
from scipy.sparse import SparseEfficiencyWarning
warnings.simplefilter('ignore', SparseEfficiencyWarning)
##################################################################
# one scenario run, getting scenario info from MPI rank
nargs = len(sys.argv)
if nargs>1:
    rho = double(sys.argv[1])
else:
    rho,nscen,nwells,nx,LB,UB = readref()
comm = MPI.COMM_WORLD; rank = comm.Get_rank()
phScenSolve(rank,rho)
