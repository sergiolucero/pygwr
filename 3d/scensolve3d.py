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
if nargs>1:    rho = double(sys.argv[1])
comm = MPI.COMM_WORLD; rank = comm.Get_rank()
phScenSolve3d(rank)
