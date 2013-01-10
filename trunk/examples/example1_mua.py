'''This script runs a standard LPA on MUA data'''

# standard libs
import os
import numpy as np
from scipy.io import loadmat
import cPickle

# pyLPA stuff
from pyLPA import LPA_Signal

# loading the underlying data matrix and define names
# LOADNAME = '/home/enorhe/Work/LPA/Data/LPAdata/R14P03S2Tones/'\
#     + 'R14P03S2TonesMUAAVG_flat'
LOADNAME = 'testdata.mat'
CASENAME = 'testcase'

DATA = loadmat(LOADNAME)
SETNAME_MUA = 'mua'
SETNAME_LFP = 'lfp'

# MUA = DATA[SETNAME_MUA]
MUA = DATA[SETNAME_MUA][SETNAME_MUA][0,0][:,:100,:]
LFP = DATA[SETNAME_LFP][SETNAME_LFP][0,0][:,:100,:]

###################################################################
# initializing the LPA_signal

nstim, ntime, nchan = MUA.shape

z_start = 0.   # depth of first electrode in mm
z_space = 0.1  # electrode spacing in mm

S_dict = {
    'mua_data' : MUA,
    'lfp_data' : LFP,
    'dt' : 1,
    'z_start' : z_start,
    'z_space' : z_space,
    'casename' : CASENAME,
    'tstim' : 15,
    'sub_at' : 'base', # 
    'verbose' : True
    }

lpa_signal = LPA_Signal(**S_dict)

####################################################################
# Initial guess

npop = 4

x0 = np.zeros(12)
x0[:npop] = [0.4, 0.7, 0.95, 1.4] # Lateral position of populations (mm)
x0[npop:2*npop] = [0.1, 0.1, 0.1, 0.1] # Width of populations
x0[2*npop:] = [0.1, 0.1, 0.1, 0.1] # Slope of populations

####################################################################
# set up arguments for solving

maxpos = z_start + (nchan - 1) * z_space
maxpopwidth = 0.7
maxslopewidth = z_space

# lower bounds
lb = npop * [z_start] + 2 * npop * [0]  
# upper bounds
ub = npop * [maxpos] + npop * [maxpopwidth] + npop * [maxslopewidth]

# put any arguments passed to the initialization of the solver here,
# see documentation for openopt.NLP for details.
init_args = {
    # 'maxIter' : 100,
    # 'maxFunEvals' : 1000,
    # 'maxTime' : 60,
    # 'maxCPUTime' : 100
    }

# There are also some additional arguments that could be passed to the
# call function of the solver, but most of this will be handled inside
# the call function for lpa_signal (like *x0* and *plot*)
solve_args = {}

# If you subclassed LPA_Signal() and wrote your own error function
# that needs additional arguments, you can put it in the *f_args*
# tuple. Remember that the first argument to the error function is
# always the parameters to be fitted, but this is handled under the
# hood. Only additional (only positional) arguments goes into *f_args*
f_args = ()

##################################################################
# Solve mua
solve_dict = {
    'init_args' : init_args,
    'solve_args' : solve_args,
    'f_args' : f_args,
    'plot' : False
    }
mode = 'mua'
solver = 'pswarm'
r, Mmat, rmat, Mphi = lpa_signal(mode, solver, x0, lb, ub, **solve_dict)

##################################################################
# Save the stuff

save_lpa = 'test_lpa.p'
save_r = 'test_r.p'

with open(save_lpa, 'w') as fid:
    lpa_signal.rmat = rmat
    lpa_signal.Mphi = Mphi
    cPickle.dump(lpa_signal, fid)
with open(save_r, 'w') as fid:
    cPickle.dump(r, fid)
