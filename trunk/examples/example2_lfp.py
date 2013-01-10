'''This script runs a standard LPA on MUA data'''

# standard libs
import os
import numpy as np
from scipy.io import loadmat
import cPickle


#############################################################
# Load results from example1_mua.py

save_lpa = 'test_lpa.p'
save_r = 'test_r.p'

with open(save_lpa, 'r') as fid:
    lpa_signal = cPickle.load(fid)
with open(save_r, 'r') as fid:
    r = cPickle.load(fid)

# rmat and MPhi does not belong in lpa_signal, just put it there for
# safekeeping.
rmat = lpa_signal.rmat
Mphi = lpa_signal.Mphi
del lpa_signal.rmat
del lpa_signal.Mphi


###########################################################
# Set up arguments for solving lfp

kernel = 'singleExp'
x0 = [50., 0.1]

lb = [np.finfo(np.double).tiny, 0]
ub = [50, 50]

# these dicts are explained in example1_mua.py
init_args = {}
solve_args = {}

f_args = (rmat, kernel)

###############################################################
# Solve lfp
solve_dict = {
    'init_args' : init_args,
    'solve_args' : solve_args,
    'f_args' : f_args,
    'plot' : False
    }

mode = 'lfp'
solver = 'scipy_cobyla'
r, Lmat, Rmat, Lphi = lpa_signal(mode, solver, x0, lb, ub, **solve_dict)

