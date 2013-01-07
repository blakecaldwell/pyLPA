'''
This file should contain most tools used for LPA.
'''

import pylab as pl

try:
    import openopt as oopt
    oopt_import = True
except ImportError:
    oopt_import = False

oopt_import_msg = "WARNING: No module named 'openopt'. LPA must " \
    "be done without using 'fmin_oopt()' method and associated "\
    "solvers. This means that only 'randw' algorithm is available."


glp_solvers = ['galileo' ,'pswarm', 'de', 'stogo', 'isres', 'mlsl']
nlp_solvers = ['ralg', 'scipy_cobyla', 'algencan', 'scipy_slsqp',
               'mma', 'auglag', 'ipopt', 'lincher', 'scipy_fmin',
               'scipy_lbfgsb']


class LPA_Signal(object):
    '''
    Main class for Laminar Population Analysis
    '''
    def __init__( self, mua_data, dt, lfp_data=None,
                  z_start=0., z_space=0.1, casename='', rneg_factor=1E6,
                  tstim=0, sub_at='base',
                  verbose=False ):
        '''
        This function ...
    
        Aguments
        --------

        Keyword arguments
        -----------------

        '''
        
        if verbose:
            msg = 'This is class *LPA_Signal* in *pyLPA* module'
            print msg
        
        if not oopt_import:
            print oopt_import_msg
            
        self.nstim, self.ntime, self.nchan = pl.asarray(mua_data).shape

        self.dt = dt
        self.z_start = z_start
        self.z_space = z_space
        self.el_coords = pl.arange(z_start, z_start+self.nchan*z_space, z_space)
        self.rneg_factor = rneg_factor
        self.tstim = tstim
        self.sub_at = sub_at
        self.verbose = verbose

        self.tstim_idx = pl.floor(self.tstim/self.dt)
        
        # create matrices and calculate variances
        self.ImportDataset(mua_data, 'MUA')
        if not lfp_data==None:
            self.ImportDataset(lfp_data, 'LFP')

    def ImportDataset( self, lpa_data, mode ):
        '''
        Reshapes matrices and comupte baseline variance and adds to attributes
        This function ...
    
        Aguments
        --------

        Keyword arguments
        -----------------


        '''
        lpadat = pl.asarray(lpa_data.copy())
        matname = mode.lower() + 'mat'
        varname = mode.lower() + 'var'

        nstim, ntime, nchan = lpadat.shape

        if mode.lower()=='lfp' and hasattr(self, '_muamat'):
            if nstim != self.nstim:
                raise Exception, 'Number of stimuli in %s and MUA data' \
                    ' does not match' % setname
            elif ntime != self.ntime+self.tstim_idx:
                raise Exception, 'Number of sampling points in %s and' \
                    ' MUA data does not match' % setname
            elif nchan != (self.maxchannel - self.minchannel +1):
                raise Exception, 'Number of channels in %s and MUA data' \
                    ' does not match' % setname
        # Apply base/mean subtraction to each stimulus and channel separately
        tmp_idx = 0
        if self.sub_at=='base':
            tmp_idx = self.tstim_idx
        elif self.sub_at=='mean':
            tmp_idx = ntime
        else:
            msg = '%s is not a valid choice for sub_at' % self.sub_at
            
        if not tmp_idx==0:
            for istim in xrange(self.nstim):
                for ichan in xrange(self.nchan):
                    lpadat[istim, :, ichan] = lpadat[istim, :, ichan] \
                        - lpadat[istim, :tmp_idx, ichan].mean()

        # reshape data to 2D
        lpamat = lpadat.reshape((self.nstim*self.ntime, self.nchan)).transpose()
        
        # Evaluate variances of stimulus evoked signal
        lpavar = lpadat[:, self.tstim_idx:, :].var()

        exec('self._'+matname+' = lpamat')
        exec('self._'+varname+' = lpavar')

    def __call__( self, mode, solver, x0, init_args={}, solve_args={},
                  f_args=(), plot=False ):
        '''
        This is where the action starts.
    
        Aguments
        --------

        Keyword arguments
        -----------------
        '''
        
        if self.verbose:
            print mode
        # Check if initial guess is provided in one or both of the argument
        # dictionaries
        redundant_args = ['x0', 'plot']
        for red_arg in redundant_args:
            if red_arg in init_args.keys():
                msg = 'Initial guess *%s* found in *init_args*. I will ' \
                    'override this value with the value *%s* passed on as '\
                    'positional/keyword argument' % (red_arg, red_arg)
                if self.verbose:
                    print msg
                del init_args[red_arg]
            if red_arg in solve_args.keys():
                msg = 'Initial guess *%s* found in *solve_args*. I will ' \
                    'override this value with the value *%s* passed on as '\
                    'positional/keyword argument' % (red_arg, red_arg)
                if self.verbose:
                    print msg
                del solve_args[red_arg]        

        if mode.lower()=='mua':
            # Check for consistent length of x0
            M_params = pl.asarray(x0).squeeze()
            if pl.mod(M_params.shape[0], 3):
                err_msg = 'Number of elements in x0 must be dividable' \
                    'by 3. Current number of elements is %d' % M_params.shape[0]
                raise Exception, err_msg
            else:
                self.npop = M_params.shape[0]/3

            # create linear constraints
            A, b = self._create_constraints()
            init_args_internal = {'A' : A, 'b' : b}
            init_args_internal.update(init_args)
            
            # Set some attributes
            self.mua_solver = solver
            self.mua_x0 = x0
            self.mua_init_args = init_args_internal
            self.mua_solve_args = solve_args            
            
            # Set the error function and arguments for it. Default is no 
            # arguments except the parameters used for fitting, but this
            # might change in a subclass
            f_fnc = self._mua_err
            # f_args = self._mua_f_args(**f_fnc_args)
            
            # Finally, put everything into a dictionary that can be 
            # passed to self._solve()
            fmin_args = {
                'solver' : self.mua_solver,
                'f_fnc' : f_fnc,
                'f_args' : f_args,
                'x0' : self.mua_x0,
                'init_args' : self.mua_init_args,
                'solve_args' : self.mua_solve_args,
                'plot' : plot}
            
            r = self._solve(fmin_args)

        elif mode.lower()=='lfp':
            self.lfp_solver = solver
            self.lfp_x0 = x0
            self.lfp_solver_args = solver_args
            # put in lfp f_fnc and args

        # Post solver processing
        return r

    def _solve( self, fmin_args ):
        '''
        This function ...
    
        Aguments
        --------

        Keyword arguments
        -----------------
        '''

        if fmin_args['solver'] == 'fmin_randw':
            r = self._fmin_randw()
        elif not oopt_import:
            raise Exception, oopt_import_msg
        elif fmin_args['solver'] in nlp_solvers+glp_solvers:
            r = self._fmin_oopt(**fmin_args)
        else:
            pass

        return r
    
    def _fmin_randw( self ):
        '''
        This function ...
    
        Aguments
        --------

        Keyword arguments
        -----------------
        '''

        pass

    def _fmin_oopt( self, f_fnc, f_args, solver, x0, init_args={},
                    solve_args={}, plot=False):
        '''
        This function ...
    
        Aguments
        --------

        Keyword arguments
        -----------------
        '''
        
        fnc = lambda x: f_fnc(x, *f_args)

        # Check what solver is in use and put initial guess argument in right
        # dictionary
        print solver
        if solver in glp_solvers:
            solve_args['plot'] = plot
            solve_args['x0'] = x0
            
            p = oopt.GLP(fnc, **init_args) # Set up openopt class
            r = p.solve(solver, **solve_args) # and solve

        elif solver in nlp_solvers:
            init_args['plot'] = plot
            init_args['x0'] = x0

            p = oopt.NLP(fnc, **init_args) # Set up openopt class
            r = p.solve(solver, **oopt_solver_args) # and solve
            
        else: 
            raise Exception, 'The solver %s is not recognized. Check' \
                'spelling!' % solver

        
        return r

    def _mua_err( self, M_params ):
        '''
        This function ...
    
        Aguments
        --------

        Keyword arguments
        -----------------
        '''
        
        Mmat, rmat = self.create_Mmat_rmat( M_params )
        ff = err_eval(self._muamat, Mmat, rmat)

        return ff
        
#    def _mua_f_args( self ):
#        '''
#        This function ...
#    
#        Aguments
#        --------
#
#        Keyword arguments
#        -----------------
#        '''
#        
#        f_args = ()
#
#        return f_args
    
    def lfp_err( self ):
        '''
        This function ...
    
        Aguments
        --------

        Keyword arguments
        -----------------
        '''

        pass

    def create_Mmat_rmat( self, M_params):
        '''
        This function ...
    
        Aguments
        --------
        
        Keyword arguments
        -----------------
        '''
        Mmat = _create_Mmat( M_params, self.el_coords )
        rmat = _create_rmat( Mmat, self._muamat, self._muavar,
                             self.rneg_factor )
        
        return Mmat, rmat
    
    def _create_constraints( self ):
        '''
        This function ...
    
        Aguments
        --------
        
        Keyword arguments
        -----------------
        '''

        npop = self.npop
    
        # Create linear constraints
        A = pl.eye(3*npop)
        b = pl.zeros(3*npop)
        for i in xrange(npop-1):
            # setting extra constraints on z0
            A[i, i+1] = -1
            # setting constraints on pop width
            A[i + npop, i] = 1
            A[i + npop, i + 1] = -1
            A[i + npop, i + npop + 1] = 1
            # no additional constraints on slope
        # Treat A[-1] separately
        A[2*npop-1, npop - 2] = 1
        A[2*npop-1, npop - 1] = -1
        # b is all zeros, since box boundaries are treated with upper
        # and lower bounds
        b[npop-1] = self.el_coords[-1]
        # Try this with and without the maxslopewidth
        b[2*npop-1:] = 0.1
        
        return A, b

    
def err_eval( lpamat, Smat, Tmat):
    '''
    This function ...
    
    Aguments
    --------
        
    Keyword arguments
    -----------------
    '''
        
    lpamat_est = pl.dot(Smat, Tmat)
    lpamat_diff = lpamat - lpamat_est
    err = pl.mean(lpamat_diff**2)/pl.mean(lpamat**2)

    return err
    
def _create_Mmat( M_params, el_coords ):
    '''
    This function ...
    
    Aguments
    --------

    Keyword arguments
    -----------------
    '''
    # Assume now that all measures are real measures in mm
    M_params = pl.asarray(M_params)
        
    nchan = el_coords.size
    
    npop = M_params.shape[0]/3
    z0=M_params[:npop]
    a=M_params[npop:2*npop]
    b=M_params[2*npop:]                
        
    # Create the Mmat matrix
    Mmat = pl.zeros((nchan, npop))
            
    # Update the basis matrix, Mmat
    for ipop in xrange(npop):
        # tmpvec is the distance between z0[ipop] and all electrode contacts
        tmpvec = abs(el_coords - z0[ipop])
        # Set all entries that are closer than a[ipop] to 1
        Mmat[pl.find(tmpvec < a[ipop]), ipop] = 1
        # Find the entries that are on the olique part of the profile and set
        # the right value
        tmpvec = tmpvec - a[ipop]
        cond1 = pl.find(tmpvec >= 0)
        cond2 = pl.find(tmpvec < b[ipop])
        isect = filter(lambda x: x in cond1, cond2)
        Mmat[isect, ipop] = 1 - 1 / b[ipop] * (tmpvec[isect])
        
    return Mmat

def _create_rmat( Mmat, muamat, muavar, rneg_factor):
    '''
    This function ...
    
    Aguments
    --------
    
    Keyword arguments
    -----------------
    '''
    
    # Estimate rmat based on pseudoinverse of Mmat
    rmat=pl.dot(pl.linalg.pinv(Mmat), muamat)
    # Rectification of rmat
    rmat[pl.find(rmat<-rneg_factor*pl.sqrt(muavar))]= \
        -rneg_factor*pl.sqrt(muavar)
    
    return rmat

            
class randw(object):
    def __init__( self ):
        pass
    def __call__( self ):
        pass
