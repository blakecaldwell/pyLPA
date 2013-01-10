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
        self.importDataset(mua_data, 'MUA')
        if not lfp_data==None:
            self.importDataset(lfp_data, 'LFP')

    def importDataset( self, lpadata, mode ):
        '''
        Reshapes matrices and comupte baseline variance and adds to attributes
        This function ...
    
        Aguments
        --------

        Keyword arguments
        -----------------


        '''

        matname = mode.lower() + 'mat'
        varname = mode.lower() + 'var'

        nstim, ntime, nchan = lpadata.shape

        if mode.lower()=='lfp' and hasattr(self, '_muamat'):
            if nstim != self.nstim:
                raise Exception, 'Number of stimuli in %s and MUA data' \
                    ' does not match' % setname
            elif ntime != self.ntime:
                raise Exception, 'Number of sampling points in %s and' \
                    ' MUA data does not match' % mode
            elif nchan != self.nchan:
                raise Exception, 'Number of channels in %s and MUA data' \
                    ' does not match' % mode
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
                    lpadata[istim, :, ichan] = lpadata[istim, :, ichan] \
                        - lpadata[istim, :tmp_idx, ichan].mean()

        # reshape data to 2D
        lpamat = self._reshapeMat( lpadata )
        #lpamat = lpadata.reshape((self.nstim*self.ntime, self.nchan)).transpose()
        
        # Evaluate variances of stimulus evoked signal
        lpavar = lpadata[:, self.tstim_idx:, :].var()

        exec('self._'+matname+' = lpamat')
        exec('self._'+varname+' = lpavar')

    def __call__( self, mode, solver, x0, lb, ub, init_args={}, solve_args={},
                  f_args=(), plot=False ):
        '''
        This is where the action starts.
    
        Aguments
        --------

        Keyword arguments
        -----------------
        '''
        
        if self.verbose:
            msg = 'Solving for %s part of signal' % mode 
            print msg
        # Check if initial guess is provided in one or both of the argument
        # dictionaries
        redundant_args = ['A', 'b', 'lb', 'ub', 'x0', 'plot']
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

        # Set temporary attributes. These are removed in the end
        attr_list = ['solver', 'x0', 'ub', 'init_args', 'solve_args']
        for attr_name in attr_list:
            try:
                exec('self.'+attr_name+' = '+attr_name+'.copy()')
            except AttributeError:
                exec('self.'+attr_name+' = '+attr_name)

        if mode.lower()=='mua':
            # Check for consistent length of x0
            M_params = pl.asarray(x0).squeeze()
            if pl.mod(M_params.shape[0], 3):
                err_msg = 'Number of elements in x0 must be dividable' \
                    'by 3. Current number of elements is %d' % M_params.shape[0]
                raise Exception, err_msg
            else:
                self.npop = M_params.shape[0]/3
            
            # Do some pre-solve operations
            f_args = self._muaWarmup(*f_args)
                
            # create linear constraints
            A, b = self._muaConstraints()            
            
            # Set the error function
            f_fnc = self._muaErr
            
        elif mode.lower()=='lfp':
            # linear constraints
            f_args = self._lfpWarmup(*f_args)
            A, b = self._lfpConstraints()
            # Error function
            f_fnc = self._lfpErr

        init_args_internal = {
            'A' : A, 'b' : b,
            'lb' : lb, 'ub' : self.ub
            }
        
        # Override self.init_args with init_args_internal
        self.init_args.update(init_args_internal)

        # Finally, put everything into a dictionary that can be 
        # passed to self._solve()
        fmin_args = {
            'solver' : self.solver,
            'f_fnc' : f_fnc,
            'f_args' : f_args,
            'x0' : self.x0,
            'init_args' : self.init_args,
            'solve_args' : self.solve_args,
            'plot' : plot
            }
        
        r = self._solve(fmin_args)

        ############################################################
        # Post solver processing
        # Find decompositions and full matrix for best parameters
        if mode.lower()=='mua':
            Smat, Tmat_tmp = self.muaDecomp( r.xf, *f_args )
        elif mode.lower()=='lfp':
            Smat, Tmat_tmp = self.lfpDecomp( r.xf, *f_args )

        Phimat_tmp = pl.dot(Smat, Tmat_tmp)
        # Reshape to 3d
        Tmat = self._reshapeMat( Tmat_tmp )
        Phimat = self._reshapeMat( Phimat_tmp )
        
        # deleting temporary attributes
        for attr_name in attr_list:
            exec('del self.'+attr_name)        
        
        return r, Smat, Tmat, Phimat

    def _reshapeMat( self, rawmat ):
        '''
        This function ...
    
        Aguments
        --------

        Keyword arguments
        -----------------
        '''
        
        if rawmat.ndim==2:
            
            # Create the 3D scores
            outmat = rawmat.transpose().reshape( 
                (self.nstim, self.ntime, rawmat.shape[0]) )
            # Create full 3D array
#            Phimat = pl.asarray(
#                [pl.dot(Smat, Tm.transpose()).transpose() for Tm in Tmat])
        elif rawmat.ndim==3:
            outmat = rawmat.reshape(
                (self.nstim*self.ntime, rawmat.shape[-1]) ).transpose()
            

        return outmat
        

    def _solve( self, fmin_args ):
        '''
        This function ...
    
        Aguments
        --------

        Keyword arguments
        -----------------
        '''

        if fmin_args['solver'] == 'fmin_randw':
            r = self._fminRandw()
        elif not oopt_import:
            raise Exception, oopt_import_msg
        elif fmin_args['solver'] in nlp_solvers+glp_solvers:
            r = self._fminOopt(**fmin_args)
        else:
            pass

        return r
    
    def _fminRandw( self ):
        '''
        This function ...
    
        Aguments
        --------

        Keyword arguments
        -----------------
        '''

        pass

    def _fminOopt( self, f_fnc, f_args, solver, x0, init_args={},
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
        if solver in glp_solvers:
            solve_args['plot'] = plot
            solve_args['x0'] = x0
            
            p = oopt.GLP(fnc, **init_args) # Set up openopt class
            r = p.solve(solver, **solve_args) # and solve

        elif solver in nlp_solvers:
            init_args['plot'] = plot
            init_args['x0'] = x0

            p = oopt.NLP(fnc, **init_args) # Set up openopt class
            r = p.solve(solver, **solve_args) # and solve
            
        else: 
            raise Exception, 'The solver %s is not recognized. Check' \
                'spelling!' % solver

        
        return r

    def _muaErr( self, M_params ):
        '''
        This function ...
    
        Aguments
        --------

        Keyword arguments
        -----------------
        '''
        
        Mmat, rmat = self.muaDecomp( M_params )
        ff = self.errEval(self._muamat, Mmat, rmat)

        return ff
        
    def _lfpErr( self, L_params, rmat, kernel ):
        '''
        This function ...
    
        Aguments
        --------

        Keyword arguments
        -----------------
        '''
        
        Lmat, Rmat = self.lfpDecomp( L_params, rmat, kernel )
        ff = self.errEval(self._lfpmat, Lmat, Rmat)

        return ff
    
    def lfpDecomp( self, L_params, rmat, kernel ):
        '''
        This function ...
    
        Aguments
        --------
        
        Keyword arguments
        -----------------
        '''

        exec('h_list = _'+kernel+'(L_params, self.dt)')
        
        Rmat = _createRmat(h_list, rmat, self.nstim, self.ntime)
        Lmat = pl.dot(self._lfpmat, pl.linalg.pinv(Rmat))

        return Lmat, Rmat

    def muaDecomp( self, M_params):
        '''
        This function ...
    
        Aguments
        --------
        
        Keyword arguments
        -----------------
        '''

        Mmat = _createMmat( M_params, self.el_coords )
        rmat = _create_rmat( Mmat, self._muamat, self._muavar,
                             self.rneg_factor )
        
        return Mmat, rmat
    
    def _muaConstraints( self ):
        '''
        This function ...
    
        Aguments
        --------
        
        Keyword arguments
        -----------------
        '''
        
        npop = self.npop
        ub = self.ub
    
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
        # Treat the last pop width separately
        A[2*npop-1, npop - 2] = 1
        A[2*npop-1, npop - 1] = -1
        # b is set for the rows where no constraints are present
        b[npop-1] = ub[npop-1]
        # Try this with and without the maxslopewidth        
        b[2*npop:] = ub[2*npop:]
        
        return A, b

    def _lfpConstraints( self ):
        '''
        This function ...
    
        Aguments
        --------
    
        Keyword arguments
        -----------------
        '''
        A, b = (None, None)
        
        return A, b

    def _muaWarmup( self, *f_args ):
        '''
        This function ...
    
        Aguments
        --------
    
        Keyword arguments
        -----------------
        '''
        return f_args

    def _lfpWarmup( self, rmat, kernel):
        '''
        This function ...
    
        Aguments
        --------
    
        Keyword arguments
        -----------------
        '''
        if rmat.ndim==3:
            rmat = self._reshapeMat( rmat )
            
        return rmat, kernel
    
    def errEval( self, lpamat, Smat, Tmat):
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

def _createLmat(  ):
    '''
    This function ...
    
    Aguments
    --------
    
    Keyword arguments
    -----------------
    '''
    pass
    

def _createRmat( h_list, rmat, nstim, ntime ):
    '''
    This function ...
    
    Aguments
    --------
    
    Keyword arguments
    -----------------
    '''
    nsplit = len(h_list)    
    npop, _ = rmat.shape

    Rmat = pl.zeros((nsplit*npop, nstim*ntime))
    
    for istim in xrange(nstim):
        for isplit in xrange(nsplit):
            for ipop in xrange(npop):
                rvec = rmat[ipop, istim*ntime:(istim+1)*ntime]
                h = h_list[isplit]
                tmp_vec= pl.convolve(rvec, h)
                tmp_vec=tmp_vec[:ntime]
                Rmat[ipop+isplit*npop,\
                         istim*ntime:(istim+1)*ntime]=tmp_vec

    return Rmat

    
def _createMmat( M_params, el_coords ):
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

def _singleExp(x, dt):
    # Create the convolution kernel
    tau=float(x[0])
    Delta=float(x[1])
    t = pl.arange(0, (10*tau), dt)
    h = pl.zeros(len(t))
    
    h = 1/tau*pl.exp(-(t-Delta)/tau)
    h[pl.find(t<Delta)]=0

    return [h]

def _singleAlpha(x, dt):
    tau = float(x[0])
    Delta = float(x[1])
    t = pl.arange(0, (10*tau), dt)
    h = pl.zeros(len(t))
    h = (t-Delta)/tau**2*pl.exp(-(t-Delta)/tau)
    h[pl.find(t<Delta)] = 0

    return [h]

def _doubleExp(x, dt):
    h1 = _singleExp(x[:2], dt)
    h2 = _singleExp(x[2:], dt)
    h = [h1[0], h2[0]]

    return h

def _doubleAlpha(x, dt):
    h = []
    h.append(_singleAlpha(x[:2], dt))
    h.append(_singleAlpha(x[2:], dt))

    return h
            
class randw(object):
    def __init__( self ):
        pass
    def __call__( self ):
        pass
