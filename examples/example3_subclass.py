from pyLPA import LPA_Signal

class LPA_mod(LPA_Signal):
    
    def __init__( self, extra_arg, *args, **kwargs):
        self.extra_arg = extra_arg
        LPA_Signal.__init__(self, *args, **kwargs)

    def _muaErr( self, M_params, extra_arg):
        pass

    def _muaWarmup( self, extra_arg):
        return extra_arg
