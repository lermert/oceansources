import numpy as np

class SincBasis(object):
    #from my newer noisi branch

    def __init__(self,K,N,**kwargs):

        self.K = K
        self.N = N
        self.basis_type='sinc_taper'
        self.freq = kwargs['freq']
        self.fmin = kwargs['fmin']
        self.fmax = kwargs['fmax']

    def basis_func(self,k,n):

        """
        Return the sinc function (reference?)
        :type k: int
        :param k: return the k'th taper
        :type n: int
        :param n: Number of samples
        """

        try:
            ix_fmin = np.argmin(np.abs(self.freq - self.fmin))
            ix_fmax = np.argmin(np.abs(self.freq - self.fmax))
            n_supp = ix_fmax - ix_fmin
        except:
            raise ValueError(message='The chosen frequency range is incompatible with the sampling rate.')

        # Number of samples between two neighbouring sinc functions
        d_n = int(round(n_supp / (self.K - 1), 0))
        # How many steps fit onto original frequency axis?
        n_steps = int(round(n / d_n, 0))
        # common x axis for the sinc basis
        x = np.linspace(0, n_steps * np.pi, n)
        # how many steps are below the support interval?
        k_0 = int(round(ix_fmin/d_n,0))
        # now shift x=0 to the desired frequency
        argu = (x / np.pi - k - k_0)
        # non-normalized sinc: divide argument by pi
        # this is necessary in the above line because numpy.sinc
        # returns the normalized sinc, and we want the non-normalized
        norm = np.sqrt((x[-1] - x[0]) / (n * np.pi))
        y = norm * np.sinc(argu)    # normalize
        return(y)