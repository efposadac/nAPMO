import pylibxc
import napmo


class Functional(object):
    """
    Electronic functionals provided by libxc
    ... complete
    """
    def __init__(self, options={}):
        super(Functional, self).__init__()
        self.options = {
            'correlation': 'lda_c_vwn',
            'exchange': 'lda_x',
            'spin': 'unpolarized',
        }

        self.options.update(options)
        self._available = pylibxc.util.xc_available_functional_names()
        self._exchange = self.options.get('exchange')
        self._spin = self.options.get('spin')

        if self.correlation in self.available:
            self._correlation = pylibxc.LibXCFunctional(self.options.get('correlation'), self._spin)
        else:
            raise NotImplementedError(self.correlation+" Functional NOT available!")

        if self.exchange in self.available:
            self._exchange = pylibxc.LibXCFunctional(self.options.get('exchange'), self._spin)
        else:
            raise NotImplementedError(self.exchange+" Functional NOT available!")

    def compute_correlation(self, rho):
        inp = {}
        inp['rho'] = rho
        ret = self._correlation.compute(inp)
        return ret['zk'], ret['vrho']

    def compute_exchange(self, rho):
        inp = {}
        inp['rho'] = rho
        ret = self._exchange.compute(inp)
        return ret['zk'], ret['vrho']

    @property
    def correlation(self):
        return self.options.get('correlation')

    @property
    def exchange(self):
        return self.options.get('exchange')

    @property
    def available(self):
        return self._available

    @property
    def spin(self):
        return self._spin


if __name__ == '__main__':
    import numpy as np
    rho = np.random.random((3))
    func = Functional()
    print(func.available)
    cene, cpot = func.compute_correlation(rho)
    eene, epot = func.compute_exchange(rho)
