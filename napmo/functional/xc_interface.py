import numpy as np

import pylibxc
import napmo


class Functional(object):
    """
    Electronic functionals provided by libxc
    ... complete
    """
    def __init__(self, symbol, options={}):
        super(Functional, self).__init__()

        _database = {
            'lda': {'correlation': 'lda_c_vwn', 'exchange': 'lda_x'}
        }

        self.options = {
            'correlation': 'lda_c_vwn',
            'exchange': 'lda_x',
            'spin': 'unpolarized',
        }

        self.options.update(options)
        self.options.update(_database[options['functional']['e-']])

        self._symbol = symbol
        self._x_factor = 0.0

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

        self.show()

    def compute_correlation(self, rho):
        inp = {}
        inp['rho'] = rho
        ret = self._correlation.compute(inp)

        if self.spin == 'polarized':
            ret['vrho'] = self._reorder_vrho(ret['vrho'])

        return ret['zk'], ret['vrho']

    def compute_exchange(self, rho):
        inp = {}
        inp['rho'] = rho
        ret = self._exchange.compute(inp)

        if self.spin == 'polarized':
            ret['vrho'] = self._reorder_vrho(ret['vrho'])

        return ret['zk'], ret['vrho']

    def show(self):
        """
        Prints information of the object.
        """
        print("\nFunctional Information:", self.symbol)
        print("-"*(24+len(self.symbol)))
        print(self.options)

    def _reorder_vrho(self, vrho):
        aux = np.zeros(vrho.shape)

        aux[0, :] = np.hstack([vrho[0, ::2], vrho[1, ::2]])
        aux[1, :] = np.hstack([vrho[0, 1::2], vrho[1, 1::2]])

        return aux

    @property
    def symbol(self):
        return self._symbol

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

    @property
    def x_factor(self):
        return self._x_factor
