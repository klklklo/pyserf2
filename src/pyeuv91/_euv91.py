import numpy as np
import xarray as xr
import pyeuv91._misc as _m

class Euv91:
    '''
    EUV91 model class
    '''
    def __init__(self):
        self.bands_dataset = _m.get_euv91_bands_dataset()
        self.lines_dataset = _m.get_euv91_lines_dataset()

        self.chromo_bands_coeffs = np.vstack((np.array(self.bands_dataset['a1'], dtype=np.float64),
                                              np.array(self.bands_dataset['b1'], dtype=np.float64) *
                                              np.array(self.bands_dataset['wchr'], dtype=np.float64)))

        self.chromo_lines_coeffs = np.vstack((np.array(self.lines_dataset['a1'], dtype=np.float64),
                                              np.array(self.lines_dataset['b1'], dtype=np.float64) *
                                              np.array(self.lines_dataset['wchr'], dtype=np.float64)))

        self.coronal_bands_coeffs = np.vstack((np.array(self.bands_dataset['a2'], dtype=np.float64),
                                               np.array(self.bands_dataset['b2'], dtype=np.float64) *
                                               np.array(self.bands_dataset['w1'], dtype=np.float64),
                                               np.array(self.bands_dataset['b2'], dtype=np.float64) *
                                               np.array(self.bands_dataset['w2'], dtype=np.float64)))

        self.coronal_lines_coeffs = np.vstack((np.array(self.lines_dataset['a2'], dtype=np.float64),
                                               np.array(self.lines_dataset['b2'], dtype=np.float64) *
                                               np.array(self.lines_dataset['w1'], dtype=np.float64),
                                               np.array(self.lines_dataset['b2'], dtype=np.float64) *
                                               np.array(self.lines_dataset['w2'], dtype=np.float64)))

    def get_e(self, f107):

        try:
            if isinstance(f107, float) or isinstance(f107, int):
                return np.array([1., f107], dtype=np.float64).reshape(1, 2)
            return np.vstack([np.array([1., x]) for x in f107], dtype=np.float64)
        except TypeError:
            raise TypeError('Only int, float or array-like object types are allowed.')

    def get_f107mod(self, f107):
        return -218.88 + 1.05453e-9 * f107

    def get_spectral_bands(self, *, f107, lya):
        pass

    def get_spectral_lines(self, *, f107, lya):
        pass

    def get_spectra(self, *, f107, lya):
        return self.get_spectral_bands(f107=f107, lya=lya), self.get_spectral_lines(f107=f107, lya=lya)