import numpy as np
import xarray as xr
import pyserf2._misc as _m

class Euv91:
    '''
    EUV91 model class
    '''
    def __init__(self):
        self.bands_dataset = _m.get_euv91_bands_dataset()
        self.lines_dataset = _m.get_euv91_lines_dataset()

        self.chromo_bands_coeffs = np.vstack((np.array(self.bands_dataset['a1'], dtype=np.float64),
                                              np.array(self.bands_dataset['b1'], dtype=np.float64) *
                                              np.array(self.bands_dataset['wchr'], dtype=np.float64))).transpose()

        self.chromo_lines_coeffs = np.vstack((np.array(self.lines_dataset['a1'], dtype=np.float64),
                                              np.array(self.lines_dataset['b1'], dtype=np.float64) *
                                              np.array(self.lines_dataset['wchr'], dtype=np.float64))).transpose()

        self.coronal_bands_coeffs = np.vstack((np.array(self.bands_dataset['a2'], dtype=np.float64),
                                               np.array(self.bands_dataset['b2'], dtype=np.float64) *
                                               np.array(self.bands_dataset['w1'], dtype=np.float64),
                                               np.array(self.bands_dataset['b2'], dtype=np.float64) *
                                               np.array(self.bands_dataset['w2'], dtype=np.float64))).transpose()

        self.coronal_lines_coeffs = np.vstack((np.array(self.lines_dataset['a2'], dtype=np.float64),
                                               np.array(self.lines_dataset['b2'], dtype=np.float64) *
                                               np.array(self.lines_dataset['w1'], dtype=np.float64),
                                               np.array(self.lines_dataset['b2'], dtype=np.float64) *
                                               np.array(self.lines_dataset['w2'], dtype=np.float64))).transpose()

    def get_e(self, f107):

        try:
            if isinstance(f107, float) or isinstance(f107, int):
                return np.array([1., f107], dtype=np.float64).reshape(1, 2)
            return np.vstack([np.array([1., x]) for x in f107], dtype=np.float64)
        except TypeError:
            raise TypeError('Only int, float or array-like object types are allowed.')

    def get_echr(self, lya):
        lya = np.array(lya, dtype=float)

        f_lya = lya * 12400 * 1.602192e-12 / 1215.67
        return np.vstack([[1., flya] for flya in f_lya])


    def get_ecor(self, f107, lya):
        f107 = np.array(f107, dtype=float)
        lya = np.array(lya, dtype=float)

        if f107.size != lya.size:
            raise Exception(f'The number of F10.7 and Lya values does not match. f107 contained {f107.size} '
                            f'elements, Lya contained {lya.size} elements.')

        f107mod = -218.88 + 1.05453*10**-9 * lya
        return np.vstack([[1., fm, f] for fm, f in zip(f107mod, f107)])


    def _check_types(self, f107, lya):
        if not isinstance(f107, (float, int, list, np.ndarray) or not isinstance(lya, (float, int, list, np.ndarray))):
            raise TypeError(f'Only float, int, list and np.ndarray. f107 was {type(f107)}, lya was {type(lya)}')

        if type(f107) != type(lya):
            raise TypeError(f'f107 and f107avg types must be equal. f107 was {type(f107).__name__}, '
                            f'lya was {type(lya).__name__}')
        return True



    def get_spectral_bands(self, *, f107, lya):

        if self._check_types(f107, lya):
            if isinstance(f107, (int, float)):
                f107 = np.array(f107, dtype=float).reshape(1,)
                lya = np.array(lya, dtype=float).reshape(1,)

            echr_data = self.get_echr(lya)
            ecor_data = self.get_ecor(f107, lya)

            echr = np.dot(self.chromo_bands_coeffs, echr_data.T)
            ecor = np.dot(self.coronal_bands_coeffs, ecor_data.T)

            print(echr[5], ecor[5])
            return echr + ecor


    def get_spectral_lines(self, *, f107, lya):

        if self._check_types(f107, lya):
            if isinstance(f107, (int, float)):
                f107 = np.array(f107, dtype=float).reshape(1, )
                lya = np.array(lya, dtype=float).reshape(1, )

            echr_data = self.get_echr(lya)
            ecor_data = self.get_ecor(f107, lya)

            echr = np.dot(self.chromo_lines_coeffs, echr_data.T)
            ecor = np.dot(self.coronal_lines_coeffs, ecor_data.T)
            return echr + ecor

    def get_spectra(self, *, f107, lya):
        return self.get_spectral_bands(f107=f107, lya=lya), self.get_spectral_lines(f107=f107, lya=lya)