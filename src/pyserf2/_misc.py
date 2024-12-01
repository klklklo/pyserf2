import functools
import xarray as xr
from importlib_resources import files


@functools.cache
def _read_coeffs(file):
    return xr.open_dataset(files('pyserf2._coeffs').joinpath(file))


def get_euv91_bands_dataset():
    return _read_coeffs('serf2_bands_dataset.nc').copy()


def get_euv91_lines_dataset():
    return _read_coeffs('serf2_lines_dataset.nc').copy()
