"""
John F. Wu
==========

Fitting infrared SEDs

Assume that the measured fluxes will be given in one of three ways:

(1) array-like with rows of [observed_wavelength, flux, flux_err]
(2) array-like with rows of ["Telescope filter", flux, flux_err]
(3) a valid filename that points to a space-delimited catalog with rows
    as specified above

The third item in each row is either the measured uncertainty on the
flux, or in cases where the second item is nan or less than 0, the third
item is treated as the root mean square noise.
    
If the wavelength unit is unspecified, assume microns. If flux density 
unit is unspecified, assume mJy.

TODO list:
==========
 - parse SPIRE measurements
"""

from astrop.constants import c
from astropy.cosmology import FlatLambdaCDM
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import os
import sys

def parse_measurements(measurements):
    """Switch method for different types of flux measurement input. 
    Returns 2d array of floats. Note that first column can be an exact
    observed central wavelength, or a number in [70, 100, 160, ...], 
    which serves as a proxy for pacs-blue, or pacs-green, etc.
    """
    if isinstance(measurements, basestring):
        try:
            measurements_array = np.genfromtxt(measurements)
        except IOError:
            sys.exit('String is not valid file of measurements.')
    elif isinstance(measurements, (list, tuple, np.ndarray)):
        measurements_array = np.array(measurements)
    else:
        sys.exit('There is a problem with the measurements.')

    # ensure that the array has 2 dimensions
    if measurements_array.shape == (3,):
        measurements_array = measurements_array.reshape(1, 3)

    # parse each row at a time and store in clean_measurements array
    clean_measurements = np.zeros_like(measurements_array, dtype='f8')
    for i, row in enumerate(measurements_array):
        try:
            observed_wavelength, flux, flux_err = row
        except:
            sys.exit('Each row must have three elements.')
        
        try:
            observed_wavelength = float(observed_wavelength)
        except ValueError:
            telescope_filter = observed_wavelength.lower()
            if 'pacs' in telescope_filter:
                if '70' in telescope_filter or 'blue' in telescope_filter:
                    observed_wavelength = 70
                elif '100' in telescope_filter or 'green' in telescope_filter:
                    observed_wavelength = 100
                elif '160' in telescope_filter or 'red' in telescope_filter:
                    observed_wavelength = 160
                else:
                    sys.exit('Incorrect PACS filter entered.')
            elif 'spire' in telescope_filter:
                # TODO later
            else:
                sys.exit('"{}"" is not supported.'.format(telescope_filter))

        clean_measurements[i, 0] = float(observed_wavelength)
        try:
            clean_measurements[i, 1] = float(flux)
            clean_measurements[i, 2] = float(flux_err)
        except ValueError:
            sys.exit('Flux and uncertainty must be floats.')

    return clean_measurements



def fit_sed(model, measurements):
    pass

if __name__ == '__main__':
    pass