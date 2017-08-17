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
 - implement models other than Kirkpatrick+15 comprehensive library
"""

from astropy.constants import c
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import os
from scipy.interpolate import interp1d
from scipy.optimize import minimize
from scipy.special import erf
import sys

root_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

pacs_filter_fnames = {
    70: os.path.join(root_dir, 'data', 'herschel_filters', 'Herschel-Pacs.blue.dat'),
    100: os.path.join(root_dir, 'data', 'herschel_filters', 'Herschel-Pacs.green.dat'),
    160: os.path.join(root_dir, 'data', 'herschel_filters', 'Herschel-Pacs.red.dat')
}

K15_SED_templates = ['AGN1', 'AGN2', 'AGN3', 'AGN4', 'Composite1', 
                     'Composite2', 'Composite3', 'Composite4', 
                     'SFG1', 'SFG2', 'SFG3']

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
    clean_measurements = np.zeros_like(measurements_array, dtype='float')
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
                pass
            else:
                sys.exit('"{}"" is not supported.'.format(telescope_filter))

        clean_measurements[i, 0] = float(observed_wavelength)
        try:
            clean_measurements[i, 1] = float(flux)
            clean_measurements[i, 2] = float(flux_err)
        except ValueError:
            sys.exit('Flux and uncertainty must be floats.')

    return clean_measurements

def model_sed(template_fname, z):
    """Given a K15 model and redshift, returns the SED in flux density 
    units via two arrays: 
    
    Returns
    -------
        waves : 1d array
            observed wavelengths in microns
        f_nu : 1d array
            observed flux density in mJy
    """
    try:
        template_sed = np.genfromtxt(template_fname, skip_header=4)
    except IOError:
        sys.exit('Something is wrong with the SED template.')
    
    # observed wavelengths [um]
    waves = template_sed[:, 0] 
    waves *= (1 + z)

    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    D_L = cosmo.luminosity_distance(z)
    
    # flux density [mJy]
    L_nu  = template_sed[:, 1]
    f_nu = (L_nu * (u.W / u.Hz) / (4 * np.pi * D_L**2)).to(u.mJy).value
    
    return waves, f_nu

def model_photometry(waves, f_nu, wavelength):
    """The wavelength is [70, 100, 160] for PACS observations. Otherwise 
    return nearest flux density for some other central wavelength.
    """

    if wavelength not in [70, 100, 160]: 
        return f_nu[np.argmin(np.abs(np.log(waves / wavelength)))]

    pacs_filter_fname = pacs_filter_fnames[wavelength]

    filter_waves, transmission = np.genfromtxt(pacs_filter_fname,
                                               unpack=True)
    filter_waves *= 1e-4 # [um]
    
    # interpolate to same wavelengths used in SED
    within_filter = (waves > min(filter_waves)) & (waves < max(filter_waves))
    wave_range = waves[within_filter]
    func_interp = interp1d(filter_waves, transmission, kind='cubic')
    
    interp_transmission = func_interp(wave_range)
    
    flux_density = np.sum([T * f for T, f in \
        zip(interp_transmission, f_nu[within_filter])]) / np.sum(interp_transmission)
        
    return flux_density

def chi_squared(normalization, modeled_fluxes, measured_fluxes, measured_uncertainties):
    """Returns the summed chi^2 for all measurements in one template SED.
    """
    modeled_fluxes *= normalization

    # for detections, get usual negative log-likelihood
    detections = 0.5 * (measured_fluxes - modeled_fluxes)**2 / measured_uncertainties**2

    # for nondetections, use survival analysis likelihood, e.g., Feigelson & Nelson (1985)
    nondetections = 0.5 * (1 + erf(modeled_fluxes / (np.sqrt(2) * measured_uncertainties)))

    return np.sum(np.where(np.isfinite(measured_fluxes), detections, nondetections))

def fit_sed(template, measurements, z, visualize=False):

    # TODO: implement models other than K15 comprehensive models
    assert template in K15_SED_templates

    template_fname = os.path.join(root_dir, 'data', 'kirkpatrick+15',
                   'Comprehensive_library', '{}.txt'.format(template))
    if not os.path.isfile(template_fname):
        sys.exit('Invalid template model entered.')

    # get and unpack redshifted wavelengths and SED
    waves, f_nu = model_sed(template_fname, z)

    # unpack measurements
    modeled_fluxes, measured_fluxes, measured_uncertainties = measurements.T

    # minimize chi-squared, with normalization (first arg) as free parameter
    opt_result = minimize(chi_squared, x0=[1.], 
        args=(modeled_fluxes, measured_fluxes, measured_uncertainties))

    if opt_result['success']:
        chi2 = opt_result['fun']
        norm = opt_result['x'][0]
        print('Template {} successful, with chi^2 = {:.2f}'.format(template, chi2))
    else:
        print('Template {} unsuccessful.'.format(template))

if __name__ == '__main__':

    template = 'AGN2'

    input_measurements = os.path.join(root_dir, 'src', 'test', 'test_measurements_a-1.txt')
    clean_measurements = parse_measurements(input_measurements)

    fit_sed(template=template, measurements=clean_measurements, z=0.87)


    

    pass