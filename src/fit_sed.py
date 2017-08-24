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
 - report integrated IR luminosity
 - write more comprehensive tests
 - parse SPIRE measurements
 - implement models other than Kirkpatrick+15 comprehensive library
"""

from astropy.constants import c as SPEED_OF_LIGHT
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

VISUALIZE = True

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
            measurements_array = np.genfromtxt(measurements, dtype=None)
        except IOError:
            sys.exit('String is not valid file of measurements.')
    elif isinstance(measurements, (list, tuple, np.ndarray)):
        measurements_array = np.array(measurements)
    else:
        sys.exit('There is a problem with the measurements.')

    # ensure that the array has 2 dimensions
    if len(measurements_array) is 1 or type(measurements_array) is np.void:
        measurements_array = np.array(measurements_array)

    # parse each row at a time and store in clean_measurements array
    clean_measurements = np.zeros((len(measurements_array), 3), 
                                  dtype='float')

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

def read_K15_template(template):
    """Reads in a K15 template SED, returning an array of 
    wavelength and corresponding array of specific luminosity."""

    template_fname = os.path.join(root_dir, 'data', 'kirkpatrick+15',
                   'Comprehensive_library', '{}.txt'.format(template))
    if not os.path.isfile(template_fname):
        sys.exit('Invalid template model entered.')

    try:
        template_sed = np.genfromtxt(template_fname, skip_header=4)
    except IOError:
        sys.exit('Something is wrong with the SED template.')
    
    # rest wavelengths, luminosity [um]
    waves = template_sed[:, 0] 
    L_nu  = template_sed[:, 1]

    return waves, L_nu


def model_sed(template, z):
    """Given a K15 model and redshift, returns the SED in flux density 
    units via two arrays: 
    
    Returns
    -------
        waves : 1d array
            observed wavelengths in microns
        f_nu : 1d array
            observed flux density in mJy
    """
    waves, L_nu = read_K15_template(template)

    # observed wavelengths
    waves *= (1 + z)

    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    D_L = cosmo.luminosity_distance(z)
    
    # flux density [mJy]
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

def chi_squared(normalization, model, data, data_err):
    """Returns the summed chi^2 for all measurements in one template SED.
    """
    model = np.array(model) * normalization

    # for detections, get usual negative log-likelihood
    detections = 0.5 * (data - model)**2 / data_err**2

    # for nondetections, use survival analysis likelihood, e.g., 
    # Feigelson & Nelson (1985)
    #nondetections = 0.5 * (1 + erf(model / (np.sqrt(2) * data_err)))

    # UPDATE: we'll instead use a normal distribution centered on 
    # zero since it gives the intuitively correct answer when fitting. 
    nondetections = 0.5 * (0 - model)**2 / data_err**2

    # For a better treatment, we may want to use the survival function
    #  for a non-normal distribution... i.e., if we know the  
    # completeness as a function of flux, we may be able to create a 
    # log-likelihood of some non-detection for any model prediction.
    

    return np.sum(np.where(np.isfinite(data), detections, nondetections))

def fit_sed(template, measurements, z, verbose=True):

    assert template in K15_SED_templates

    # get and unpack redshifted wavelengths and SED
    waves, f_nu = model_sed(template, z)

    # unpack measurements and then model what they should be
    measured_waves, measured_fluxes, measured_uncertainties = measurements.T

    modeled_fluxes = np.array([model_photometry(waves, f_nu, wave) for wave in measured_waves])

    # minimize chi-squared, with normalization (first arg) as free parameter
    opt_result = minimize(chi_squared, x0=[1.], 
        args=(modeled_fluxes, measured_fluxes, measured_uncertainties))

    if opt_result['success']:
        chi2 = opt_result['fun']
        norm = opt_result['x'][0]
        if verbose:
            print('Template {} successful, with chi^2 = {:.2f}'.format(template, chi2))

        return chi2, norm
    else:
        if verbose:
            print('Template {} unsuccessful.'.format(template))

        return np.nan, np.nan

def infrared_luminosity(template, norm=1.):
    """Calculates the IR luminosity from 8 to 1000 microns, assuming some 
    normalizing factor, norm (default 1). Returns log_10(L_IR/L_sun)."""

    waves, L_nu = read_K15_template(template)

    wavelength_range = np.logical_and((waves >= 8), (waves <= 1000))

    # integrate L_nu over dnu
    freqs = SPEED_OF_LIGHT.to(u.micron / u.s).value / waves[wavelength_range] # Hz
    delta_freqs= freqs[:-1] - freqs[1:]
    L_IR = np.sum(dnu* l_nu for dnu, l_nu in zip(delta_freqs, L_nu[wavelength_range])) * norm

    return L_IR * (u.W).to(u.Lsun)

def find_best_template(input_measurements, z, library=K15_SED_templates, 
                       visualize=True, ax=None, verbose=True):
    """Executes the entire pipeline, attempting to fit all templates to the 
    measured data. 

    If visualize is True and a matplotlib axis object is provided, then the
    best fits will be plotted on the given axis.
    """
    clean_measurements = parse_measurements(input_measurements)

    # record chi^2 of successful templates
    chi_squareds = np.zeros_like(library, dtype=float) * np.nan
    normalizations = np.zeros_like(library, dtype=float) * np.nan
    infrared_luminosities = np.zeros_like(library, dtype=float) * np.nan

    for i, template in enumerate(library):
        chi2, norm = fit_sed(template=template, measurements=clean_measurements, z=z, 
                             verbose=verbose)

        chi_squareds[i] = chi2
        normalizations[i] = norm
        infrared_luminosities[i] = infrared_luminosity(template, norm)

    model_order = np.argsort(chi_squareds)
    lowest_chi2 = np.nanmin(chi_squareds)

    best_template = library[np.nanargmin(chi_squareds)]
    best_L_IR = infrared_luminosities[np.nanargmin(chi_squareds)]
    best_norm = normalizations[np.nanargmin(chi_squareds)]

    if VISUALIZE:
        # plot up to the top 5 successful templates
        if ax is None:
            fig, ax = plt.subplots(1, 1, figsize=(8, 5))

        NUM_BEST_TEMPLATES = 3
        for arg, ls in zip(model_order[:NUM_BEST_TEMPLATES], ['-', '--', ':', '-.']):

            template = K15_SED_templates[arg]
            chi2 = chi_squareds[arg]
            norm = normalizations[arg]
            L_IR = infrared_luminosity(best_template, norm)

            waves, f_nu = model_sed(template, z)
            
            if np.isnan(chi2):
                continue

            log_L_IR_text = r'$\log (L_{{\rm IR}}/L_\odot)={:.2f}$'.format(np.log10(L_IR))

            ax.plot(waves, f_nu * norm, alpha=(0.5 + lowest_chi2 / (2 * chi2)), ls=ls, c='k',
                label='{} ($\chi^2={:.2f}$, {:s})'.format(template[:3]+template[-1], 
                                                          chi2, log_L_IR_text))

        # plot measurements (including upper limits)
        for measured_wave, measured_flux, measured_uncertainty in clean_measurements:
            if np.isfinite(measured_flux):
                ax.errorbar(measured_wave, measured_flux, measured_uncertainty, 
                    marker='o', color='black', ls='', zorder=10)
            else:
                ax.errorbar(measured_wave, 3 * measured_uncertainty, yerr=measured_uncertainty,
                            color='black', uplims=True)
        
        # aesthetics (if within a subplot, then don't show these)
        ax.set_xscale('log')
        ax.set_yscale('log')

        if ax is None:
            ax.set_xlabel(r'Observed wavelength [$\mu$m]', fontsize=12)
            ax.set_ylabel(r'Flux density [mJy]', fontsize=12)

        # legend
        ax.legend(frameon=False, loc='lower center')


    return best_template, best_L_IR, best_norm

if __name__ == '__main__':

    z = 0.87
    input_measurements = os.path.join(root_dir, 'src', 'test', 'test_measurements_a-3.txt')

    best_template, best_L_IR, best_norm = find_best_template(input_measurements, z, visualize=True)