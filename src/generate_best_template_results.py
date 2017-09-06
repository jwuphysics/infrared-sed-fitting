"""
John F. Wu

Given some infrared measurements and uncertainties, determine the
best fit teplate SED and uncertainty in the total luminosity.
"""

import fit_sed
import matplotlib.pyplot as plt
import numpy as np
import os

root_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

PLOT_RESULTS = True
VERBOSE = False

def assign_alma_wavelengths_by_redshift(z):
    if z > 1.0:
        return 1300
    elif z > 0.8:
        return 1178.4
    else:
        return np.nan

if __name__ == '__main__':

    # get data
    data_fname = os.path.join(root_dir, 'results', 'lascar_cluster_measurements.cat')
    lascar = np.genfromtxt(data_fname, usecols=(0, 3, 5, 6, 8, 9, 17, 18), dtype=None)

    fig, axes = plt.subplots(4, 5, figsize=(20, 16), sharex=True, sharey=True)

    for ax, g in zip(axes.flat, lascar):
        name, z, S_100, S_100_err, S_160, S_160_err, S_alma, S_alma_err = g

        if np.isnan(S_alma_err):
            measurements = [[100, S_100, S_100_err],
                            [160, S_160, S_160_err]]
        else:
            wave_alma = assign_alma_wavelengths_by_redshift(z)
            measurements = [[100, S_100, S_100_err],
                            [160, S_160, S_160_err],
                            [wave_alma, S_alma, S_alma_err]]

        template, L_IR, norm = fit_sed.find_best_template(measurements, z,
                                NUM_BEST_TEMPLATES=3, visualize=PLOT_RESULTS,
                                ax=ax, verbose=VERBOSE)

        # get the model wavelengths and fluxes
        waves, f_nu = fit_sed.model_sed(template, z)

        # unpack and clean measurements
        measured_waves, measured_fluxes, measured_uncertainties = np.array(measurements).T
        modeled_fluxes = np.array([fit_sed.model_photometry(waves, f_nu, wave) for wave in measured_waves])

        if PLOT_RESULTS:
            # aesthetics
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_xlim(8, 2e3)
            ax.set_ylim(3e-4, 3e1)
            ax.text(0.05, 0.9, name, transform=ax.transAxes, va='bottom', ha='left', fontsize=12)

            if ax is axes[3, 0]:
                ax.set_xlabel(r'Observed wavelength [$\mu$m]', fontsize=12)
                ax.set_ylabel(r'Flux density [mJy]', fontsize=12)

    if PLOT_RESULTS:
        # do some tweaking
        fig.subplots_adjust(wspace=0, hspace=0)

        # save results figure
        result_figname = os.path.join(root_dir, 'results', 'lascar_best_templates.pdf')
        fig.savefig(result_figname, dpi=100)
