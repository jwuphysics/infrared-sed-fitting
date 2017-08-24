import sys
sys.path.append('src')
import matplotlib.pyplot as plt
import numpy as np
import os
import fit_sed

root_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

if __name__ == '__main__':

    # get data
    data_fname = os.path.join(root_dir, 'results', 'lascar_cluster_measurements.cat')
    lascar = np.genfromtxt(data_fname, usecols=(0, 3, 5, 6, 8, 9, 17, 18), dtype=None)


    def assign_alma_wavelengths_by_redshift(z):                                                             
        if z > 1.0:                                                                                       
            return 1300         
        elif z > 0.8:                                                                               
            return 1178.4
        else:                                                 
            return nan  

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
        
        template, L_IR, norm = fit_sed.find_best_template(measurements, z, ax=ax, verbose=False)

        # aesthetics
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(8, 2e3)
        ax.set_ylim(3e-4, 3e1)
        ax.text(0.05, 0.9, name, transform=ax.transAxes, va='bottom', ha='left', fontsize=12)

        if ax is axes[3, 0]:
            ax.set_xlabel(r'Observed wavelength [$\mu$m]', fontsize=12)
            ax.set_ylabel(r'Flux density [mJy]', fontsize=12)

        print('{} {} {} {}'.format(name, template, L_IR, norm))

    # do some tweaking
    fig.subplots_adjust(wspace=0, hspace=0)

    # save results figure
    result_figname = os.path.join(root_dir, 'results', 'lascar_herschel_sample.pdf')
    fig.savefig(result_figname, dpi=100)