import sys
sys.path.append('src')
import matplotlib.pyplot as plt
import numpy as np
import fit_sed

# get data
lascar = np.genfromtxt('results/lascar_cluster_measurements.cat', usecols=(3, 5, 6, 8, 9, 17, 18), dtype=None)


def assign_alma_wavelengths_by_redshift(z):                                                             
    if z > 1.0:                                                                                       
        return 1300         
    elif z > 0.8:                                                                               
        return 1178.4
    else:                                                 
        return nan  

fig, axes = plt.subplots(4, 5, figsize=(20, 16), sharex=True, sharey=True)                              

for i, [ax, g] in enumerate(zip(axes.flat, lascar)):
    z, S_100, S_100_err, S_160, S_160_err, S_alma, S_alma_err = g
    if np.isnan(S_alma_err):
        fit_sed.find_best_template([[100, S_100, S_100_err], [160, S_160, S_160_err]], z, ax=ax)
    else:
        wave_alma = assign_alma_wavelengths_by_redshift(z)
        fit_sed.find_best_template([[100, S_100, S_100_err], [160, S_160, S_160_err], [wave_alma, S_alma, S_alma_err]], z, ax=ax)
