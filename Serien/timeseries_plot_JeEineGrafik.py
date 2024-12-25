#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
timeseries_plot.py

Creates a file list for wavelength calibrated 1d_spectra (fits- or csv-format)
of a time series.

Plots the spectra and save the graph as PNG.

Release 28.3.2021
@author: Lothar Schanne
"""

import numpy as np
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import glob

plt.switch_backend('Qt5Agg')  # Interaktives GrafikBackend einschalten

# Create file list. All spectra in one (sub)folder.
files = input('Path and name of the file (use wildcards) : ')
filelist = glob.glob(files)

# Alphabetical sorting. If the name is correct, this results in a
# temporal order
filelist.sort()

# Print the list
print('\nList of spectra: \n')
print('Number of spectra: ', len(filelist), '\n')


# Read header and flux
for i in range(len(filelist)):
    fig = plt.figure()
    print(filelist[i], ':')

    # für fits-Spektren:
    flux, header = fits.getdata(filelist[i], header=True)
    # Generate the spectrum sections and save them as fit:
    step = header['CDELT1']
    if 'CRPIX1' not in header:
        header['CRPIX1'] = 1
    refpix = header['CRPIX1']
    wave_erstesPix = header['CRVAL1'] - step*(refpix - 1)

    wave = np.zeros(header['NAXIS1'], dtype=float)
    for k in range(header['NAXIS1']):
        wave[k] = wave_erstesPix + k * step

    # für csv-Spectren:
    # spec = ascii.read(filelist[i])
    # print(filelist[i], spec.colnames)
    # wave = spec['WAVE']
    # flux = spec['FLUX']

    plt.plot(wave, flux, 'k-', label=filelist[i], linewidth=.6)
    plt.grid(True)
    plt.title(filelist[i], fontsize=12)  # fontsize anpassen !!!
    plt.xlabel('Wavelength in Angström')
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.ylim(0.8, 2.1)  # anpassen oder auskommentieren !!!!!!!!!!!!!!!!!!!!!
    # plt.ylim(0., 4)  # anpassen oder auskommentieren !!!!!!!!!!!!!!!!!!!!!
    # Legend can be adjusted, e.g. font size (fontsize)
    # plt.legend(bbox_to_anchor=(0., 1.02, 1., .402), loc=3,
    #            ncol=2, mode="expand", borderaxespad=1., fontsize=6)

    fig.savefig(filelist[i].rsplit('.fit')[0] +'.png', format='png')  # anpassen
    plt.pause(.01)
    plt.close()

print('Ende des Programms')
plt.close('all')
