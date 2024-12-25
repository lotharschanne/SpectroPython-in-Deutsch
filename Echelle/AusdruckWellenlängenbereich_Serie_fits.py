#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 15:38:10 2024

@author: lothar
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import glob


plt.switch_backend('Qt5Agg')  # Interaktives GrafikBackend einschalten

# Create file list for objects. All spectra in one (sub)folder.
files = input('Pfad und Name der normierten Ordnungen (use wildcards) : ')
filelist_obj = glob.glob(files)
filelist_obj.sort()

# Print the list
print('\nList of spectra: \n')
print('Number of spectra: ', len(filelist_obj), '\n')


for i in range(len(filelist_obj)):
    flux, header_obj = fits.getdata(filelist_obj[i], header=True)
    # Generate the spectrum sections and save them as fit:
    step_obj = header_obj['CDELT1']
    if 'CRPIX1' not in header_obj:
        header_obj['CRPIX1'] = 1
    refpix = header_obj['CRPIX1']
    wave_erstesPix_obj = header_obj['CRVAL1'] - step_obj*(refpix - 1)
    wave_obj = np.zeros(header_obj['NAXIS1'], dtype=float)
    for k in range(header_obj['NAXIS1']):
        wave_obj[k] = wave_erstesPix_obj + k * step_obj
    print('Spektrum: ', filelist_obj[i])
    print(f'Wellenlänge des ersten Pixels: {wave_erstesPix_obj:.2f}')
    print(f'Wellenlänge des letzten Pixels:  {wave_obj[-1]:.2f}')
    print('Anzahl der Pixel: ', header_obj['NAXIS1'])
    print('Schrittweite', step_obj)
    print()
