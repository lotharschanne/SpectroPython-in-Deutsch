#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Einlesen einer Serie von 1d-fits-Dateien und smoothen des fluxes mit einer
bestimmten Fensterbreite. Abspeichern der geglätteten fits.

Created on Fri Dec 10 16:10:46 2021

@author: lothar
"""

from PyAstronomy import pyasl
import numpy as np
import matplotlib.pylab as plt
from astropy.io import fits
import glob

plt.switch_backend('Qt5Agg')  # Interaktives GrafikBackend einschalten

# Create file list. All spectra in one (sub)folder.
files = input('Path and name of the file (use wildcards) : ')
filelist = glob.glob(files)

obj = input('Geben Sie den Namen des Objekts ein: ')
namenszusatz = input('Geben Sie den Namenszusatz für die abzuspeichernden Spektren ein: ')

# Aalphabetical sorting. If the name is correct, this results in a
# temporal order
filelist.sort()

# Print the list
print('\nList of spectra: \n')
print('Number of spectra: ', len(filelist), '\n')

# Das smoothing-Fenster muss eine ungerade Zahl sein, anpassen
Fenster = 11

# Read header and flux
for i in range(len(filelist)):
    print(filelist[i], ':')
    flux, header = fits.getdata(filelist[i], header=True)
    flux_smoothed = pyasl.smooth(flux, Fenster, 'flat')
    newfile = obj + '_' + str(header['JD']) + '_' + namenszusatz + '.fits'
    fits.writeto(newfile, flux_smoothed, header,
                 overwrite=True, output_verify='silentfix')
