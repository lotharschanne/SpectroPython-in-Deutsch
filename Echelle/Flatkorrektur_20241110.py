#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Die Ordnungen eines Echellespektrums im fits-Format werden eingelesen (Objekt 
und Flats).
Überprüft, ob die Wellenlänge des ersten Pixels, Schrittweite und Anzahl der Pixel
für Flat und Objekt übereinstimmen und wenn ja, teilt die Objektordnung durch 
das zugehörige Flat.
Es kann gewählt werden, ob einfach nur das Objektspektrum durch das Flatspektrum 
geteilt wird oder alternativ dieses flatkorrigierte Spektrum wieder mit einem
geglätteten Flat multipliziert wird (und deshalb wieder die Form der 
Blaze-Funktion bekommt).

20241111
@author: Lothar Schanne
"""

import numpy as np
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import glob
from scipy.interpolate import make_smoothing_spline


def wavecalc(flux, header):
    """Berechnet die Wellenlängen aus den Headerdaten der fits-Datei."""
    wave = np.zeros(header["NAXIS1"])
    if "CRPIX1" not in header:
        header["CRPIX1"] = 1
    header["CRVAL1"] = header["CRVAL1"] + \
        (1 - header["CRPIX1"]) * header["CDELT1"]
    for i in np.arange(header["NAXIS1"]):
        wave[i] = header["CRVAL1"] + i * header["CDELT1"]
    return wave, flux


plt.switch_backend('Qt5Agg')  # Interaktives GrafikBackend einschalten

#  **************  Teilung durch das Flat: ***********************
# Create file list for objects. All spectra in one (sub)folder.
files = input('Path and name of the OBJECT files (use wildcards) : ')
filelist_obj = glob.glob(files)
filelist_obj.sort()

# Create file list for flats All spectra in one (sub)folder.
files = input('Path and name of the FLAT files (use wildcards) : ')
filelist_flat = glob.glob(files)
# Alphabetical sorting. If the name is correct, this results in a temporal order
filelist_flat.sort()

frage = input('Möchten Sie Objektspektrum durch Flatspektrum teilen (1) oder \
zusätzlich nach diesem teilen wieder mit dem geglätteten Flat multiplizieren (2)? \
Geben Sie 1 oder 2 ein: ')

# Print the list
print('Number of spectra: ', len(filelist_obj), '\n')
print('Number of flats: ', len(filelist_flat), '\n')

# Read headers and fluxes
for i in range(len(filelist_obj)):
    print()
    print(filelist_obj[i], ':')

    # für Object fits-Spektren:
    flux_obj, header_obj = fits.getdata(filelist_obj[i], header=True)

    wave_obj, flux_obj = wavecalc(flux_obj, header_obj)

    print(f'Objekt Wellenlänge des ersten Pixels: {wave_obj[0]:.2f}')
    print('Anzahl der Pixel: ', header_obj['NAXIS1'])
    print('Schrittweite', header_obj['CDELT1'])

    # für Flat fits-Spektren:
    flux_flat, header_flat = fits.getdata(filelist_flat[i], header=True)
    # Vermeidung von Nullen im Flux
    flux_mask = flux_flat == 0
    flux_flat[flux_mask] = 1e-10

    wave_flat, flux_flat = wavecalc(flux_flat, header_flat)

    smoothing = 100  # Rigizität des Splines, evtl. anpassen !!!!!!!!!!!!!!!!
    flatfunc = make_smoothing_spline(
        wave_flat, flux_flat, lam=smoothing)

    print(f'Flat Wellenlänge des ersten Pixels: {wave_flat[0]:.2f}')
    print(f'Flat Wellenlänge des letzten Pixels: {wave_flat[-1]:.2f}')
    print('Anzahl der Pixel: ', header_flat['NAXIS1'])
    print('Schrittweite', header_flat['CDELT1'])

    if (wave_flat[0] == wave_obj[0]
        and header_flat['CDELT1'] == header_obj['CDELT1']
            and header_flat['NAXIS1'] == header_obj['NAXIS1']):

        # Korrektur durch flat
        if frage == '1':
            flux_flatkorrigiert = flux_obj / flux_flat
            flux_flatkorrigiert = flux_flatkorrigiert / flux_flatkorrigiert.mean()

        if frage == '2':
            flux_flatkorrigiert = (flux_obj / flux_flat) * flatfunc(wave_obj)
            flux_flatkorrigiert = flux_flatkorrigiert / flux_flatkorrigiert.mean()

        fits.writeto(filelist_obj[i].rsplit('.')[
            0] + '_flatkorrigiert.fits', flux_flatkorrigiert, header_obj,
            overwrite=True, output_verify='silentfix')
