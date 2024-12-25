#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Liest eine Zeitserie von 1d-Spektren im fits-Format ein.
Berechnet jeweils das Minimum und Maximum von Wellenlänge und Flux im
Spektrum und schreibt diese Werte in eine ASCII-dat-Datei (Trennzeichen tab).

Stand 20231214
@author: Lothar Schanne
"""

import numpy as np
from astropy.io import fits, ascii
import glob


# Create file list. All spectra in one (sub)folder.
files = input("Path and name of the file (use wildcards) : ")
filelist = glob.glob(files)

# Aalphabetical sorting. If the name is correct, this results in a
# temporal order
filelist.sort()

# Print the list
print("\nList of spectra: \n")
print("Number of spectra: ", len(filelist), "\n")

minimumwave = np.zeros(len(filelist))
maximumwave = np.zeros(len(filelist))
minimumflux = np.zeros(len(filelist))
maximumflux = np.zeros(len(filelist))
JD = np.zeros(len(filelist))

# Read header and flux
for i in range(len(filelist)):
    print(i, filelist[i])
    flux, header = fits.getdata(filelist[i], header=True)
    # Generate the spectrum sections and save them as fit:
    step = header["CDELT1"]
    refpix = header["CRPIX1"]
    # JD der Beobachtung, zur Beschriftung der Spektren verwendet:
    if 'JD' in header:
        JD[i] = header['JD']
    elif 'MJD-OBS' in header:
        JD[i] = header['MJD-OBS']
    elif 'DATE-OBS' in header:
        JD[i] = header['DATE-OBS']
    else:
        print('\nKein Beobachtungsdatum im header')

    wave_erstesPix = header["CRVAL1"] - step * (refpix - 1)

    wave = np.zeros(header["NAXIS1"], dtype=float)
    for k in range(header["NAXIS1"]):
        wave[k] = wave_erstesPix + k * step

    minimumflux[i] = flux.min()
    maximumflux[i] = flux.max()
    minimumwave[i] = wave[flux.argmin()]
    maximumwave[i] = wave[flux.argmax()]


# Abspeichern als ascii-Datei
ascii.write(
    [filelist, JD, minimumwave, minimumflux, maximumwave, maximumflux],
    "MinMax_anJD" + ".dat",
    overwrite=True,
    names=["Spektrum", "JD", 'WellenlangeMinimum', 'FluxMinimum',
           'WellenlaengeMaximum', 'FluxMaximum'],
    format="tab",
)
