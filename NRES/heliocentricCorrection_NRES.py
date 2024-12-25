#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Das Skript liest die baryzentrische Korrektur BC einer Serie von NRES-Spektren 
im Fits-Format ein.
Korrigiert die eingelesenen Spektren um die baryzentrische Korrektur und 
speichert sie als fits mit dem Namenszusatz BCcorrected ab.

@author: Lothar Schanne
Stand 20241211
"""

from PyAstronomy import pyasl
import numpy as np
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table
import glob
from astroplan import FixedTarget


c = 299792458  # Lichtgeschwindigkeit in m/s

# Create file list. Spectra in a (sub)folder.
files = input('Path and name of the spectra (use wildcards) : ')
filelist = glob.glob(files)

# Alphabetical sorting. If the name is correct (e.g. 20180922-xyz.fit),
# this results in a temporal order.
filelist.sort()

# Print the list
print('\nList of spectra: \n')
print('Number of spectra: ', len(filelist), '\n')


# Processing the filelist, reading flux and header:
for i in range(len(filelist)):
    flux, header = fits.getdata(filelist[i], header=True)

    bc_faktor = 1 + header['BARYCORR'] / c

    print()
    print(filelist[i])
    print("Barycentric correction [m/s]: ", header['BARYCORR'])
    print("Alte Anfangswellenlänge CRVAL1: ", header["CRVAL1"])
    print('Alte Stepweite ', header['CDELT1'])

    # Schreiben des RV-korrigierten Spektrums in fits-file
    header["CRVAL1"] = header["CRVAL1"] * bc_faktor
    header['CDELT1'] = header['CDELT1'] * bc_faktor
    newfile = filelist[i].rsplit(".", 1)[0] + "_BCcorrected.fits"
    fits.writeto(newfile, flux, header, overwrite=True,
                 output_verify="ignore")

    print("Neue Anfangswellenlänge CRVAL1: ", header["CRVAL1"])
    print('Neue Stepweite ', header['CDELT1'])
