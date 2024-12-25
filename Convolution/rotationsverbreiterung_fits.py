#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Einlesen einer Serie von 1d-Spektren im fits-Format.

Nach Eingabe des limb-darkening-Koeffizients und der Rotationsgeschwindigkeit
werden die Spektren rotationsverbreitert.
Die convolvierte fits-Dateien werden abgespeichert, wobei die Dateinamen mit dem
Begriff 'rotverbreitert' ergänzt werden.

Stand 20220105
@author: Lothar Schanne
"""

import matplotlib.pyplot as plt
from astropy.io import fits
import glob
import numpy as np
from PyAstronomy import pyasl

# Fileliste erstellen. Spektren in einem (Unter)Ordner.
files = input("Pfad und Name der Spektren (nutze wildcards) : ")
filelist = glob.glob(files)

# Alphabetisches Sortieren. Bei richtiger Namensgebung ergibt das eine
# zeitliche Ordnung
filelist.sort()

epsilon = float(
    input('Geben Sie den limb-darkening coefficient ein, zwischen 0 und 1: '))
vsini = float(
    input('Geben Sie die projizierte Rotationsgeschwindigkeit in km/s ein: '))

# Ausdruck der Liste
print("\nSpektrenliste: \n")
print("Anzahl der Spektren: ", len(filelist), "\n")

for i in range(len(filelist)):
    flux, header = fits.getdata(filelist[i], header=True)
    step = header['CDELT1']
    if 'CRPIX1' in header:
        refpix = header['CRPIX1']
    else:
        header['CRPIX1'] = 1
        refpix = 1
    wave_erstesPix = header['CRVAL1'] - step*(refpix - 1)

    wave = np.zeros(header['NAXIS1'], dtype=float)
    for k in range(header['NAXIS1']):
        wave[k] = wave_erstesPix + k * step

    # Rotationsverbreiterung
    flux_broadened = pyasl.rotBroad(wave, flux, epsilon, vsini)

    name = filelist[i].rsplit('.', 1)[0]+'_rotverbreitert'+'.fits'

    # Speichern des convolved
    fits.writeto(name, flux_broadened, header,
                 overwrite=True, output_verify='silentfix')
