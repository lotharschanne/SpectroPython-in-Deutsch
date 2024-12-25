#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Umwandeln einer Serie von wellenlängenkalibrierten 1d-Spektren im fits-Format in
Textformat .csv (Komma-separiert) und .dat-Format (tab-separiert).
Mit Spaltenüberschriften 'WAVE' und 'FLUX'.

Stand 20180815
@author: lothar schanne
"""

import numpy as np
from astropy.io import fits
from astropy.io import ascii
import glob

# Fileliste erstellen. Spektren in einem (Unter)Ordner.
files = input('Pfad und Name der Spektren (nutze wildcards) : ')
filelist = glob.glob(files)

# Alphabetisches Sortieren. Bei richtiger Namensgebung ergibt das eine
# zeitliche Ordnung
filelist.sort()

# Ausdruck der Liste
print('\nSpektrenliste: \n')
print(filelist)
print('Anzahl der Spektren: ', len(filelist), '\n')
print('Bitte warten. Berechnungen laufen.')

for i in range(len(filelist)):
    #   Einlesen von Header und Daten
    flux, header = fits.getdata(filelist[i], header=True)

    #   Prüfung auf die nötigen header-Einträge
    print('\nSpektrum: ', filelist[i])
    print('Ausgabe der zur Wellenlängenberechnung nötigen Headereinträge:')
    if 'NAXIS' in header:
        print('Dimension, NAXIS:                        ', header['NAXIS'])
    else:
        print('Das ist kein 1d-Spektrum !')
    if 'NAXIS1' in header:
        nax = header['NAXIS1']
        print('Anzahl der Werte (Abszisse), NAXIS1:     ', nax)
    else:
        print('NAXIS1 fehlt im header !')
    if 'CRVAL1' in header:
        crval = header['CRVAL1']
        print('Anfangs-Wellenlänge, CRVAL1:             ', crval)
    else:
        print('CRVAL1 fehlt im header !')
    if 'CDELT1' in header:
        cdel = header['CDELT1']
        print('Schrittweite der Wellenlänge, CDELT1:    ', cdel)
    else:
        print('CDELT1 fehlt im header !')
    if 'CRPIX1' not in header:
        header['CRPIX1'] = 1

    #   Erzeugen von numpy-Arrays mit den Wellenlängen und Fluxes des Spektrums
    wave = np.ones(nax, dtype=float)
    for j in range(nax):
        wave[j] = crval + (j - header['CRPIX1'] + 1) * cdel
    ascii.write([wave, flux], filelist[i].split('.')[0]+'.csv', overwrite=True,
                names=['WAVE', 'FLUX'], format='csv')
    ascii.write([wave, flux], filelist[i].split('.')[0]+'.dat', overwrite=True,
                names=['WAVE', 'FLUX'], format='tab')

print('Ende des Programms')
