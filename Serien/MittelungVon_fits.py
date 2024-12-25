#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
1d-Spektren im fits-Format.
Sie müssen alle gleiches header['CRVAL1'], header['CDELT1'], header['NAXIS1'], header['CRPIX1'])
haben.
Diese Variablen werden für jedes Spektrum ausgedruckt (Kontrolle). Mit der
Eingabe von 'y' werden alle Fluxe der Spektrenserie gemittelt und das mittlere 
Spektrum als fit abgespeichert.

Created on Fri Oct 16 17:47:58 2020

@author: lothar
"""

import numpy as np
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
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
print('\nAnzahl der Spektren: ', len(filelist), '\n')


print('Headereinträge')
print('CRVAL1', 'CDELT1', 'NAXIS1', 'CRPIX1')

for i in range(len(filelist)):
    flux, header = fits.getdata(
        filelist[i], header=True, ignore_missing_end=True)

    print(header['CRVAL1'], header['CDELT1'],
          header['NAXIS1'], header['CRPIX1'])

frage = input('Sollen die Spektren gemittelt werden? Dann y eingeben: ')

if frage == 'y':
    fluxsumme = np.zeros(int(header['NAXIS1']))
    for i in range(len(filelist)):
        flux, header = fits.getdata(
            filelist[i], header=True, ignore_missing_end=True)
        fluxsumme += flux
    fluxsumme = fluxsumme / len(filelist)
    fits.writeto('mean_'+files, fluxsumme, header, overwrite=True,output_verify='silentfix')
