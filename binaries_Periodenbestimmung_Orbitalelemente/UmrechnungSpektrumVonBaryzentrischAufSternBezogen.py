#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Erforderlich:
    Spektren als fits
    Tabelle RVs und Phasen des Sterns berechnet mit den Orbitparametern des Orbits

@author: lothar
"""


import numpy as np
from astropy.io import fits, ascii
import glob
from PyAstronomy import pyasl


# Fileliste erstellen. Spektren in einem (Unter)Ordner.
files = input("Pfad und Name der fits-Differenz-Spektren (nutze wildcards) : ")
filelist = glob.glob(files)

# Alphabetisches Sortieren. Bei richtiger Namensgebung ergibt das eine
# zeitliche Ordnung
filelist.sort()

# Ausdruck der Liste
# print("\nSpektrenliste: \n")
# print(filelist, end='\n')
print("Anzahl der Spektren: ", len(filelist), "\n")

tabelle = ascii.read(
    'theoretischeRVsAlgolABCzuJDZeitpunktenAusSpektren.csv', format='csv')  # anpassen !!!!
print('Spaltennamen der Tabelle:')
print(tabelle.colnames)


for i in range(len(filelist)):
    flux, header = fits.getdata(filelist[i], header=True)
    for j in range(len(tabelle)):
        if tabelle['JD'][j] == header['JD']:
            # falls eine andere RV erwünscht ist anpassen !!!
            rv_A = tabelle['RV_A'][j]
            wave = np.zeros(header["NAXIS1"], dtype=float)
            for k in range(len(wave)):
                if 'CRPIX1' not in header:
                    header['CRPIX1'] = 1
                wave[k] = header["CRVAL1"] + \
                    (k + 1 - header["CRPIX1"]) * header["CDELT1"]
            # Shift the spectrum
            flux_rv, wave_rv = pyasl.dopplerShift(
                wave, flux, -rv_A, edgeHandling="firstlast")
            header["CRVAL1"] = wave[0]
            header['CRPIX1'] = 1
            fits.writeto(filelist[i].rsplit('.fit')[0]+'auf_Stern_bezogen.fit',
                         flux_rv, header, overwrite=True,
                         output_verify="silentfix")  # Dateiname anpassen !!!
            ascii.write([wave, flux_rv],
                        filelist[i].rsplit('.fit')[0]+'auf_Stern_bezogen.csv',
                        overwrite=True,
                        names=['WAVE', 'FLUX'],
                        format="csv")  # Dateiname anpassen !!!
        else:
            pass
