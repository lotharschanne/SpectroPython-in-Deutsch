#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Berechnet für eine auf das Kontinuum normierte Zeitserie im fits-Format
die Äquivalentweite einer in Bereiche aufgeteilten Linie.
Eingabe der für alle Spektren der Serie geltenden Integrationsgrenzen 
(in Angström) über eine Liste, die in Zeile 25 an den jeweiligen Fall 
anzupassen ist.
Die EW-Berechnung setzt voraus, dass für alle Spektren der Serie die gleichen
Wellenlängenintervalle für die Berechnung der Integrale verwendet werden können,
also keine wesentlichen RV-Änderungen stattfinden. Deshalb grundsätzlich 
baryzentrisch korrigierte Spektren benutzen.
Die ermittelten EW's werden in einer asci-Datei (Format csv) gespeichert.

Created 2023-04-26

@author: lothar schanne
"""

import numpy as np
from astropy.io import fits
import glob

###############################################################################
# Liste der Wellenlängen, die als Bereichsgrenzen für die EW-Integration dienen
WL = [6500.1, 6501.2, 6502.3]  # bitte anpassen
##############################################################################

# Fileliste erstellen. Spektren in einem (Unter)Ordner.
files = input("Pfad und Name der Spektren (nutze wildcards) : ")
filelist = glob.glob(files)

# Alphabetisches Sortieren. Bei richtiger Namensgebung ergibt das eine
# zeitliche Ordnung
filelist.sort()

# Ausdruck der Liste
# print("\nSpektrenliste: \n")
# print(filelist)
print("Anzahl der Spektren: ", len(filelist), "\n")

fileobj = open('EWs.csv', 'w')
fileobj.write('Liste der Bereichsgrenzen in Angström:,')
fileobj.write(str(WL)+'\n')
fileobj.write('Spektrumname,')
fileobj.write('JD,')
fileobj.write('Äquivalentweiten in Angström\n')


# Abarbeiten der filelist
for k in np.arange(len(filelist)):
    file = filelist[k]
    sp = fits.open(filelist[k], ignore_missing_end=True)
    try:
        JD = sp[0].header['JD']
    except:
        JD = None

    fileobj.write(file)
    fileobj.write(',')
    fileobj.write(str(JD))
    fileobj.write(',')
    # Generation of arrays with the wavelengths and fluxes of the spectrum
    flux = np.array(sp[0].data)
    wave = np.ones(sp[0].header["NAXIS1"], dtype=float)

    for i in np.arange(sp[0].header["NAXIS1"]):
        wave[i] = (
            sp[0].header["CRVAL1"]
            + (i - sp[0].header["CRPIX1"] + 1) * sp[0].header["CDELT1"]
        )

    # Close the fits-file:
    sp.close()

    # Abarbeiten der Integrationsbereiche:
    for m in np.arange(len(WL)-1):
        EW = 0
        for n in np.arange(len(wave)):
            if wave[n] >= WL[m] and wave[n] <= WL[m+1]:
                EW += (1 - flux[n])*sp[0].header['CDELT1']
        fileobj.write(str(EW))
        fileobj.write(',')
    fileobj.write('\n')

fileobj.close()
