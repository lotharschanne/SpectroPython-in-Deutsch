#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fileliste erstellen für  wellenlängenkalibrierte 1d_Spektren einer Zeitreihe im
fits-Format. Abspeichern der filelist mit Angaben zum Ende und Beginn der
Wellenlängenskala. Ausdrucken des gemeinsamen Wellenlängenbereichs.

20220125
@author: lothar schanne
"""

import numpy as np
from astropy.io import fits, ascii
import glob

# Fileliste erstellen. Spektren in einem (Unter)Ordner.
files = input("Pfad und Name der Spektren (nutze wildcards) : ")
filelist = glob.glob(files)

# Alphabetisches Sortieren. Bei richtiger Namensgebung ergibt das eine
# zeitliche Ordnung
filelist.sort()

# Ausdruck der Liste
print("\nSpektrenliste: \n")
print(filelist)
print("Anzahl der Spektren: ", len(filelist), "\n")

begin = np.zeros(len(filelist))
end = np.zeros(len(filelist))
print('Spektrum          begin   end')
# Einlesen von flux und header:
for i in range(len(filelist)):
    flux, header = fits.getdata(filelist[i], header=True)
    begin[i] = header["CRVAL1"]
    end[i] = header["CRVAL1"] + header["CDELT1"] * header["NAXIS1"]
    print(filelist[i], begin[i].round(2), end[i].round(2))

ascii.write(
    [filelist, begin.round(2), end.round(2)],
    "Wellenlaengenbereiche.dat",
    names=["Spektrum", "Beginn", "Ende"],
    format="tab",
    overwrite=True,
)

min = 0
max = 10000
for i in range(len(filelist)):
    flux, header = fits.getdata(filelist[i], header=True)
    if header["CRVAL1"] > min:
        min = header["CRVAL1"]
    if header["CRVAL1"] + header["CDELT1"] * header["NAXIS1"] < max:
        max = header["CRVAL1"] + header["CDELT1"] * header["NAXIS1"]
print("\nGemeinsamer Wellenlängenbereich: ", min, max)
