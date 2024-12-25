#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
badPixFilter.py

Fileliste erstellen für  wellenlängenkalibrierte 1d_Spektren einer Serie im
fit-Format.
Ersatz von einzelnen Pixelwerten, die oberhalb eines Grenzwertes (badpixel) sind,
durch den Mittelwert der Nachbarpixel und abspeichern aller korrigierten
Spektren als sonst unverändertes fit. Name um _badPixRemoved ergänzt.

Das Skript eignet sich nur für auf das Kontinuum normierte Spektren.

Created on Fri Apr  9 18:09:23 2021

@author: lothar
"""
import numpy as np
from astropy.io import fits
import glob

# Fileliste erstellen.
files = input("Pfad und Name der Spektren (nutze wildcards) : ")
filelist = glob.glob(files)

# Ausdruck der Liste
print("\nSpektrenliste: \n")
print(filelist)
print("Anzahl der Spektren: ", len(filelist), "\n")

grenzwert = float(
    input(
        "Geben Sie den Grenzwert ein, oberhalb dessen die Fluxwerte durch die der Nachbarn ersetzt werden: "
    )
)

# Einlesen von flux und header:
for i in range(len(filelist)):
    flux, header = fits.getdata(filelist[i], header=True)
    for j in range(2, len(flux) - 2):
        if flux[j] > grenzwert:
            flux[j] = (flux[j - 2] + flux[j - 1] +
                       flux[j + 1] + flux[j + 2]) / 4
    # Abspeichern des fit
    filename = filelist[i].rstrip(".fit") + "_badPixRemoved.fit"
    fits.writeto(filename, flux, header, overwrite=True,
                 output_verify="silentfix")
