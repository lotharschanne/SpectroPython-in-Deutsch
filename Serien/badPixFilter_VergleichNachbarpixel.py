#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
badPixFilter.py

Fileliste erstellen für  wellenlängenkalibrierte 1d_Spektren einer Serie im
fit-Format.
Ersatz von einzelnen Pixelwerten, die die Nachbarpixel umd einen Faktor
überschreiten,durch den Mittelwert der Nachbarpixel und abspeichern
aller korrigierten Spektren als sonst unverändertes fit. Name um _
badPixRemoved ergänzt.

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


# Einlesen von flux und header:
for i in range(len(filelist)):
    flux, header = fits.getdata(filelist[i], header=True)
    for j in range(2, len(flux) - 2):
        if (flux[j-2] + flux[j-1]) < flux[j]/1:
            if (flux[j+2]+flux[j+1]) < flux[j]/1:
                flux[j] = (flux[j - 2] + flux[j + 2]) / 2
    # Abspeichern des fit
    filename = filelist[i].rstrip(".fits") + "_badPixRemoved.fits"
    fits.writeto(filename, flux, header, overwrite=True,
                 output_verify="silentfix")
