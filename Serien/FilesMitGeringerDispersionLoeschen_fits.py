#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fileliste erstellen für  wellenlängenkalibrierte 1d_Spektren einer Zeitreihe im
fits-Format.  Die Dispersionen (Schrittweiten, CDELT1) der Spektren werden 
überprüft. Wenn die Schrittweite einen in Zeile 40 festgelegten Wert überschreitet 
wird das betreffende Spektrum im Arbeitsverzeichnis gelöscht.

@author: lothar schanne
20241218
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import glob
import os

plt.switch_backend('Qt5Agg')  # Interaktives GrafikBackend einschalten


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
print('Bitte warten. Berechnungen laufen.')

print('Dispersion > 0.2 haben:')
# Einlesen von flux und header:
for i in range(len(filelist)):
    flux, header = fits.getdata(filelist[i], header=True)
    if header['CDELT1'] > .2:
        print(i, filelist[i])
        os.remove(filelist[i])  # Dieser Befehl löscht die Spektrenfiles
