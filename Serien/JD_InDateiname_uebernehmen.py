#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Das Skript schreibt das Beobachtungsdatum als JD an den Anfang der Dateinamen
einer Serie von Spektren und speichert die fits mit dem neuen Namen als fits ab.

20240201
@author: lothar schanne
"""
from astropy.io import ascii, fits
import glob


# Fileliste erstellen. Spektren in einem (Unter)Ordner.
files = input("Pfad und Name der Spektren (nutze wildcards) : ")
filelist = glob.glob(files)

# Alphabetisches Sortieren. Bei richtiger Namensgebung ergibt das eine
# zeitliche Ordnung
filelist.sort()

# Ausdruck der Liste
# print("\nSpektrenliste: \n")
# print(filelist, end='\n')
print("Anzahl der Spektren: ", len(filelist), "\n")


for i in range(len(filelist)):
    flux, header = fits.getdata(filelist[i], header=True)
    JD = header['JD']
    print('Dateiname: ', filelist[i])
    name = 'JD_'+str(JD)+'_'+filelist[i]
    print('neuer Dateiname: ', name)
    fits.writeto(name, flux, header, overwrite=True)
