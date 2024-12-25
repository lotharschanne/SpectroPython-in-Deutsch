#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Das Skript schreibt das Beobachtungsdatum als JD und die Anfangswellenlänge 
an den Anfang der Dateinamen einer Serie von Spektren und speichert mit dem 
neuen Namen als fits ab.

20241208
@author: lothar schanne
"""
from astropy.io import ascii, fits
import glob


def anfangswellenlaenge_calc(header):
    """
    Berechnet die Wellenlänge des ersten Pixels aus den Headerdaten der fits-Datei
    """
    if "CRPIX1" not in header:
        header["CRPIX1"] = 1.
    header['CRVAL1'] = header["CRVAL1"] + \
        (1 - header["CRPIX1"]) * header["CDELT1"]
    header["CRPIX1"] = 1.
    return header


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
    JD = '{:.2f}'.format(header['MJD-OBS'])
    header = anfangswellenlaenge_calc(header)

    print('Dateiname: ', filelist[i])
    name = 'JD_'+str(JD)+'_Anfangswellenlänge_' + \
        str(int(header['CRVAL1']))+'.fits'
    print('neuer Dateiname: ', name)
    fits.writeto(name, flux, header, overwrite=True)
