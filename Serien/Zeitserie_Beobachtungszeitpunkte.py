#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Das Skript liest eine Serie von 1d-Spektren im fits-Format ein
und gibt den Beobachtungszeitpunkt jeden Spektrums in der Konsole aus.
Außerdem wird der Wellenlängenbereich jedes Spektrums in der Konsole ausgedruckt.
Am Ende wird eine ascii-Datei mit den Spektrennamen und dem Beobachtungszeitpunkt
abgespeichert.

Created on Fri Oct 16 17:47:58 2020

@author: lothar schanne
"""

import numpy as np
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import glob

# Fileliste erstellen. Spektren in einem (Unter)Ordner.
files = input("Pfad und Name der Spektren (nutze wildcards) : ")
filelist = glob.glob(files)

# Alphabetisches Sortieren. Bei richtiger Namensgebung ergibt das eine
# zeitliche Ordnung
filelist.sort()

# Ausdruck der Liste
print("\nAnzahl der Spektren: ", len(filelist), "\n")
print("\nSpektrenliste: \n")
print(filelist)


date_obs = list(range(len(filelist)))
jd = list(range(len(filelist)))

print("\nHeadereinträge")

# Beobachtungszeitpunkt einlesen und ausgeben
print("\nSpektrum", "DATE-OBS", 'JD')

for i in range(len(filelist)):
    flux, header = fits.getdata(
        filelist[i], header=True, ignore_missing_end=True)
    if "DATE-OBS" in header:
        date_obs[i] = header["DATE-OBS"]
    elif "DATE_OBS" in header:
        date_obs[i] = header["DATE_OBS"]
    elif "MJD" in header:
        date_obs[i] = header["MJD"]
    elif "JD-MID" in header:
        date_obs[i] = header["JD-MID"]
    else:
        print("Kein Beobachtunsgzeitpunkt im Header")

    if "JD" in header:
        jd[i] = header['JD']
    else:
        print("Kein JD im Header")

    print(filelist[i], date_obs[i], jd[i])

# Spektraler Bereich berechnen und ausgeben
print("    Spektrumdatei    ", "Anfang   ", "Ende")
for i in range(len(filelist)):
    flux, header = fits.getdata(
        filelist[i], header=True, ignore_missing_end=True)

    Anfang = header["CRVAL1"]
    Ende = header["CRVAL1"] + header["NAXIS1"] * header["CDELT1"]
    print(filelist[i], f"{Anfang:.1f}, {Ende:.1f}")

# Abspeichern als ascii-Datei
ascii.write(
    [filelist, date_obs, jd],
    "BeobachtungszeitpunkteUndJD" + ".dat",
    overwrite=True,
    names=["Spektrum", "DATE-OBS", 'JD'],
    format="tab",
)
