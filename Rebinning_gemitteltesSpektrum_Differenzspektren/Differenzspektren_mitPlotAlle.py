#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Differenzspektren.py

Liest tab-Spektren einer Serie mit den Spalten "WAVE" und "FLUX ein
und bildet Differenspektren zum angegebenen Bezugs-Spektrum, die dann als
tab-Datei gespeichert werden. Alle verwendeten Spektren müssen den gleichen
Wellenlängenbereich und die gleiche Schrittweite haben.
Alle Differenzspektren werden abgespeichert und in je einem Plot grafisch
dargestellt, die auch abgespeichert werden.

20231115

@author: lothar
"""
import numpy as np
from astropy.io import ascii
import matplotlib.pyplot as plt
import glob


# Fileliste erstellen. Spektren in einem (Unter)Ordner.
files = input("Pfad und Name der dat-Spektren (nutze wildcards) : ")
filelist = glob.glob(files)

# Alphabetisches Sortieren. Bei richtiger Namensgebung ergibt das eine
# zeitliche Ordnung
filelist.sort()

# Ausdruck der Liste
print("\nSpektrenliste: \n")
print(filelist, end='\n')
print("Anzahl der Spektren: ", len(filelist), "\n")

mean_name = input("Name der dat-Datei für das Bezugs-Spektrum: ")
mean_spec = ascii.read(mean_name, format="tab")
print("Verwendetes Bezugs-Spektrum", mean_name)

for i in range(len(filelist)):
    spec = ascii.read(filelist[i], format="tab")
    diff = np.zeros(len(spec["WAVE"]))
    for j in range(len(spec["WAVE"])):
        diff[j] = spec["FLUX"][j] - mean_spec["FLUX"][j]
    filename = filelist[i].rsplit(".", 1)[0] + "_diffspektrum" + ".dat"
    ascii.write(
        [spec["WAVE"], diff[:]],
        filename,
        overwrite=True,
        names=["WAVE", "FLUX"],
        format="tab",
    )

    if i > 0:
        fig = plt.figure(figsize=(10, 10))
        plt.title(filename, fontsize=10)
        plt.xlabel("Wavelength [Angstroem]")
        plt.ylabel("Differenz Flux")
        plt.ylim(-.3, .45)  # anpassen !!!!!!!!!!!!!!!!!!
        plt.grid(True)
        plt.plot(spec["WAVE"], diff)
        plt.pause(.5)
        # fig.savefig("Differenzspektrum" +
        #             filelist[i].rsplit('.', 1)[0] + ".png", format='png')
        plt.close()

        plt.figure(0)
        plt.plot(spec["WAVE"], diff)

    elif i == 0:
        fig = plt.figure(0, figsize=(10, 10))
        plt.title(filename, fontsize=10)
        plt.xlabel("Wavelength [Angstroem]")
        plt.ylabel("Differenz Flux")
        plt.ylim(-.3, .45)  # anpassen !!!!!!!!!!!!!!!!!!
        plt.grid(True)
        plt.plot(spec["WAVE"], diff)
        plt.pause(.5)
        # fig.savefig("Differenzspektrum" +
        #             filelist[i].rsplit('.')[0] + ".png", format='png')

plt.figure(0)
fig.savefig("Differenzspektren_alle" + ".png", format='png')
