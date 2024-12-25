#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Einlesen einer Serie von heliozentrisch korrigierten Spektren im fits-Format,
Darstellung einer wählbaren Linie im Geschwindigkeitsraum in zwei Plots:
    1. alle Spektren übereinander geplottet.
    2. alle Spektren mit einem wählbaren offset übereinander versetzt geplottet.
Die plots können als pdf abgespeichert werden.

20230216
lothar schanne
"""

import numpy as np
from astropy.io import fits
import glob
import matplotlib.pylab as plt

# lokaler Modul, muss im gleichen Verzeichnis stehen wie das Skript (oder das
# Verzeichnis des Moduls im Pythonpath, damit Python ihn findet)
import Linienlisten


plt.switch_backend('Qt5Agg')  # Interaktives GrafikBackend einschalten


# Fileliste erstellen. Spektren in einem (Unter)Ordner.
files = input("Pfad und Name der Spektren (nutze wildcards) : ")
filelist = glob.glob(files)

# Alphabetisches Sortieren. Bei richtiger Namensgebung ergibt das eine
# zeitliche Ordnung
filelist.sort()

# Ausdruck der Liste
print("\nSpektrenliste: \n")
print("Anzahl der Spektren: ", len(filelist), "\n")


# Auswahl der Linie
linie, wellenlaenge = Linienlisten.linienauswahl()

bereich = float(input('Geben Sie den darzustellenden Geschwindigkeitsbereich um \
die gewählte Linie in km/s ein: '))

# Eingabe der Systemgeschwindigkeit:
systemgeschwindigkeit = float(
    input("Geben Sie die Systemgeschwindigkeit in km/s ein: "))
systemgeschwindigkeit = systemgeschwindigkeit / 299710 * wellenlaenge

# Parameter für den Abstand zwischen den Spektren im offset Plot fig1:
offset = float(input("Please enter the desired offset: "))
obj = input("Please enter the object name: ")


# Abarbeiten der Spektrenliste:
for i in range(len(filelist)):
    flux, header = fits.getdata(filelist[i], header=True)
    if "CRPIX1" in header:
        refpix = int(header["CRPIX1"])
    else:
        refpix = 1
    step = float(header["CDELT1"])
    lambda0 = float(header["CRVAL1"]) - (refpix - 1) * step
    JD = header['JD']

    wave = np.zeros(header['NAXIS1'])
    for j in range(header['NAXIS1']):
        wave[j] = lambda0 + j * step

    index_wellenlaenge = int(
        (wellenlaenge - lambda0 + systemgeschwindigkeit) / step)

    vstep = step / wellenlaenge * 299710
    pixbereich = int(bereich / vstep)

    fluxbereich = flux[(
        index_wellenlaenge - pixbereich):(index_wellenlaenge + pixbereich)]
    wellenlaengenbereich = wave[(
        index_wellenlaenge - pixbereich):(index_wellenlaenge + pixbereich)]
    geschwindigkeitsbereich = (wellenlaengenbereich - wellenlaenge) / \
        wellenlaenge * 299710

    plt.figure(1, figsize=(7, 10))
    plt.plot(geschwindigkeitsbereich, fluxbereich, linewidth=.5)
    plt.xlabel('Geschwindigkeit relativ zur Laborwellenlänge in km/s')
    plt.ylabel('relative Intensität')
    plt.title(obj + ', ' + linie)
    plt.pause(.01)

    plt.figure(2, figsize=(7, 10))
    plt.plot(geschwindigkeitsbereich, fluxbereich +
             i * offset, "-", label=JD, linewidth=1)
    # Beschriftung der einzelnen Spektren:
    plt.text(geschwindigkeitsbereich[1], 1.0 +
             i * offset, format(JD, '.2f'), ha="left", size=7)
    plt.xlabel('Geschwindigkeit relativ zur Laborwellenlänge in km/s')
    plt.ylabel('relative Intensität')
    plt.title(obj + ', ' + linie)
    plt.axvline(x=0, color='k', linewidth=.05)
    plt.grid(visible=True, axis='x')
    plt.pause(.01)

frage = input(
    'Möchten Sie die Grafiken abspeichern als pdf? Wenn ja "y" eingeben: ')
if frage == 'y':
    plt.figure(1)
    plt.savefig('Overplot_' + obj + '_'+linie + '.pdf')
    plt.savefig('Overplot_' + obj + '_'+linie + '.png')
    plt.figure(2)
    plt.savefig('OverplotMitOffset' + obj + '_'+linie + '.pdf')
    plt.savefig('OverplotMitOffset' + obj + '_'+linie + '.png')

plt.close(fig=1)
plt.close(fig=2)
print('\nProgramm ist beendet')
