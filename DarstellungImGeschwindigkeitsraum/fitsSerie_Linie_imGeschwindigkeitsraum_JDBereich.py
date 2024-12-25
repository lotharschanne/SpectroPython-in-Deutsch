#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Einlesen einer Serie von heliozentrisch korrigierten Spektren im fits-Format,
Eingabe eines Zeitraums als JD Anfang und JD Ende.
Darstellung einer wählbaren Linie im Geschwindigkeitsraum in zwei Plots:
    Erstens alle Spektren übereinander geplottet.
    Zweitens alle Spektren mit einem wählbaren offset übereinander versetzt 
    geplottet.
Die plots können als pdf und png abgespeichert werden.

20230216

@author: lothar
"""

import numpy as np
from astropy.io import fits, ascii
import glob
import matplotlib.pylab as plt

# lokaler Modul, muss im gleichen Verzeichnis stehen wie das Skript (oder das
# Verzeichnis des Moduls im Pythonpath, damit Python ihn findet)
import Linienlisten

plt.switch_backend('Qt5Agg')  # Interaktives GrafikBackend einschalten
plt.style.use('seaborn-whitegrid')

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

# JD-Bereich der Serie
flux, header = fits.getdata(filelist[0], header=True)
print('JD des ersten Spektrums', header['JD'])
flux, header = fits.getdata(filelist[-1], header=True)
print('JD des letzten Spektrums', header['JD'])

# Eingabe des Bereichs der Beobachtungszeitpunkte als JD
JD_anfang = float(input('Geben Sie den Anfang des Bereichs der \
Beobachtungszeitpunkte (JD) ein: '))
JD_ende = float(input('Geben Sie das Ende des Bereichs der \
Beobachtungszeitpunkte (JD) ein: '))

bereich = float(input('Geben Sie den darzustellenden Geschwindigkeitsbereich um \
die gewählte Linie in km/s ein: '))

# Eingabe der Systemgeschwindigkeit:
systemgeschwindigkeit = float(
    input("Geben Sie die Systemgeschwindigkeit in km/s ein: "))
systemgeschwindigkeit = systemgeschwindigkeit / 299710 * wellenlaenge

# Parameter für den Abstand zwischen den Spektren im offset Plot fig1:
offset = float(input("Please enter the desired offset: "))
obj = input("Please enter the object name: ")

k = -1
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
    JD_float = float(JD)

    if JD_float >= JD_anfang and JD_float <= JD_ende:
        k += 1
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

        plt.figure(1, figsize=(7, 10))  # Overplot
        plt.plot(geschwindigkeitsbereich, fluxbereich, linewidth=.5)
        plt.xlabel('Geschwindigkeit relativ zur Laborwellenlänge in km/s')
        plt.ylabel('relative Intensität')
        plt.title(obj + ', ' + linie)
        plt.pause(.01)

        plt.figure(2, figsize=(7, 10))  # Overplot mit offset
        plt.plot(geschwindigkeitsbereich, fluxbereich +
                 k * offset, "-", label=JD, linewidth=1)
        # Beschriftung der einzelnen Spektren:
        plt.text(geschwindigkeitsbereich[1], 1.0 +
                 k * offset, format(JD, '.2f'), ha="left", size=7)
        plt.xlabel('Geschwindigkeit relativ zur Laborwellenlänge in km/s')
        plt.ylabel('relative Intensität')
        plt.title(obj + ', ' + linie)
        plt.axvline(x=0, color='k', linewidth=.05)
        plt.grid(visible=True, axis='x')
        plt.pause(.01)


frage = input(
    'Möchten Sie die Grafiken abspeichern als pdf und png? Wenn ja "y" eingeben: ')
if frage == 'y':
    # plt.figure(1)
    # plt.savefig('Overplot_' + obj + '_' + linie + '_JD_Bereich_' + str(JD_anfang) +
    #             '_' + str(JD_ende) + '.pdf')
    # plt.savefig('Overplot_' + obj + '_' + linie + '_JD_Bereich_' + str(JD_anfang) +
    #             '_' + str(JD_ende) + '.png')
    plt.figure(2)
    plt.savefig('OverplotMitOffset' + obj + '_'+linie + '_JD_Bereich_' +
                str(JD_anfang) + '_' + str(JD_ende) + '.pdf')
    plt.savefig('OverplotMitOffset' + obj + '_'+linie + '_JD_Bereich_' +
                str(JD_anfang) + '_' + str(JD_ende) + '.png')

plt.close(fig=1)
plt.close(fig=2)
print('\nProgramm ist beendet')
