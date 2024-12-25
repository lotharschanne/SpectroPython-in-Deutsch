#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
timeseries_einPlot_mitOffset.py

Liest eine Zeitserie von 1d-Spektren im fits-Format ein. Es wird ein JD-Bereich
gewählt, der geplottet werden soll.
Plottet die Spektren mit einem Offset übereinander und speichert den plt als
.png und .pdf ab.

Stand 20221105
@author: Lothar Schanne
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import glob


# Create file list. All spectra in one (sub)folder.
files = input("Path and name of the file (use wildcards) : ")

filelist = glob.glob(files)

# Aalphabetical sorting. If the name is correct, this results in a
# temporal order
filelist.sort()

# Print the list
print("\nList of spectra: \n")
print("Number of spectra: ", len(filelist), "\n")

JD = np.zeros(len(filelist))
for i in range(len(filelist)):
    flux, header = fits.getdata(filelist[i], header=True)
    JD[i] = header['JD']
print('JD-Bereich von ', JD.min(), ' bis ', JD.max())

# Eingabe des Bereichs der Beobachtungszeitpunkte als JD
JD_anfang = float(input('Geben Sie den Anfang des Bereichs der \
Beobachtungszeitpunkte (JD) ein: '))
JD_ende = float(input('Geben Sie das Ende des Bereichs der \
Beobachtungszeitpunkte (JD) ein: '))

# Parameter für den Abstand zwischen den Spektren im zweiten Plot
offset = float(input("Please enter the desired offset: "))
# obj = input("Please enter the object name: ")

obj = 'Algol JD' + str(JD_anfang) + ' bis ' + str(JD_ende)

fig = plt.figure(1, figsize=(7, 10))
xlim_a = 4330  # Anpassen !!!!!!
xlim_e = 4350  # Anpassen !!!!!!
plt.xlim(xlim_a, xlim_e)
zaehler = 0
# Read header and flux
for i in range(len(filelist)):
    flux, header = fits.getdata(filelist[i], header=True)
    JD = header['JD']
    JD_float = float(JD)
    if JD_float >= JD_anfang and JD_float <= JD_ende:
        print(filelist[i], ":")
        step = header["CDELT1"]
        refpix = header["CRPIX1"]
        wave_erstesPix = header["CRVAL1"] - step * (refpix - 1)
        zusatz = max(flux) - 1
        wave = np.zeros(header["NAXIS1"], dtype=float)
        for k in range(header["NAXIS1"]):
            wave[k] = wave_erstesPix + k * step
        plt.plot(wave, flux + zaehler * offset, "-", linewidth=1)
        # Beschriftung der einzelnen Spektren:
        plt.text(xlim_a, 1 + zaehler * offset, JD, ha="left", size=8)
        # plt.xlim(6300,6500)
        plt.pause(0.1)
        zaehler += 1
    else:
        pass


# Customize plot properties (grid, labels, etc.):
# plt.ylim(0.2, 1.0 + (zaehler * offset) + zusatz)
plt.title("Time Series of " + obj)
plt.grid(True)
plt.xlabel("Wavelength in Angström")
plt.ylabel('normierter Flux')
plt.pause(0.1)


frage = input(
    'Möchten Sie die Grafik abspeichern? Wenn ja "y" eingeben: ')
if frage == 'y':
    plt.savefig('ofsettedOverplot_' + obj + '.pdf')
    plt.savefig('ofsettedOverplot_' + obj + '.png')

e = input('Zum beenden des Programms die Taste e drücken.')
if e == 'e':
    plt.close('all')
    print('Ende des Programms')
