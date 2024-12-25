#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
timeseries_einPlot_mitOffset_Wellenlaengenbereich.py

Liest eine Zeitserie von 1d-Spektren im fits-Format ein.
Es ist ein Wellenlängenbereich wählbar.
Plottet die Spektrenausschnitte mit einem Offset übereinander und speichert
den plt als .png und .pdf ab.

Stand 20221105
@author: Lothar Schanne
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import glob

plt.switch_backend('Qt5Agg')  # Interaktives GrafikBackend einschalten

# Create file list. All spectra in one (sub)folder.
files = input("Path and name of the files (use wildcards) : ")
filelist = glob.glob(files)

# Aalphabetical sorting. If the name is correct, this results in a
# temporal order
filelist.sort()


# Print the list
print("\nList of spectra: \n")
print("Number of spectra: ", len(filelist), "\n")


# Eingabe des abgebildeten Wellenlängenbereichs
lambda_anfang = float(input('Geben Sie den Anfang des abzubildenden \
Wellenlängenbereichs an: '))
lambda_ende = float(input('Geben Sie das Ende des abzubildenden \
Wellenlängenbereichs an: '))

# Parameter für den Abstand zwischen den Spektren im zweiten Plot
offset = float(input("Please enter the desired offset: "))
obj = input("Please enter the object name: ")

fig = plt.figure(1, (7, 12))
zaehler = 0

# Read header and flux
for i in range(len(filelist)):
    print(filelist[i])
    flux, header = fits.getdata(filelist[i], header=True)
    # Generate the spectrum sections:
    step = header["CDELT1"]
    refpix = header["CRPIX1"]

    wave_erstesPix = header["CRVAL1"] - step * (refpix - 1)

    wave = np.zeros(header["NAXIS1"], dtype=float)
    wave_bereich = np.array([])
    flux_bereich = np.array([])
    for k in range(header["NAXIS1"]):
        wave[k] = wave_erstesPix + k * step
        if wave[k] >= lambda_anfang and wave[k] <= lambda_ende:
            wave_bereich = np.hstack([wave_bereich, wave[k]])
            flux_bereich = np.hstack([flux_bereich, flux[k]])
    plt.plot(wave_bereich, flux_bereich + zaehler * offset, "k-", linewidth=1.)
    plt.text(wave_bereich[1], 1.0 + zaehler *
             offset, filelist[i], ha="left", size=7)
    zaehler = zaehler + 1


# Customize plot properties (grid, labels, etc.):
# plt.ylim(0.2, 1.0 + zaehler * offset + zusatz)
plt.title("Time Series of " + obj)
plt.grid(True)
plt.xlabel("Wavelength in Angström")
# fig.savefig(obj + '_'+ str(lambda_anfang)+'_'+ str(lambda_ende) +
#             '.pdf', format='pdf')
fig.savefig(obj + '_' + str(lambda_anfang)+'_' + str(lambda_ende) +
            '.png', format='png')

plt.pause(.1)

print('Zum beenden des Programms in das zuletzt geöffnete Diagramm klicken.')
plt.waitforbuttonpress(-1)
plt.close('all')
