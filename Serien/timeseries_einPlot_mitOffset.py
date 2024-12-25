#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Liest eine Serie von 1d-Spektren im fits-Format ein.

Plottet die Spektren mit einem Offset übereinander und speichert den plot als
.png und .pdf ab

Stand 20221105
@author: Lothar Schanne
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import glob

plt.switch_backend('Qt5Agg')  # Interaktives GrafikBackend einschalten

# Create file list. All spectra in one (sub)folder.
files = input("Path and name of the fits-files (use wildcards) : ")
filelist = glob.glob(files)

# Aalphabetical sorting. If the name is correct, this results in a
# temporal order
filelist.sort()

# Print the list
print("\nList of spectra: \n")
print("Number of spectra: ", len(filelist), "\n")

# Parameter für den Abstand zwischen den Spektren im zweiten Plot
offset = float(input("Please enter the desired offset: "))
obj = input("Please enter the object name: ")

fig = plt.figure(1, figsize=(7, 10))

# Read header and flux
for i in range(len(filelist)):
    print(filelist[i], ":")
    flux, header = fits.getdata(filelist[i], header=True)
    # Generate the spectrum sections and save them as fit:
    step = header["CDELT1"]
    if 'CRPIX1' in header:
        refpix = header["CRPIX1"]
    else:
        refpix = 1
    # JD der Beobachtung, zur Beschriftung der Spektren verwendet:
    if 'JD' in header:
        JD = header['JD']
    elif 'MJD-OBS' in header:
        JD = header['MJD-OBS']
    elif 'DATE-OBS' in header:
        JD = header['DATE-OBS']
    else:
        print('\nKein Beobachtungsdatum im header')

    wave_erstesPix = header["CRVAL1"] - step * (refpix - 1)
    if i == len(filelist) - 1:
        zusatz = max(flux) - 1
    wave = np.zeros(header["NAXIS1"], dtype=float)
    for k in range(header["NAXIS1"]):
        wave[k] = wave_erstesPix + k * step
    plt.plot(wave, flux + i * offset, "k-", label=JD, linewidth=.3)
    # Beschriftung der einzelnen Spektren:
    plt.text(wave[1], flux[1] + i * offset, JD, ha="left", size=5)
    # plt.xlim(6300,6500)
    plt.pause(.1)


# Customize plot properties (grid, labels, etc.):
plt.ylim(0., 1.0 + len(filelist) * offset + zusatz)
plt.title("Time Series of " + obj)
plt.grid(True)
plt.xlabel("Wavelength in Angström")
plt.pause(0.1)

# fig.savefig("timeseries.png", format='png')
fig.savefig("timeseries.pdf", format='pdf')

# print('Zum beenden des Programms in das zuletzt geöffnete Diagramm klicken.')
# plt.waitforbuttonpress(-1)
# plt.close('all')
