#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
timeseries_einPlot_mitOffset.py

Liest eine Zeitserie von 1d-Spektren im tab-Format ein.

Plottet die Spektren mit einem Offset übereinander und speichert den plot als
.png und .pdf ab

Stand 20221105
@author: Lothar Schanne
"""

import numpy as np
from astropy.io import ascii
import matplotlib.pyplot as plt
import glob

plt.switch_backend('Qt5Agg')  # Interaktives GrafikBackend einschalten

# Create file list. All spectra in one (sub)folder.
files = input("Path and name of the file (use wildcards) : ")
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

# Read wave and flux
for i in range(len(filelist)):
    print(filelist[i], ":")
    spectrum = ascii.read(filelist[i], guess=True)
    if i == len(filelist) - 1:
        zusatz = max(spectrum['FLUX']) - 1
    plt.plot(spectrum['WAVE'], spectrum['FLUX'] +
             i * offset, "k-", linewidth=.6)
    # Beschriftung der einzelnen Spektren:
    plt.text(spectrum['WAVE'][1], spectrum['FLUX'][1] +
             i * offset,  filelist[i], ha="left", size=8)

    plt.pause(.05)


# Customize plot properties (grid, labels, etc.):
# plt.xlim(6520, 6600)
plt.ylim(-0.2, 1.2 + len(filelist) * offset + zusatz)
plt.title("Time Series of " + obj)
plt.grid(True)
plt.xlabel("Wavelength in Angström")
plt.ylabel('Flux')

plt.pause(0.1)
fig.savefig("timeseries.pdf", format='pdf')
fig.savefig("timeseries.png", format='png')
