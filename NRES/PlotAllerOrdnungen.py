#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Liest die Ordnungen von NRES-Spektren im fits-Format ein.

Plottet die Spektren und speichert den Plot als .png und .pdf ab

Stand 20241211
@author: Lothar Schanne
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import glob
from PyAstronomy import pyasl

plt.switch_backend('Qt5Agg')  # Interaktives GrafikBackend einschalten

# Create file list. All spectra in one (sub)folder.
files = input("Path and name of the fits-files (use wildcards) : ")
filelist = glob.glob(files)

filelist.sort()

# Print the list
print("\nList of spectra: \n")
print("Number of spectra: ", len(filelist), "\n")


fig = plt.figure(1, figsize=(20, 10))

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

    wave_erstesPix = header["CRVAL1"] - step * (refpix - 1)
    wave = np.zeros(header["NAXIS1"], dtype=float)
    for k in range(header["NAXIS1"]):
        wave[k] = wave_erstesPix + k * step

    plt.plot(wave, flux, linewidth=1.)
    # Beschriftung der einzelnen Spektren:
    plt.text(wave[100], flux[i], filelist[i].rsplit('.')[0], ha="left", size=7)
    # plt.xlim(6300,6500)
    plt.pause(.05)

plt.title('Ordnungen für ' + filelist[0].rsplit('.', 1)[0])
plt.grid(True)
# plt.ylim((0., 1.2))
plt.xlabel("Wavelength in Angström")
plt.pause(0.1)

fig.savefig("timeseries.png", format='png')
fig.savefig("timeseries.pdf", format='pdf')

print('Ende des Programms')
