#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SNR_fits_mehrereBereiche.py
Lädt ein 1d-Spektrum im fits-Format und plottet es.
Auswahl von Spektrenausschnitten durch 2 Mausklicks und anschließendem
zweimaligen Drücken des Escape-Buttons.
Die SNR der Bereiche werden ausgedruckt und in eine kommaseparierte csv-Datei
gespeichert.

Stand 20221105
@author: Lothar Schanne
"""

import numpy as np
from astropy.io import fits, ascii
import matplotlib.pyplot as plt


plt.switch_backend('Qt5Agg')  # Interaktives GrafikBackend einschalten


# Path and name of the fits-file
file = input("Path and name of the fits-file: ")

#   Read header and data (flux)
flux, header = fits.getdata(file, header=True)

print("Minimum and Maximum in flux: ", flux.min(), "  ", flux.max())


#   Check for the necessary header entries
nax = header["NAXIS1"]
crval = header["CRVAL1"]
cdel = header["CDELT1"]
if "CRPIX1" not in header:
    header["CRPIX1"] = 1

#   Generation of a numpy array with the wavelengths of the spectrum
wave = np.ones(nax, dtype=float)
crval = crval + (1 - header["CRPIX1"]) * cdel
for i in range(nax):
    wave[i] = crval + i * cdel
# The wave list contains the wavelengths of the pixels.
# In the list flux the corresponding intensities.

# Plot spectrum
fig = plt.figure(1, figsize=(14, 10))
plt.plot(wave, flux)
plt.xlabel("Wavelength [Angstroem]")
plt.ylabel("ADU")
plt.title("Spectrum " + file)
plt.grid(True)
plt.pause(10)  # Zeit zum Vergrößern des gewünschten Spektrumausschnitts

# Interactive setting of the base points begin and end for SNR estimation
# Press the ESC key twice to end the interactivity
# plt.waitforbuttonpress()

snr = []
Anfang = []
Ende = []
eingabe = "y"
while eingabe == "y":
    print("nach Markieren der beiden Wellenlängen im Diagramm per linkem ")
    print("Mausklick zweimal escape-Taste drücken: ")
    pts = []
    pts = np.asarray(plt.ginput(n=-1, timeout=-1))
    if plt.waitforbuttonpress():
        pass
    plt.plot(pts[:, 0], pts[:, 1], "o", markersize=6)

    a = pts[0, 0]
    b = pts[1, 0]

    aindex = int((a - crval) / cdel)
    bindex = int((b - crval) / cdel)

    newflux = flux[aindex:bindex]
    newwave = wave[aindex:bindex]

    # Calculation of SNR:
    SNR = newflux.mean() / newflux.std()
    snr.append(SNR)
    Anfang.append(a)
    Ende.append(b)
    print("SNR zwischen %.2f und %.2f  = %.1f" % (a, b, SNR))
    print("\nNeuen Bereich eingeben? Dafür y eingeben.")
    print("ZurProgrammweiterführung einen anderen Buchstaben eingeben: ")
    eingabe = input()
    if eingabe != "y":
        break

# Abspeichern als ascii-Datei
ascii.write(
    [Anfang, Ende, snr],
    file + "_SNR.csv",
    overwrite=True,
    names=["Anfang", "Ende", "SNR"],
    format="csv",
)

SNR = 0
for i in range(len(snr)):
    SNR += snr[i]
SNR_mean = SNR / len(snr)
print("Mittleres SNR = ", SNR_mean)

plt.close('all')
