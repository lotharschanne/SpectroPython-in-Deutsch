#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Einlesen eines synthetischen Spektrums in Form einer ascii-Tabelle mit 2 Spalten 
WAVE und FLUX überschrieben.
Berechnung eines wählbaren Wellenlängenausschnitts. Dieser Bereich wird dann mit
einer wählbaren Schrittweite rebinned und anschließend noch zusätzlich mit einer
wählbaren FWHM (Apparateprofil) gefaltet.
Geplottet wird der gewählte rebinnte Wellenlängenbereich und zusätzlich das
gefaltete Spektrum.
Der rebinnte Flux-Ausschnitt des ursprünglichen ascii-Datei und der
rebinnte und convolvierte Flux wird zusammen mit der Wellenlänge in je einer
zweispaltigen ascii-Tabelle (spacer = Komma) und in je einer fits-Datei
abgespeichert.

Stand 20231113
@author: lothar
"""

import matplotlib.pyplot as plt
from astropy.convolution import Gaussian1DKernel, convolve
from astropy.io import ascii, fits
from PyAstronomy.pyasl import binningx0dt


plt.switch_backend('Qt5Agg')  # Interaktives GrafikBackend einschalten

file = input("Pfad und Filebezeichnung eingeben: ")
table = ascii.read(file)

begin = float(
    input("Geben Sie den Beginn des gewünschten Wellenlängenbereichs in Angström\
 ein: "))
end = float(
    input("Geben Sie das Ende des gewünschten Wellenlängenbereichs ein: "))

fwhm = float(
    input("Geben Sie die gewünschte FWHM des Apparateprofils in Angström ein: ")
)
step = float(
    input("Geben Sie die gewünschte Schrittweite in Angström des \
Ergebnisspektrums ein: "))

newtable = table[table["WAVE"] >= begin]
newtable = newtable[newtable["WAVE"] <= end]

# Binning
data, d = binningx0dt(
    newtable["WAVE"],
    newtable["FLUX"],
    dt=step,
    x0=newtable["WAVE"].min(),
    removeEmpty=False,
)
Data = data.T
wave = Data[0]
flux = Data[1]

#   Convolve mit astropy.convolution
kernel = Gaussian1DKernel(stddev=fwhm / 2.3 / step)
convoluted = convolve(flux, kernel, normalize_kernel=True, boundary="extend")

# Grafik
# plt.style.use('seaborn-whitegrid')
fig, ax = plt.subplots(2)
plt.xlabel("Wellenlänge [Angström]", fontsize=5)
fig.suptitle("Pollux-Spektrum " + file, fontsize=5)
ax[0].plot(wave, flux, linewidth=0.2)
# ax[0].set_xlim(4000, 6700)
# ax[0].set_ylim(0.2, 1.1)
ax[0].tick_params(axis="both", labelsize=5)
ax[1].plot(wave, convoluted, linewidth=0.2)
ax[1].tick_params(axis="both", labelsize=5)
# ax[1].set_xlim(4000, 6700)
# ax[1].set_ylim(0.2, 1.1)
name = file.rstrip(".csv")  # anpassen !!!
plt.savefig(file + ".png")
plt.savefig(name + ".pdf")

# Speichern des beschnittenen ascii, FLUX als dat und fits

# nicht convolviert
name = (
    file.rstrip(".csv") + "_" + str(int(begin)) +  # anpassen !!!!
    "_" + str(int(end)) + "_rebinned.dat"
)
ascii.write([wave, flux], name, overwrite=True,
            names=["WAVE", "FLUX"], format="tab")

name = (
    file.rstrip(".csv")  # anpassen !!!
    + "_"
    + str(int(begin))
    + "_"
    + str(int(end))
    + "_rebinned.fits"
)

header = fits.Header()
header["SIMPLE"] = "T"
header["BITPIX"] = -32
header["NAXIS"] = 1
header["CRVAL1"] = wave[0]
header["NAXIS1"] = len(wave)
header["CDELT1"] = step
header["CUNIT1"] = "Angstrom"
header["CTYPE1"] = "Wavelength"
header["CRPIX1"] = 1

fits.writeto(
    name, flux, header, overwrite=True, output_verify="silentfix",
)

# convolviert
name = (
    file.rstrip(".csv")   # abpassen !!!
    + "_"
    + str(int(begin))
    + "_"
    + str(int(end))
    + "_rebinned_convolved.csv"
)

ascii.write(
    [wave, convoluted], name, overwrite=True, names=["WAVE", "FLUX"], format="csv"
)

name = (
    file.rstrip(".csv")   # anpassen
    + "_"
    + str(int(begin))
    + "_"
    + str(int(end))
    + "_rebinned_convolved.fits"
)

header = fits.Header()
header["SIMPLE"] = "T"
header["BITPIX"] = -32
header["NAXIS"] = 1
header["CRVAL1"] = wave[0]
header["NAXIS1"] = len(wave)
header["CDELT1"] = step
header["CUNIT1"] = "Angstrom"
header["CTYPE1"] = "Wavelength"
header["CRPIX1"] = 1

fits.writeto(
    name, convoluted, header, overwrite=True, output_verify="silentfix",
)
