#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
terrLinien_RVKorrektur_RBF_1d_spectra_fits_series.py

Input: Einlesen einer Spektrenserie im fits-Format. Diese dürfen nicht
heliozentrisch korrigiert sein.

Das Skript korrigiert die Kalibrierung durch Bestimmung des Linienminimums von
bekannten terrestrischen Linien und rechnet die Spektren mit der ermittelten
mittleren RV um, soweit sie betragmäßig 3 km/s überschreitet.
Abspeichern der korrigierten Spektren als fits. Zeigen der modellierten
terrestrischen Linien als Grafik. Ausdrucken der ermittelten RV's. Schreiben
der RV's, des Mittelwerts und der Standardabweichung in ein ascii-file.

@author: Lothar Schanne
Stand 20221116
"""

import numpy as np
from astropy.io import fits, ascii
import glob
from PyAstronomy import pyasl
import matplotlib.pylab as plt
from scipy.interpolate import Rbf


# Liste der verwendeten terrestrischen Linien:
Linien = [
    6543.912,
    6552.646,
    6574.88]
print("Verwendete terrestrische Linien: ",
      Linien)

object = input('Geben Sie den Namen des Objekts ein (ohne blanks): ')

# Create file list. Spectra in a (sub)folder.
files = input("Path and name of the spectra (use wildcards) : ")
filelist = glob.glob(files)

# Alphabetical sorting. If the name is correct (e.g. 20180922-xyz.fit),
# this results in a temporal order.
filelist.sort()

# Print the list
print("\nList of spectra: \n")
print("Number of spectra: ", len(filelist), "\n")


miniflux = np.zeros([len(filelist), len(Linien)])
miniwave = np.zeros([len(filelist), len(Linien)])
RV = np.zeros([len(filelist), len(Linien)])
RV_mean = np.zeros(len(filelist))
RV_std = np.zeros(len(filelist))


# Processing the filelist:
for i in range(len(filelist)):
    flux, header = fits.getdata(filelist[i], header=True)

    print("\n", filelist[i])
    print("Original begin of wavelength scale, CRVAL1: ", header["CRVAL1"])

    if "NAXIS" in header:
        print("Dimension, NAXIS:                        ", header["NAXIS"])
    else:
        print("Das ist kein 1d-Spektrum !")
    if "NAXIS1" in header:
        nax = header["NAXIS1"]
        print("Anzahl der Werte (Abszisse), NAXIS1:     ", nax)
    else:
        print("NAXIS1 fehlt im header !")
    if "CRVAL1" in header:
        crval = header["CRVAL1"]
        print("Anfangs-Wellenlänge, CRVAL1:             ", crval)
    else:
        print("CRVAL1 fehlt im header !")
    if "CRPIX1" in header:
        crpix = header["CRPIX1"]
        print("Referenzpixel, CRPIX1:             ", crpix)
    else:
        print("CRPIX1 fehlt im header !")
        crpix = 1
    if "CDELT1" in header:
        cdel = header["CDELT1"]
        print("Schrittweite der Wellenlänge, CDELT1:    ", cdel)
    else:
        print("CDELT1 fehlt im header !")

    #   Erzeugen von numpy-Arrays mit den Wellenlängen und Fluxes des Spektrums
    wave = np.ones(nax, dtype=float)
    for k in range(nax):
        wave[k] = crval + (k - crpix + 1) * cdel

    # Berechnen des mittleren Kalibrierfehlers anhand der terrestrischen Linien

    for m in range(len(Linien)):
        wellenlaenge = Linien[m]
        index_wellenlaenge = int((wellenlaenge - wave[0]) / cdel)
        # **********************************************************
        # Breite des Suchintervalls, Breite in Angström, anpassen.
        suchintervall = int(1 / cdel)
        # ************************************************************
        intervallflux = np.zeros(suchintervall)
        intervallwave = np.zeros(suchintervall)
        for j in range(suchintervall):
            intervallflux[j] = flux[j +
                                    index_wellenlaenge - int(suchintervall / 2)]
            intervallwave[j] = (
                wave[0] + (j + index_wellenlaenge -
                           int(suchintervall / 2)) * cdel
            )

        # Radial basis function (RBF) über das Spektrum im Intervall
        # smooth anpassen, 0. = Funktion geht durch alle Punkte, >0. = ausgleich
        rbf = Rbf(intervallwave, intervallflux, smooth=0.1)

        newwaveintervall = np.arange(
            intervallwave[0], intervallwave[-1], cdel / 10)
        newfluxinterpolated = rbf(newwaveintervall)
        miniflux[i, m] = newfluxinterpolated.min()
        # linienminimum_wave[i] = intervallwave[0] + \
        #     newfluxinterpolated.argmin()*step/10
        miniwave[i, m] = newwaveintervall[newfluxinterpolated.argmin()]

        RV[i, m] = (miniwave[i, m] - wellenlaenge) / wellenlaenge * 299710

        # Plotten
        fig = plt.figure()
        plt.plot(intervallwave, intervallflux, "o")
        plt.plot(newwaveintervall, newfluxinterpolated, "-")
        plt.plot(miniwave[i, m], miniflux[i, m], "or")
        # plt.savefig(filelist[i].rsplit(".", 1)[0] + "_" + Linien[m] + "_RBF.png")

        # RV
        RV[i, m] = (miniwave[i, m] - wellenlaenge) / wellenlaenge * 299710
        print("RV: ", RV[i, m].round(2))

    RV_mean[i] = RV[i].mean()
    RV_std[i] = RV[i].std()
    print('Mittlere RV: ', RV_mean[i], '\nStandardabweichung RV: ', RV_std[i])

    # Shift the spectrum unter der Bedingung, dass |RV_mean| > 3 km/s ist
    if abs(RV_mean[i]) >= 3.:
        flux_rv, wave_rv = pyasl.dopplerShift(
            wave, flux, -RV_mean[i], edgeHandling="firstlast")

        # Sxhreiben des RV-korrigierten Spektrums in ascii-file
        # ascii.write(
        #     [wave, flux_rv],
        #     filelist[i].rstrip(".fits") + "_perRBF_RVcorrected.dat",
        #     overwrite=True,
        #     names=["WAVE", "FLUX"],
        #     format="tab",
        # )

        # Schreiben des RV-korrigierten Spektrums in fits-file
        header["CRVAL1"] = wave[0]
        header["CRPIX1"] = 1
        header["NAXIS1"] = len(wave)
        newfile = filelist[i].rsplit(".", 1)[0] + "_perRBF_RVcorrected.fits"
        fits.writeto(newfile, flux_rv, header, overwrite=True,
                     output_verify="silentfix")

        print("Neue Anfangswellenlänge CRVAL1: ", header["CRVAL1"], "\n\n")

print('\nMittlere RVs der Spektren aus den terrestrischen Linien: ', RV_mean, '\n')

# Speichern der RV's als ascii-Datei (csv):
ascii.write(
    [filelist, RV[:, 0], RV[:, 1], RV[:, 2], RV_mean, RV_std],
    object + "_terrLinienKorrigiert" + ".csv",
    overwrite=True,
    names=[
        "Spektrum",
        "RV[0]", 'RV[1]', 'RV[2]',
        "RV_mean", 'RV_std'
    ],
    format="csv",
)
