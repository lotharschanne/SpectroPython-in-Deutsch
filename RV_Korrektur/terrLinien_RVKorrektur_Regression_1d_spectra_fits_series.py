#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
terrLinien_RVKorrektur_1d_spectra_fits_series.py

Input: Einlesen einer Spektrenserie im fits-Format. Diese dürfen nicht
heliozentrisch korrigiert sein.

Das Skript korrigiert die Kalibrierung durch Bestimmung des Linienminimums von
bekannten terrestrischen Linien per Regression und rechnet die Spektren mit der
ermittelten mittleren RV um.
Abspeichern der korrigierten Spektren als fits. Zeigen der modellierten
terrestrischen Linien als Grafik. Ausdrucken der ermittelten RV's. Schreiben
der RV's, des Mittelwerts und der Standardabweichung in ein ascii-file.

@author: Lothar Schanne
Stand 20221117
"""

import numpy as np
from astropy.io import fits, ascii
import glob
from PyAstronomy import pyasl
import matplotlib.pylab as plt

#################################################################
# Regression, Grad anpassen (2 oder 4 oder 6)
grad = 2
#################################################################

# Liste der wählbaren terrestrischen Linien:
Linien = [
    6543.912,
    6552.646,
    6574.88]
print("Verwendete terrestrische Linien: ",
      Linien)

object = input('Geben Sie den Namen des Objekts ein: ')

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
        suchintervall = int(3 / cdel)
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

        linienminimum_wave_index = intervallflux.argmin()

        ##########################################################################
        # Das folgende Intervall anpassen, bei verrauschten Linien breit (-20,+21)
        # bei nicht verrauschten, schmalen Linien klein (-2,+3)
        ##########################################################################
        linienminimum_intervall = np.arange(
            linienminimum_wave_index - 1, linienminimum_wave_index + 2
        )

        # Wellenlängen und Flux direkt um das Linienminimum
        wl = intervallwave[linienminimum_intervall]
        fl = intervallflux[linienminimum_intervall]

        model = np.poly1d(np.polyfit(wl, fl, grad))

        polyline = np.linspace(wl[0], wl[-1], 100)
        modflux = model(polyline)

        # Linienminimum berechnen
        miniflux[i, m] = modflux.min()
        miniwave = polyline[modflux.argmin()]
        print('')

        # Plotten
        fig = plt.figure()
        plt.plot(intervallwave, intervallflux, "o-")
        plt.plot(polyline, modflux)
        plt.plot(miniwave, miniflux[i, m], "o", color="black")
        plt.title(filelist[i] + 'Linie' + str(Linien[m]))
        plt.xlabel("Wellenlänge [Angström]")
        plt.ylabel("relative Intensität")
        # plt.savefig(filelist[i].rsplit(".", 1)[0] + "_" + Linie + "_Regression.png")

        # RV
        RV[i, m] = (miniwave - wellenlaenge) / wellenlaenge * 299710
        print("RV: ", RV[i, m].round(2))

    RV_mean[i] = RV[i].mean()
    RV_std[i] = RV[i].std()

    # Shift that spectrum
    flux_rv, wave_rv = pyasl.dopplerShift(
        wave, flux, -RV_mean, edgeHandling="firstlast")

    ascii.write(
        [wave, flux_rv],
        filelist[i].rsplit(".", 1)[0] + "_RVcorrected.dat",
        overwrite=True,
        names=["WAVE", "FLUX"],
        format="tab",
    )

    # Schreiben des RV-korrigierten Spektrums in fits-file
    header["CRVAL1"] = wave[0]
    header["CRPIX1"] = 1
    header["NAXIS1"] = len(wave)
    newfile = filelist[i].rsplit(".", 1)[0] + "_RVcorrected.fit"
    fits.writeto(newfile, flux_rv, header, overwrite=True,
                 output_verify="silentfix")

    print("Neue Anfangswellenlänge CRVAL1: ", header["CRVAL1"], "\n\n")
print('\nMittlere RVs der Spektren aus den terrestrischen Linien: \n', RV_mean, '\n')

# Speichern der RV's als ascii-Datei (csv):
ascii.write(
    [filelist, RV[:, 0], RV[:, 1], RV[:, 2], RV_mean, RV_std],
    object + "_terrLinienKorrigiert_Regression" + ".csv",
    overwrite=True,
    names=[
        "Spektrum",
        "RV[0]", 'RV[1]', 'RV[2]',
        "RV_mean", 'RV_std'
    ],
    format="csv",
)

print('Zum beenden des Programms in das zuletzt geöffnete Diagramm klicken.')
plt.waitforbuttonpress(-1)
plt.close('all')
