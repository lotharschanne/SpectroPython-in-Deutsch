#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RVKorrektur_1d_spectra_fits_series.py

The script corrects a series of 1d-spectra in fits format of an object
by an manually editetd RV [km/s].
Writes fit-files and ascii-files of the RV-corrected spectra into the
working directory. The respective RV is noted in the the header of the
generated fit.

@author: Lothar Schanne
Stand 20211225
"""

import numpy as np
from astropy.io import fits, ascii
import glob
from PyAstronomy import pyasl


# Create file list. Spectra in a (sub)folder.
files = input("Path and name of the spectra (use wildcards) : ")
filelist = glob.glob(files)

# Alphabetical sorting. If the name is correct (e.g. 20180922-xyz.fit),
# this results in a temporal order.
filelist.sort()

# Print the list
print("\nList of spectra: \n")
print("Number of spectra: ", len(filelist), "\n")

rv = float(input("Geben Sie die Korrekturgeschwindigkeit ein in km/s: "))

# Processing the filelist, reading flux and header:
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

    # Shift that spectrum
    flux_rv, wave_rv = pyasl.dopplerShift(
        wave, flux, rv, edgeHandling="firstlast")

    ascii.write(
        [wave, flux_rv],
        filelist[i].rstrip(".fit") + "_RVcorrected.dat",
        overwrite=True,
        names=["WAVE", "FLUX"],
        format="tab",
    )
    # Schreiben des RV-korrigierten Spektrums in fits-file
    header["CRVAL1"] = wave[0]
    header["CRPIX1"] = 1
    header["NAXIS1"] = len(wave)
    newfile = filelist[i].rstrip(".fit") + "_RVcorrected.fit"
    header["RV_COR"] = (rv, "km/s, corrected")
    fits.writeto(newfile, flux_rv, header, overwrite=True,
                 output_verify="silentfix")

    print("Neue Anfangswellenlänge CRVAL1: ", header["CRVAL1"], "\n\n")
