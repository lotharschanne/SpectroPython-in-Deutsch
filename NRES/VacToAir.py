#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Das Skript wandelt die Wellenlängenskala von im Vakuum gemessenen Spektren um
in die Wellenlängen, wie sie in der Luft gemessen worden wären.

Version: 20241207

@author: lothar schanne
"""

import numpy as np
from astropy.io import fits
import glob
from PyAstronomy import pyasl


def anfangswellenlaenge_calc(header):
    """
    Berechnet die Wellenlänge des ersten Pixels aus den Headerdaten der fits-Datei
    """
    if "CRPIX1" not in header:
        header["CRPIX1"] = 1.
    header['CRVAL1'] = header["CRVAL1"] + \
        (1 - header["CRPIX1"]) * header["CDELT1"]
    header["CRPIX1"] = 1.
    return header


# Create file list. All spectra in one (sub)folder.
files = input("Path and name of the fits-files (use wildcards) : ")
filelist = glob.glob(files)

filelist.sort()

# Print the list
print("\nList of spectra: \n")
print("Number of spectra: ", len(filelist), "\n")

for i in range(len(filelist)):
    flux, header = fits.getdata(filelist[i], header=True)

    header = anfangswellenlaenge_calc(header)

    # # Umwandlung Wellenlängen von Vakuum zu Luft
    header['CRVAL1'] = pyasl.vactoair2(header['CRVAL1'])

    print('Dateiname: ', filelist[i])
    name = filelist[i].rstrip('.',1)[0] + '_air.fits'
    print('neuer Dateiname: ', name)
    fits.writeto(name, flux, header, overwrite=True)
