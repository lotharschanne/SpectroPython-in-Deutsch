#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 14:04:06 2024

@author: lothar
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import glob
from PyAstronomy import pyasl


def wavecalc(flux, header):
    """Berechnet die Wellenlängen aus den Headerdaten der fits-Datei."""
    wave = np.zeros(header["NAXIS1"])
    if "CRPIX1" not in header:
        header["CRPIX1"] = 1
    header["CRVAL1"] = header["CRVAL1"] + \
        (1 - header["CRPIX1"]) * header["CDELT1"]
    for i in np.arange(header["NAXIS1"]):
        wave[i] = header["CRVAL1"] + i * header["CDELT1"]
    return wave, flux


plt.switch_backend('Qt5Agg')  # Interaktives GrafikBackend einschalten


# Create file list for objects. All spectra in one (sub)folder.
files = input('Path and name of the OBJECT files (use wildcards) : ')
filelist = glob.glob(files)
filelist.sort()

print('Number of spectra: ', len(filelist), '\n')

for i in range(len(filelist)):
    flux, header = fits.getdata(filelist[i], header=True)
    wave, flux = wavecalc(flux, header)
    print('\nObjektspektrum: ', filelist[i])

    iin, iout = pyasl.slidingPolyResOutlier(wave, flux, points=10, deg=2, stdlim=2.5,
                                            mode='above', dx=5)

    # What about the outliers
    print("Number of outliers: ", len(iout))
    print("Indices of outliers: ", iout)

    plt.figure()
    plt.plot(wave, flux)
    plt.scatter(wave[iout], flux[iout], c='r')

    # number, indices = pyasl.pointDistGESD(flux, 20, alpha=0.001)
    # # What about the outliers
    # print("Number of outliers: ", number)
    # print("Indices of outliers: ", indices)

    # plt.figure()
    # plt.plot(wave, flux)
    # plt.scatter(wave[indices], flux[indices], c='r')
