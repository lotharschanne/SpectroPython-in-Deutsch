#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Berechnet für eine Serie von 1d-fits-Spektren mittels des beta*sigma-Verfahrens
das SNR jedes Spektrums. Speichert die SNR in einer Datei (kommasepariert, .csv).


Stand 20221105

@author: lothar
"""

from PyAstronomy import pyasl
import numpy as np
from astropy.io import fits, ascii
import glob

# Fileliste erstellen. Spektren in einem (Unter)Ordner.
files = input('Pfad und Name der Spektren (nutze wildcards) : ')
filelist = glob.glob(files)

# Alphabetisches Sortieren. Bei richtiger Namensgebung ergibt das eine
# zeitliche Ordnung
filelist.sort()

# Ausdruck der Liste
print('\nSpektrenliste: \n')
print('Anzahl der Spektren: ', len(filelist), '\n')

fluxmean = np.zeros(len(filelist))
nstd = np.zeros(len(filelist))
nstdstd = np.zeros(len(filelist))
snr = np.zeros(len(filelist))

for i in range(len(filelist)):
    flux, header = fits.getdata(filelist[i], header=True)
    # Estimate noise using robust estimate
    beq = pyasl.BSEqSamp()
    # Define order of approximation (use larger values such as 2,3, or 4 for
    # faster varying or less well sampled data sets; also 0 is a valid order)
    N = 1
    # Define 'jump parameter' (use larger values such as 2,3, or 4 if correlation
    # between adjacent data point is suspected)
    j = 3
    # Estimate noise assuming equidistant sampling (often a good approximation even
    # if data are not strictly equidistant) and robust estimation (often advantageous
    # in working with real data)
    nstd[i], nstdstd[i] = beq.betaSigma(flux, N, j, returnMAD=True)

    fluxmean[i] = flux.mean()
    snr[i] = fluxmean[i] / nstd[i]

    print('\nSpektrum: ', filelist[i])
    print("Estimated noise std = %7.5f +/- %7.5f" % (nstd[i], nstdstd[i]))
    print('SNR = ', int(snr[i]))
    print('für die Parameter N = ', N, ' und j = ', j)

ascii.write([filelist, fluxmean, nstd, nstdstd, snr], files+'_SNR'+'.csv', overwrite=True,
            names=['Spektrum', 'mittlerer Flux', 'noise', 'noise_std', 'SNR'], format='csv')
