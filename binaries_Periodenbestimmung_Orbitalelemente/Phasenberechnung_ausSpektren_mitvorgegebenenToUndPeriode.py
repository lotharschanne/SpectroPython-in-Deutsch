#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Das Skript liest einen Katalog von fits-1d-Spektren ein. 
Mit den im Skript vorgegebenen (und evtl. anzupassenden) Parametern (Periode 
und T0) werden dann die Beobachtunsgszeitpunkte gefaltet, und die Ergebnisse 
(Phase, JD) in einem ascii-File (Format csv) gespeichert.
Die Spektren müssen einen Headereintrag namens 'JD' enthalten.

Created on 20231219

@author: lothar
"""
from PyAstronomy import pyasl
import matplotlib.pyplot as plt
import pandas as pd
from astropy.io import fits
import numpy as np
import glob

plt.switch_backend('Qt5Agg')  # Interaktives GrafikBackend einschalten

Periode = 680.168  # anpassen !!!!!!!!!!!!!!!!!!!!!!
T0 = 2454433.2     # anpassen !!!!!!!!!!!!!!!!!!!!!!

# Create file list. All spectra in one (sub)folder.
files = input('Path and name of the fits-files (use wildcards) : ')
filelist = glob.glob(files)

# Alphabetical sorting. If the name is correct, this results in a
# temporal order
filelist.sort()

# Print the list
print('\nList of spectra: \n')
print('Number of spectra: ', len(filelist), '\n')
jd = np.zeros(len(filelist))

for i in range(len(filelist)):
    flux, header = fits.getdata(filelist[i], header=True)
    try:
        jd[i] = header['JD']
    except:
        print('Spektrum '+filelist[i]+'hat keinen JD-Eintrag')
        jd[i] = None
        pass

# Folding (Phase)
phases = pyasl.foldAt(jd, Periode, T0)


data = pd.DataFrame({'Spektrum': filelist, 'JD': jd, 'Phase': phases})
datafile = data.to_csv('Phasen_JD_AusSpektren_perToUndPeriode.csv')
