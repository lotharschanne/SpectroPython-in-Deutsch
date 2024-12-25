#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Extrahiert aus einer Serie von NRES-Dateien (-e92-1d.fits)
fits-Dateien der 68 Ordnungen, wellenlängenkalibriert und normiert, und speichert
sie im fits-Format ab.
Alle mit einheitlichem binning (0,05 Angström).

@author: lothar
20241207
"""

# notwendige Bibliotheken/Module
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io import ascii
from PyAstronomy.pyasl import binningx0dt
import glob

# Fileliste erstellen. Spektren in einem (Unter)Ordner.
files = input('Pfad und Name der Spektren (nutze wildcards) : ')
filelist = glob.glob(files)

# Alphabetisches Sortieren. Bei richtiger Namensgebung ergibt das eine
# zeitliche Ordnung
filelist.sort()

# Ausdruck der Liste
print('\nSpektrenliste: \n')
print(filelist)
print('Anzahl der Spektren: ', len(filelist), '\n')

for j in np.arange(len(filelist)):

    HDU = fits.open(filelist[j], ignore_missing_end=True)
    data = HDU[1].data
    HEADER = HDU[0].header
    sciencefiber = HEADER['SCIFIBER']

    for i in np.arange(len(data)):
        if data[i]['fiber'] == sciencefiber:
            flux = data[i]['normflux']
            wave = data[i]['wavelength']
            newflux = flux[wave != 0]
            newwave = wave[wave != 0]

            # neues Binning
            newbinned_data, dt = binningx0dt(
                newwave, newflux, x0=newwave[0], dt=0.05)

            hdr = fits.Header()
            hdr['NAXIS1'] = len(newbinned_data)
            hdr['CRVAL1'] = newbinned_data[0, 0]
            hdr['CDELT1'] = dt
            hdr['CRPIX1'] = 1
            hdr['DATE-OBS'] = (HEADER['DATE-OBS'],
                               'UTC] Start date and time of observation')
            hdr['MJD-OBS'] = (HEADER['MJD-OBS'],
                              'Start date/time (Modified Julian Date)')
            hdr['OBJECT'] = HEADER['OBJECT']
            hdr['SITE'] = (HEADER['SITE'], 'Site name')
            hdr['SITEID'] = (
                HEADER['SITEID'], 'ID code of the Observatory site')
            hdr['TIMESYS'] = (HEADER['TIMESYS'], 'Time system used')
            hdr['DATE-OBS'] = (HEADER['DATE-OBS'],
                               '[UTC] Start date and time of observation')
            hdr['EXPTIME'] = (HEADER['EXPTIME'], '[s] Exposure length')
            hdr['USERID'] = (HEADER['USERID'], 'User ID')
            hdr['SNR'] = (
                HEADER['SNR'], 'Signal-to-noise ratio/res element @ 5180 Angstr')
            hdr['RV'] = (
                HEADER['RV'], 'Radial Velocity in Barycentric Frame [m/s]')
            hdr['RVERR'] = (HEADER['RVERR'],
                            'Radial Velocity Uncertainty [m/s]')
            hdr['BARYCORR'] = (
                HEADER['BARYCORR'], ' Barycentric Correction Applied to the RV [m/s]')
            hdr['TCORR'] = (HEADER['TCORR'],
                            'Exposure Mid-Time (Barycentric Julian Date)')

            hdu = fits.PrimaryHDU(header=hdr)
            hdu.data = newbinned_data[:, 1]
            hdul = fits.HDUList([hdu])

            hdul.writeto(filelist[j].rsplit('.', 1)[0]+'_order_' +
                         str(data[i]['order'])+'_norm.fits', overwrite=True)

    HDU.close()

print('Ende des Programs')
