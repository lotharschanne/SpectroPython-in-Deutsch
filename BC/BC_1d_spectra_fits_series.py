#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BC_1d_spectra_fits_series.py

Das Skript berechnet die heliozentrische Korrektur BC und das korrigierte JD als
HJD einer Serie von 1d-Spektren im Fits-Format eines Objekts.
Schreibt eine ascii-Datei (Format-Tab) mit den BC's und HJD's in das Arbeits-
verzeichnis.
Korrigiert auf Wunsch die eingelesenen Spektren um die berechnete
heliozentrische Korrektur BC und speichert sie als fits ab (und nach
auskommentieren der Zeilen 226 bis 233 auch als asci-Tabelle).

Wichtig: Koordinaten des Beobachters anpassen !!!
Die Sternkoordinaten können im Skript manuell angepasst oder aus dem Internet
eingelesen werden.

@author: Lothar Schanne
Stand 20221117
"""

from PyAstronomy import pyasl
import numpy as np
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table
import glob
from astroplan import FixedTarget


# ************* Koordinaten des Observatory in folgender Form: *************
# longitude = 289.5967661 in Grad oder [grad, min, sec]
# latitude = -24.62586583 in Grad oder [grad, min, sec]
# altitude = 2635.43 in Meter
# oder Eingabe der Koordinaten im hexagesimalen System Grad Minuten Sekunden
# in Form einer Liste:
# longitude = [7, 28, 39.0]
# latitude = [49, 28, 31.0]

# ***************  BITTE AUF EIGENE KOORDINATEN ÄNDERN ***************
# Falls das versäumt wird sind die baryzentrischen Korrekturen falsch

# Koordinaten vom Berthold
# longitude = +7.4775
# latitude = 49.47527777778
# altitude = 200
# longitude = [7, 28, 39.0]
# latitude = [49, 28, 31.0]


# Koordinaten vom Bernd Bitnar
# longitude = +13.708425
# latitude = 51.003557
# altitude = 270

# # Koordinaten vom Ulrich Waldschläger
# longitude = +13.533
# latitude = 52.533
# altitude = 100

# Koordinaten vom Uwe Zurmühl
# longitude = +9.9036
# latitude = +52.1993
# altitude = 100

# Koordinaten des Wise Observatory in Israel
# longitude = +34.76333333
# latitude = 30.59583333
# altitude = 875

# Koordinaten von Siegfried Hold
longitude = +15.68461111
latitude = 47.00161111
altitude = 460

# Koordinaten VEGA-Sternwarte Haunsberg, Salzburg, Herbert
# latitude = 47.923971244662944
# longitude = 13.00748535806445
# altitude = 790


# Wichtig !!!!!!!!!!!!!!!! :
# esign : int, optional, {-1,0,1}
# Explicit sign with -1 representing negative sign, +1 representing positive
# sign, and 0 indicating no explicit sign specification. The explicit sign is
# necessary if negative southern coordinates are specified but d is 0 and,
# thus, cannot carry the sign.

if type(longitude) == list:
    longitude = pyasl.dmsToDeg(longitude[0], longitude[1], longitude[2])

if type(latitude) == list:
    latitude = pyasl.dmsToDeg(latitude[0], latitude[1], latitude[2], esign=+1)

# ********************************************************************

# ***************  BITTE AUF EIGENE KOORDINATEN ÄNDERN ***************
# ************ Eingabe der Koordinaten des Sterns in folgender Form: *******
# Falls das versäumt wird sind die baryzentrischen Korrekturen falsch
# ra2000 = 030.20313477 in Grad
# dec2000 = -12.87498346 in Grad
# oder als Koordinaten-String RA DEC in der Form "hh mmm ss +dd mm ss"

# # Koordinaten von del Cep
# ra2000 = 337.29277083333335
# dec2000 = +58.415198
# coord = "22 29 10.265 +58 24 54.714"

# Koordinaten von gam Cyg
# ra2000 = 305.55708
# dec2000 = +40.2566

# Koordinaten von Betageuze
# ra2000 = 88.7929583
# dec2000 = 7.40705555

# Koordinaten von theta1 Ori C
# ra2000 = 83.81858333
# dec2000 = -5.389694444

# Koordinaten von 7 And
# ra2000 = 348.1375
# dec2000 = 49.40620275
# coord = "23 12 33 +49 24 22.3299"

# Koordinaten von gam Cyg
# ra2000 = 305.55708
# dec2000 = 40.25666666

# Koordinaten von OX Aurigae
# ra2000 = 103.25
# dec2000 = 38.86916
# coord = "06 53 01.41099 +38 52 08.9353"

# if type(coord) == str:
#     ra2000, dec2000 = pyasl.coordsSexaToDeg(coord)

# Einlesen der Sternkoordinaten über das Internet
Frage = input(
    'Möchten Sie die Sternkoordinaten im Internet suchen lassen? Dann "y" eingeben: ')
if Frage == 'y':
    star = input('Geben Sie den Namen des Objektsterns ein: ')
    ra2000 = FixedTarget.from_name(star).ra.value
    dec2000 = FixedTarget.from_name(star).dec.value
# ********************************************************************


# Create file list. Spectra in a (sub)folder.
files = input('Path and name of the spectra (use wildcards) : ')
filelist = glob.glob(files)

# Alphabetical sorting. If the name is correct (e.g. 20180922-xyz.fit),
# this results in a temporal order.
filelist.sort()

# Print the list
print('\nList of spectra: \n')
print('Number of spectra: ', len(filelist), '\n')

frage = input(
    'Möchten Sie die Spektren baryzentrisch korrigieren? Dann y eingeben: ')


corr = np.zeros(len(filelist))
hjd = np.zeros(len(filelist))
obs_time = np.zeros(len(filelist))

# Processing the filelist, reading flux and header:
for i in range(len(filelist)):
    flux, header = fits.getdata(filelist[i], header=True)
    # Check for the necessary header entries
    if 'JD' in header:
        obs_time[i] = float(header["JD"])
    elif "JD-OBS" in header:
        obs_time[i] = float(header["JD-OBS"])
    elif "JD_OBS" in header:
        obs_time[i] = float(header["JD_OBS"])
    elif "MJD-OBS" in header:
        mjd = header["MJD-OBS"]
        obs_time[i] = mjd + 2400000.5
    elif "BAS_MJD" in header:
        obs_time[i] = float(header["BAS_MJD"]) + 2400000.5
    else:
        print("Es ist kein Beobachtungszeitpunkt im Header von ", filelist[i])
        break

    print('\nSpectrum and observation date: ', filelist[i], obs_time[i])

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

    #   Erzeugen von numpy-Arrays mit den Wellenlängen des Spektrums
    wave = np.ones(nax, dtype=float)
    for k in range(nax):
        wave[k] = crval + (k - crpix + 1) * cdel

    # Calculate heliocentric correction (debug=True show
    # various intermediate results)
    corr[i], hjd[i] = pyasl.helcorr(longitude, latitude, altitude,
                                    ra2000, dec2000, obs_time[i], debug=False)
    print("Heliocentric correction [km/s]: ", corr[i])
    print("Heliocentric Julian day: ", hjd[i])

    # Shift the spectrum
    if frage == 'y':
        flux_bc, wave_bc = pyasl.dopplerShift(
            wave, flux, corr[i], edgeHandling="firstlast")

        # Schreiben des RV-korrigierten Spektrums in ein ascii-file
        # ascii.write(
        #     [wave, flux_bc],
        #     filelist[i].rstrip(".fits") + "_BCcorrected.dat",
        #     overwrite=True,
        #     names=["WAVE", "FLUX"],
        #     format="tab",
        # )

        # Schreiben des RV-korrigierten Spektrums in fits-file
        header["CRVAL1"] = wave[0]
        header["CRPIX1"] = 1
        header["NAXIS1"] = len(wave)
        newfile = filelist[i].rstrip(".fits") + "_BCcorrected.fits"
        fits.writeto(newfile, flux_bc, header, overwrite=True,
                     output_verify="ignore")

        print("Neue Anfangswellenlänge CRVAL1: ", header["CRVAL1"], "\n\n")

data = Table([filelist, obs_time, hjd, corr], names=[
             "Spektrum", "Beobachtungsdatum [JD]", 'HJD', 'BC'])
ascii.write(data, "BC_Liste" + ".dat", overwrite=True, format="tab")
