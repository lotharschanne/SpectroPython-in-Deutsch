#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Liest einen Spektrenkatalog ein (Zeitserie im fits-Format).
Berechnet aus dem Beobachtungszeitpunkt
und den (anzupassenden) Koordinaten des Beobachters und Objekts die
heliozentrische Korrektur, zeigt die gewählte Linie in einem Plot,
das Minimum der gewählten Linie ist anzuklicken und bestimmt die
heliozentrisch korrigierte Radialgeschwindigkeit RV_bc.
Gibt ermittelte Daten als csv-ascii-Dateien (Komma-separiert) aus.

Stand 20221105
@author: lothar
"""

import numpy as np
from astropy.io import fits, ascii
import glob
import matplotlib.pylab as plt
from PyAstronomy import funcFit as fuf
from PyAstronomy import pyTiming as pyt
from PyAstronomy import pyasl
from astroplan import FixedTarget
from astropy.coordinates import SkyCoord
import astropy.units as u
from time import sleep

# lokaler Modul, muss im gleichen Verzeichnis stehen wie das Skript (oder das
# Verzeichnis des Moduls im Pythonpath, damit Python ihn findet)
import Linienlisten


plt.switch_backend('Qt5Agg')  # Interaktives GrafikBackend einschalten


# ************* Koordinaten des Observatory in folgender Form: *************
# longitude = 289.5967661 in Grad oder [grad, min, sec]
# latitude = -24.62586583 in Grad oder [grad, min, sec]
# altitude = 2635.43 in Meter

# ***************  BITTE AUF EIGENE KOORDINATEN ÄNDERN ***************
# Falls das versäumt wird sind die baryzentrischen Korrekturen falsch

# Koordinaten vom Berthold
# longitude = +7.4775
# latitude = 49.47527777778
# altitude = 200
# longitude = [7, 28, 39.0]
# latitude = [49, 28, 31.0]


# Koordinaten des Wise Observatory in Israel
# longitude = +34.76333333
# latitude = 30.59583333
# altitude = 875

# Koordinaten von Siegfried Hold
longitude = +15.68461111
latitude = 47.00161111
altitude = 380


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

Frage = input(
    'Möchten Sie die Sternkoordinaten im Internet suchen lassen? Dann "y" eingeben: ')
if Frage == 'y':
    star = input('Geben Sie den Namen des Objektsterns ein: ')
    ra2000 = FixedTarget.from_name(star).ra.value
    dec2000 = FixedTarget.from_name(star).dec.value
# ********************************************************************

# Fileliste erstellen. Spektren in einem (Unter)Ordner.
files = input("Pfad und Name der Spektren (nutze wildcards) : ")
filelist = glob.glob(files)

# Alphabetisches Sortieren. Bei richtiger Namensgebung ergibt das eine
# zeitliche Ordnung
filelist.sort()

# Ausdruck der Liste
print("\nSpektrenliste: \n")
print("Anzahl der Spektren: ", len(filelist), "\n")


# Auswahl der Linie oder Wellenlänge
linie, laborwellenlaenge = Linienlisten.linienauswahl()

# Eingabe der Systemgeschwindigkeit
sys = input(
    "Möchten Sie eine Systemgeschwindigkeit eingeben? Zahl eingeben in km/s oder n: "
)
if sys == "n":
    sysRVlambda = 0.0
else:
    sysRVlambda = float(sys) * laborwellenlaenge / 299710

# Abarbeiten der filelist, Einlesen von flux und header, Auswahl des Flux um die Linie:
hjd = np.zeros(len(filelist))
corr = np.zeros(len(filelist))
obs_time = np.zeros(len(filelist))
RV = np.zeros(len(filelist))
RV_bc = np.zeros(len(filelist))
apex = np.zeros(len(filelist))

for i in range(len(filelist)):
    flux, header = fits.getdata(filelist[i], header=True)
    if "JD" in header:
        obs_time[i] = header["JD"]
    elif "JD-OBS" in header:
        obs_time[i] = header["JD-OBS"]
    elif "JD_OBS" in header:
        obs_time[i] = header["JD_OBS"]
    elif "MJD-OBS" in header:
        mjd = header["MJD-OBS"]
        obs_time[i] = mjd + 2400000.5
    elif "BAS_MJD" in header:
        obs_time[i] = header["BAS_MJD"]
    else:
        print("!!!!! Kein Beobachtungszeitpunkt im Header  !!!!")

    # Berechnung der heliozentrischen Korrektur:
    corr[i], hjd[i] = pyasl.helcorr(
        longitude, latitude, altitude, ra2000, dec2000, obs_time[i], debug=False
    )

    print("\n" + filelist[i] + ":")
    print("Beobachtungszeitpunkt: ", obs_time[i])
    print("Barycentric correction [km/s]: ", corr[i])
    # print("Heliocentric Julian day: ", hjd[i])
    lambda0 = header["CRVAL1"]
    if "CRPIX1" in header:
        refpix = header["CRPIX1"]
    else:
        refpix = 1
    step = header["CDELT1"]

    wave = np.ones(header["NAXIS1"], dtype=float)
    header["CRVAL1"] = header["CRVAL1"] + (1 - refpix) * step
    for j in range(header["NAXIS1"]):
        wave[j] = header["CRVAL1"] + j * step + sysRVlambda

    # Plot the spectrum
    fig = plt.figure()
    plt.plot(wave, flux)
    plt.xlabel("Wavelength [Angström]", fontsize="medium")
    plt.ylabel("Flux")
    # plt.ylim(0.7, 1.03)
    plt.title("Spectrum " + filelist[i])
    plt.grid(True)
    plt.xticks(fontsize="small")
    plt.yticks(fontsize="small")
    plt.xlim(laborwellenlaenge + sysRVlambda - 5,
             laborwellenlaenge + sysRVlambda + 5)

    print("Klicke das Minimum der Linie und danach die Basis (das Kontinuum) an: ")
    pts = []
    pts = np.asarray(plt.ginput(n=2, timeout=-1))

    RV[i] = (pts[0, 0] - laborwellenlaenge) / laborwellenlaenge * 299710
    RV_bc[i] = RV[i] + corr[i]
    apex[i] = pts[0, 1] / pts[1, 1]
    print("\nbaryzentrisch (und systemisch) korrigierte RV: ",
          RV_bc[i], "Apex: ", apex[i])
    plt.pause(.1)


# # Plot der RV's
# fig=plt.figure()
# plt.plot(obs_time, RV,'bo', markersize=1)
# plt.plot(obs_time, RV_bc,'r+', markersize=1)

# Abspeichern als ascii-Datei
ascii.write(
    [filelist, obs_time, RV_bc],
    linie + "_RV_bc" + ".csv",
    overwrite=True,
    names=['Spectrum',  "JD",  "RV_bc"],
    format="csv",
)
ascii.write(
    [filelist, obs_time, apex],
    linie + "_apex" + ".csv",
    overwrite=True,
    names=['Spectrum', "JD", "apex"],
    format="csv",
)
# Speichern von obs_time und corr in ascii-file
ascii.write(
    [filelist, obs_time, corr],
    linie + "_Tabelle_obs_time_bc.csv",
    overwrite=True,
    names=['Spectrum', "JD", "BARYCORR"],
    format="csv",
)

fr = input('Wenn Sie fertig mit der Betrachtung der Grafiken sind und das\
 Programm beenden möchten, drücken sie die Eingabe-Taste.')
plt.close('all')
print('\nProgramm ist beendet')
