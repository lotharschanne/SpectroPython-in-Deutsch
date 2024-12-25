#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Liest Spektrenkatalog ein (Zeitserie von normierten 1d-Spektren im
fits-Format). Berechnet aus dem Beobachtungszeitpunkt und den (anzupassenden)
Koordinaten des Beobachters und Objekts die heliozentrische Korrektur,
fitted die gewählte Linie per Spline
und bestimmt aus dem heliozentrisch korrigierten Minimum die
heliozentrisch korrigierte Radialgeschwindigkeit RV.
Plottet und speichert alle fittings und gibt ermittelte Daten als ascii-Dateien
(Komma-separiert, als .csv) aus.

Stand 20220408

@author: lothar
"""

import numpy as np
from astropy.io import fits, ascii
import glob
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from PyAstronomy import pyasl
from scipy import interpolate
from astroplan import FixedTarget
from astropy.coordinates import SkyCoord
import astropy.units as u

# lokaler Modul, im Ordner "Bausteine" zu finden,
# muss im gleichen Verzeichnis stehen wie das Skript (oder das
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

# Einlesen der Sternkoordinaten über das Internet
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
linie, wellenlaenge = Linienlisten.linienauswahl()

frage_bary = input(
    'Möchten Sie die RV baryzentrisch korrigieren? Dann Eingabe von "j"')

frage_grafikenspeichern = input(
    'Wenn Sie die Grafiken speichern möchten, geben Sie j ein: ')

# Abarbeiten der filelist, Einlesen von flux und header, Auswahl des Flux um die Linie:

# Definition von Variablen
hjd = np.zeros(len(filelist))
corr = np.zeros(len(filelist))
obs_time = np.zeros(len(filelist))
linienminimum_flux = np.zeros(len(filelist))
linienminimum_wave = np.zeros(len(filelist))
linienminimum_bc = np.zeros(len(filelist))
RV = np.zeros(len(filelist))
RV_bc = np.zeros(len(filelist))


for i in range(len(filelist)):
    flux, header = fits.getdata(filelist[i], header=True)

    if "JD-OBS" in header:
        obs_time[i] = float(header["JD-OBS"])
    elif 'JD' in header:
        obs_time[i] = float(header["JD"])
    elif "JD_OBS" in header:
        obs_time[i] = float(header["JD_OBS"])
    elif "MJD-OBS" in header:
        mjd = header["MJD-OBS"]
        obs_time[i] = mjd + 2400000.5
    elif "BAS_MJD" in header:
        obs_time[i] = float(header["BAS_MJD"])
    else:
        print("Es ist kein Beobachtungszeitpunkt im Header von ",
              filelist[i])
        break
    # Berechnung der heliozentrischen Korrektur und Zeit:
    if frage_bary == 'j':
        corr[i], hjd[i] = pyasl.helcorr(
            longitude, latitude, altitude, ra2000, dec2000, obs_time[i], debug=False
        )
    else:
        corr[i] = 0
        hjd[i] = 0

    print("\n" + filelist[i] + ":")
    print("Beobachtungszeitpunkt: ", obs_time[i])
    print("Barycentric correction [km/s]: ", corr[i])

    if "CRPIX1" in header:
        refpix = int(header["CRPIX1"])
    else:
        refpix = 1

    step = float(header["CDELT1"])

    lambda0 = float(header["CRVAL1"]) - (refpix - 1) * step

    index_wellenlaenge = int((wellenlaenge - lambda0) / step)

    wave = np.zeros(header['NAXIS1'], dtype=float)
    for k in range(header['NAXIS1']):
        wave[k] = lambda0 + k * step

    # **********************************************************
    # Breite des Suchintervalls, Breite in Angström, anpassen.
    suchintervall = int(8 / step)
    # ************************************************************
    intervallflux = np.zeros(suchintervall)
    intervallwave = np.zeros(suchintervall)
    for j in range(suchintervall):
        intervallflux[j] = flux[j +
                                index_wellenlaenge - int(suchintervall / 2)]
        # intervallwave[j] = (
        #     lambda0 + (j + index_wellenlaenge - int(suchintervall / 2)) * step
        intervallwave[j] = wave[j +
                                index_wellenlaenge - int(suchintervall / 2)]

    # Spline über das Spektrum im Intervall
    tck = interpolate.splrep(intervallwave, intervallflux, s=0.0001)
    newwaveintervall = np.arange(
        intervallwave[0], intervallwave[-1], step / 50)
    newfluxinterpolated = interpolate.splev(newwaveintervall, tck, der=0)
    linienminimum_flux[i] = newfluxinterpolated.min()
    # linienminimum_wave[i] = intervallwave[0] + \
    #     newfluxinterpolated.argmin()*step/10
    linienminimum_wave[i] = newwaveintervall[newfluxinterpolated.argmin()]
    linienminimum_bc[i] = linienminimum_wave[i] * (1 + corr[i] / 299710)

    RV[i] = (linienminimum_wave[i] - wellenlaenge) / wellenlaenge * 299710
    print('\n' + filelist[i])
    print('unkorrigierte RV: ', RV[i])
    # RV nur bc-korrigiert, sysv nicht berücksichtigt:
    RV_bc[i] = (linienminimum_bc[i] - wellenlaenge) / wellenlaenge * 299710
    print("baryzentrisch korrigierte RV: ", RV_bc[i])

    # Plotten
    fig = plt.figure()
    plt.title(filelist[i] + ', JD: ' + str(obs_time[i]), fontsize=7)
    plt.plot(intervallwave, intervallflux, "o")
    plt.plot(newwaveintervall, newfluxinterpolated, "-")
    plt.plot(linienminimum_wave[i], linienminimum_flux[i], "o", color="red")
    plt.pause(1)
    if frage_grafikenspeichern == 'j':
        plt.savefig(filelist[i].rstrip(".fit") + "_" + linie + "_Spline.png")
    plt.close()

# # Plot der RV's
# fig = plt.figure()
# plt.plot(obs_time, RV, 'bo', markersize=1)

# Abspeichern als ascii-Datei
ascii.write(
    [filelist, obs_time, corr, RV, RV_bc, linienminimum_flux],
    linie + "_RV_perSpline" + ".csv",
    overwrite=True,
    names=[
        "      Spektrum     ",
        "      JD      ",
        "     BC     ",
        "     RV     ",
        "      RV_bc     ",
        "      miniFlux",
    ],
    format="csv",
)


# ascii.write([obs_time, apex3], linie+'_apex'+'.dat', overwrite=True,
#             names=['JD', 'apex'], format='tab')


# # Speichern von obs_time und corr in ascii-file
# ascii.write([obs_time, corr], linie+'_Tabelle_obs_time_bc.dat', overwrite=True,
#             names=['JD', 'BARYCORR'], format='tab')

# # Speichern von obs_time und jhd in ascii-file
# ascii.write([obs_time, hjd], linie+'_Tabelle_obs_time_hjd.dat', overwrite=True,
#             names=['JD', 'HJD'], format='tab')


plt.close('all')
print('\nProgramm ist beendet')
