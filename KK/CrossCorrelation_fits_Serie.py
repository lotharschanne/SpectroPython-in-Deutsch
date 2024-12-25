#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Kreuzkorrelation
abgeleitet von einem Beispiel in PyAstronomy
https://www.hs.uni-hamburg.de/DE/Ins/Per/Czesla/PyA/PyA/pyaslDoc/aslDoc/
crosscorr.html

Es wird eine Kreuzkorrelation einer Serie von target-Spektren bzgl. eines
template-Spektrums durchgeführt. Beide liegen als fits vor. Die Objektspektren 
können, müssen aber nicht heliozentrisch korrigiert sein.
Die RV's und wahlweise heliozentrisch korrigierten RV's werden ausgedruckt und
in einer Datei abgespeichert.
Falls die baryzentrische korrigierten RV's berechnet werden sollen, müssen
die Beobachter- und die Objektkoordinaten im Skript angepasst werden. 
Standardmäßig werden die Sternkoordinaten aber durch Eingabe des Objektnamens 
(z.B. alp Ori) aus dem Internet übernommen.

Stand 20220408

@author: Lothar Schanne
"""


from PyAstronomy import pyasl
import numpy as np
import matplotlib.pylab as plt
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table
import glob
import time
from astroplan import FixedTarget
from astropy.coordinates import SkyCoord
import astropy.units as u

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
altitude = 200
# Hier Koordinaten in Form einer Liste
longitude = [7, 28, 39.0]
latitude = [49, 28, 31.0]


# Koordinaten des Wise Observatory in Israel
# longitude = +34.76333333
# latitude = 30.59583333
# altitude = 875

# Koordinaten von Siegfried Hold
# longitude = +15.68461111
# latitude = 47.00161111
# altitude = 380


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

# Koordinaten von Polaris (alp UMi)
# coord = "02 31 49 +89 15 51"

# if type(coord) == str:
#     ra2000, dec2000 = pyasl.coordsSexaToDeg(coord)

# Einlesen der Sternkoordinaten über das Internet
star = input('Geben Sie den Namen des Objektsterns ein: ')
Frage = input(
    'Möchten Sie die Sternkoordinaten im Internet suchen lassen? Dann "y" eingeben: ')
if Frage == 'y':
    ra2000 = FixedTarget.from_name(star).ra.value
    dec2000 = FixedTarget.from_name(star).dec.value
# ********************************************************************


# Template auswählen
# Pfad und Name des templates
tfile = input("Pfad und Filebezeichnung des template eingeben: ")

#   Einlesen von Header und Daten (Flux vom template in tf,
#   Header des template in theader gespeichert)
tf, theader = fits.getdata(tfile, header=True)
print("Minimum und Maximum im Template [ADU]: ", tf.min(), "  ", tf.max())

tnax = theader["NAXIS1"]
tcrval = theader["CRVAL1"]
tcdel = theader["CDELT1"]

tw = np.ones(tnax, dtype=float)
tcrval = tcrval + (1 - theader["CRPIX1"]) * tcdel
for i in range(tnax):
    tw[i] = tcrval + i * tcdel


# Zu korrelierende Spektren (targets) auswählen
# Pfad und Name der targets
# Create file list.
files = input("\nPath and name of the target spectra (use wildcards) : ")
filelist = glob.glob(files)

# Sort alphabetically. If the spectrum files are named correctly, this results
# in a temporal order.
filelist.sort()

# Printout of the list for control purposes.
print("\nList of target spectra:")
print(filelist)
print("\nNumber of target spectra: ", len(filelist), "\n")
time.sleep(1.)
frage_bary = input(
    'Möchten Sie die RV baryzentrisch korrigieren? Dann Eingabe von "j": ')
print("Bitte warten. Berechnung läuft\n")

RV = np.zeros(len(filelist))
RV_bc = np.zeros(len(filelist))
hjd = np.zeros(len(filelist))
corr = np.zeros(len(filelist))
obs_time = np.zeros(len(filelist))

for i in range(len(filelist)):
    f, header = fits.getdata(filelist[i], header=True)
    # print('Minimum und Maximum im'+filelist[i]+' [ADU]: ', f.min(), '  ', f.max())
    nax = header["NAXIS1"]
    crval = header["CRVAL1"]
    cdel = header["CDELT1"]
    if "JD" in header:
        obs_time[i] = float(header["JD"])
    else:
        if "JD-OBS" in header:
            obs_time[i] = float(header["JD-OBS"])
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
    if frage_bary == 'j':
        # Berechnung der heliozentrischen Korrektur und Zeit:
        corr[i], hjd[i] = pyasl.helcorr(
            longitude, latitude, altitude, ra2000, dec2000, obs_time[i],
            debug=False
        )

    #   Erzeugen eines numpy-Arrays mit den Wellenlängen des target
    w = np.ones(nax, dtype=float)
    crval = crval + (1 - header["CRPIX1"]) * cdel
    for j in range(nax):
        w[j] = crval + j * cdel

    # Cross-Correlation ausführen.
    # Das RV-range ist (Parameter 1) - bis (Parameter 2) km/s in
    # Schritten von (Parameter 3) km/s.
    # Die ersten und letzten (Parameter 4) Datenpunkte sind unberücksichtigt.
    # Muss angepasst werden, damit die Spektren noch überlappen können bei der
    # maximalen Verschiebung (Parameter 2)
    rv, cc = pyasl.crosscorrRV(
        w, f, tw, tf, -200.0, 200.0, 0.1, mode="doppler", skipedge=10
    )

    # Maximum der cross-correlation function finden
    maxind = np.argmax(cc)
    RV[i] = rv[maxind]
    RV_bc[i] = RV[i] + corr[i]

    # Plot CroCo-Funktion
    fig = plt.figure()
    plt.plot(rv, cc, "b-")
    plt.plot(rv[maxind], cc[maxind], "ro")
    plt.title("Kreuzkorrelationsfunktion " + filelist[i], size=8)
    plt.grid(True)
    plt.xlabel("km/s")
    plt.pause(1)
    # plt.savefig(filelist[i]+'_CroCo.png')
    plt.close()

print("Spektrum                   ", '  JD    ',
      '     HZK   ', " RV [km/s] ", "  RV_bc [km/s]")
for i in range(len(filelist)):
    print(filelist[i], obs_time[i], corr[i], RV[i], RV_bc[i])

data = Table([filelist, obs_time, corr, RV, RV_bc], names=[
             "Spektrum", "JD", 'baryz. Korrektur', "RV", 'RV_bc'])
ascii.write(data, "KK_RV_Liste_" + star + ".dat", overwrite=True, format="tab")


plt.close('all')
print("Ende des Programs")
