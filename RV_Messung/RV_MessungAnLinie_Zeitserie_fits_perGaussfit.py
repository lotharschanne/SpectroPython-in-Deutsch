#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Liest Spektrenkatalog ein (Zeitserie von normierten (!!!)  1d-Spektren im
fits-Format). Berechnet aus dem Beobachtungszeitpunkt und den (anzupassenden)
Koordinaten des Beobachters und Objekts die heliozentrische Korrektur,
fitted die gewählte Linie per Gaussfit dreimal (gesamt und innerer Bereich)
und bestimmt aus dem heliozentrisch korrigierten Minimum die
heliozentrisch korrigierte Radialgeschwindigkeit RV.
Plottet und speichert alle fittings und gibt ermittelte Daten als ascii-Dateien
(Komma-separiert, als .csv) aus.

Stand 20221105

@author: lothar
"""

import numpy as np
from astropy.io import fits, ascii
import glob
import matplotlib.pyplot as plt
from PyAstronomy import funcFit as fuf
from PyAstronomy import pyasl
from astroplan import FixedTarget


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
altitude = 200
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

# if type(coord) == str:
#     ra2000, dec2000 = pyasl.coordsSexaToDeg(coord)

# Einlesen der Sternkoordinaten über das Internet
Frage = input(
    'Möchten Sie die Sternkoordinaten im Internet suchen lassen? Dann "y" eingeben')
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


# Eingabe der Systemgeschwindigkeit
systemgeschwindigkeit = float(
    input("Geben Sie eine Systemgeschwindigkeit in km/s  ein: ")
)
systemgeschwindigkeit = systemgeschwindigkeit / 299710 * wellenlaenge

frage_bary = input(
    'Möchten Sie die RV baryzentrisch korrigieren? Dann Eingabe von "j"')

frage_grafikenspeichern = input(
    'Wenn Sie die Grafiken speichern möchten, geben Sie j ein: ')

# Abarbeiten der filelist, Einlesen von flux und header, Auswahl des Flux um die Linie:
hjd = np.zeros(len(filelist))
corr = np.zeros(len(filelist))
obs_time = np.zeros(len(filelist))
linienminimum1 = np.zeros(len(filelist))
linienminimum1_bc = np.zeros(len(filelist))
linienminimum2 = np.zeros(len(filelist))
linienminimum2_bc = np.zeros(len(filelist))
linienminimum3 = np.zeros(len(filelist))
linienminimum3_bc = np.zeros(len(filelist))
RV1 = np.zeros(len(filelist))
RV2 = np.zeros(len(filelist))
RV3 = np.zeros(len(filelist))
RV_diff1 = np.zeros(len(filelist))
RV_diff2 = np.zeros(len(filelist))
apex3 = np.zeros(len(filelist))


for i in range(len(filelist)):
    flux, header = fits.getdata(filelist[i], header=True)
    if frage_bary == 'j':
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
        corr[i], hjd[i] = pyasl.helcorr(
            longitude, latitude, altitude, ra2000, dec2000, obs_time[i], debug=False
        )
        print("\n" + filelist[i] + ":")
        print("Beobachtungszeitpunkt: ", obs_time[i])
        print("Barycentric correction [km/s]: ", corr[i])
    else:
        corr[i] = 0

    lambda0 = header["CRVAL1"]
    if "CRPIX1" in header:
        refpix = header["CRPIX1"]
    else:
        refpix = 1
    step = header["CDELT1"]
    index_wellenlaenge = int(
        (wellenlaenge - lambda0 + systemgeschwindigkeit) / step + refpix - 1
    )

    # **********************************************************
    # Breite des Suchintervalls, Breite in Angström, anpassen.
    suchintervall1 = int(2 / step)
    # ************************************************************
    intervallflux1 = np.zeros(suchintervall1)
    intervallwave1 = np.zeros(suchintervall1)
    for j in range(suchintervall1):
        intervallflux1[j] = flux[j +
                                 index_wellenlaenge - int(suchintervall1 / 2)]
        intervallwave1[j] = (
            lambda0 + (j + index_wellenlaenge - int(suchintervall1 / 2)) * step
        )

    fig = plt.figure()
    plt.plot(intervallwave1, intervallflux1, "k:")

    intervall1_flux_min = intervallflux1.min()
    intervallwave1_min = intervallwave1[intervallflux1.argmin()]

    gf1 = fuf.GaussFit1d()

    # Schätzwerte für das Gaußfitting:
    # A: - für Absorption, + für Emission
    gf1["A"] = -(1 - intervall1_flux_min)
    gf1["sig"] = 1.5  # auch der Wert muss an die Linienbreite angepasst werden
    gf1["off"] = 1.0
    gf1["mu"] = intervallwave1_min
    gf1["lin"] = 0.0

    gf1.thaw(["A", "sig", "off", "mu"])
    gf1.setRestriction({"A": [None, 0]})
    gf1.fit(intervallwave1, intervallflux1)
    # gf1.parameterSummary()

    linienminimum1[i] = gf1["mu"]
    if linienminimum1[i] <= wellenlaenge - 4:
        break
    if linienminimum1[i] >= wellenlaenge + 4:
        break
    linienminimum1_bc[i] = linienminimum1[i] * (
        1 + corr[i] / 299710
    )  # heliozentrische Korrektur
    RV1[i] = (linienminimum1_bc[i] - wellenlaenge) / \
        linienminimum1_bc[i] * 299710

    plt.plot(intervallwave1, gf1.model, "r--")

    # zweiter Durchlauf, eingegrenzt auf Intervall 2*sigma
    mu_index = int((gf1["mu"] - lambda0) / step + refpix - 1)
    suchintervall2 = 2 * abs(int(gf1["sig"] / step))
    intervallflux2 = np.zeros(suchintervall2)
    intervallwave2 = np.zeros(suchintervall2)
    for j in range(len(intervallflux2)):
        intervallflux2[j] = flux[j + mu_index - suchintervall2 // 2]
        intervallwave2[j] = lambda0 + \
            (j + mu_index - suchintervall2 // 2) * step
    gf2 = fuf.GaussFit1d()
    gf2.setRestriction({"A": [None, 0]})
    gf2["A"] = gf1["A"]
    gf2["sig"] = gf1["sig"]
    gf2["off"] = gf1["off"]
    gf2["mu"] = gf1["mu"]
    gf2.thaw(["A", "sig", "off", "mu"])
    gf2.fit(intervallwave2, intervallflux2)
    # gf2.parameterSummary()
    linienminimum2[i] = gf2["mu"]
    if linienminimum2[i] <= wellenlaenge - 4:
        break
    if linienminimum2[i] >= wellenlaenge + 4:
        break
    linienminimum2_bc[i] = linienminimum2[i] * (1 + corr[i] / 299710)
    RV2[i] = (linienminimum2_bc[i] - wellenlaenge) / \
        linienminimum2_bc[i] * 299710

    plt.plot(intervallwave2, gf2.model, "b--")

    RV_diff1[i] = RV2[i] - RV1[i]

    # dritter Durchlauf, eingegrenzt auf Intervall 2*sigma
    mu_index = int((gf2["mu"] - lambda0) / step + refpix - 1)
    suchintervall3 = abs(int(1.2 * gf2["sig"] / step))
    intervallflux3 = np.zeros(suchintervall3)
    intervallwave3 = np.zeros(suchintervall3)
    for k in range(len(intervallflux3)):
        intervallflux3[k] = flux[k + mu_index - suchintervall3 // 2]
        intervallwave3[k] = lambda0 + \
            (k + mu_index - suchintervall3 // 2) * step
    gf3 = fuf.GaussFit1d()
    gf3.setRestriction({"A": [None, 0]})
    gf3["A"] = gf2["A"]
    gf3["sig"] = gf2["sig"]
    gf3["off"] = gf2["off"]
    gf3["mu"] = gf2["mu"]
    gf3.thaw(["A", "sig", "off", "mu"])
    gf3.fit(intervallwave3, intervallflux3)
    gf3.parameterSummary()
    linienminimum3[i] = gf3["mu"]
    if linienminimum3[i] <= wellenlaenge - 4:
        break
    if linienminimum3[i] >= wellenlaenge + 4:
        break
    apex3[i] = gf3.evaluate(gf3["mu"])
    linienminimum3_bc[i] = linienminimum3[i] * (1 + corr[i] / 299710)
    RV3[i] = (linienminimum3_bc[i] - wellenlaenge) / \
        linienminimum3_bc[i] * 299710

    plt.plot(intervallwave3, gf3.model, "g--")
    plt.title(filelist[i] + " Beobachtungszeitpunkt " + str(obs_time[i]))
    plt.plot(linienminimum3[i], apex3[i], "o", color="black")
    plt.pause(.1)
    if frage_grafikenspeichern == 'j':
        plt.savefig(filelist[i].rsplit(".", 1)[0] + "_" + linie + ".png")
    RV_diff2[i] = RV3[i] - RV2[i]

    print('(baryzentrisch korrigierte) RV: ', RV3[i])


# # Plot der RV's
# fig=plt.figure()
# plt.plot(obs_time, RV1,'bo', markersize=1)
# plt.plot(obs_time, RV2,'r+', markersize=1)
# plt.plot(obs_time, RV3,'g+', markersize=1)

# Abspeichern als ascii-Datei
# ascii.write([obs_time, RV1], linie+'_RV1'+'.dat', overwrite=True,
#                 names=['JD', 'RV'], format='tab')

# ascii.write([obs_time, RV2], linie+'_RV2'+'.dat', overwrite=True,
#                 names=['JD', 'RV'], format='tab')

ascii.write(
    [filelist, obs_time, RV1, RV2, RV3],
    linie + "_RV_perGaussfit" + ".csv",
    overwrite=True,
    names=['Spectrum', "JD", "RV1", "RV2", "RV3"],
    format="csv",
)

# ascii.write([obs_time, apex3], linie+'_apex'+'.dat', overwrite=True,
#             names=['JD', 'apex'], format='tab')

# ascii.write([obs_time, linienminimum3_bc], linie+'_heliozentrisch_korrigierte_Linienminima' +
#             '.dat', overwrite=True, names=['JD', 'WAVE'], format='tab')

# # Speichern von obs_time und corr in ascii-file
# ascii.write([obs_time, corr], linie+'_Tabelle_obs_time_bc.dat', overwrite=True,
#             names=['JD', 'BARYCORR'], format='tab')

# # Speichern von obs_time und jhd in ascii-file
# ascii.write([obs_time, hjd], linie+'_Tabelle_obs_time_hjd.dat', overwrite=True,
#             names=['JD', 'HJD'], format='tab')

fr = input('Wenn Sie fertig mit der Betrachtung der Grafiken sind und das\
 Programm beenden möchten, drücken sie die Eingabe-Taste.')
plt.close('all')
print('\nProgramm ist beendet')
