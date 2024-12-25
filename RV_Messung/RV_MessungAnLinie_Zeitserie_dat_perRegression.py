#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Liest Spektrenkatalog ein (Zeitserie von normierten 1d-Spektren im
tab-Format). Fitted die gewählte Linie per Regression nten Grades
und bestimmt aus dem Minimum die Radialgeschwindigkeit RV.
Plottet alle fittings und gibt ermittelte Daten als ascii-Dateien
(Komma-separiert, als .csv) aus.

Stand 20221105

@author: lothar schanne
"""

import numpy as np
from astropy.io import ascii
import glob
from PyAstronomy import pyasl
import matplotlib.pyplot as plt
from astroplan import FixedTarget
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroplan import FixedTarget
from astropy.coordinates import SkyCoord
import astropy.units as u

# lokaler Modul, im Ordner "Bausteine" zu finden,
# muss im gleichen Verzeichnis stehen wie das Skript (oder das
# Verzeichnis des Moduls im Pythonpath, damit Python ihn findet)
import Linienlisten

plt.style.use('seaborn-whitegrid')
plt.switch_backend('Qt5Agg')  # Interaktives GrafikBackend einschalten


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

frage_grafikenspeichern = input(
    'Wenn Sie die Grafiken speichern möchten, geben Sie j ein: ')


# # Abarbeiten der filelist, Einlesen von flux und header, Auswahl des Flux um die Linie:

# Definition von Variablen
RV = np.zeros(len(filelist))
apex = np.zeros(len(filelist))
linien_minwave = np.zeros(len(filelist))
miniwave = np.zeros(len(filelist))
miniflux = np.zeros(len(filelist))


for i in range(len(filelist)):
    ts = ascii.read(filelist[i], format="tab")

    # **********************************************************
    # Breite des Suchintervalls, Breite in Angström, anpassen.
    suchintervall = 3
    # ************************************************************

    intervall = np.extract(
        abs(ts.columns[0] - wellenlaenge) < suchintervall, ts)
    a = np.zeros(len(intervall))
    b = np.zeros(len(intervall))

    for j in range(len(intervall)):
        a[j], b[j] = intervall[j]
    # apex[i] = np.min(b)
    linien_argmin = np.argmin(b)
    # absolutes Linienminimum (pixelweise)
    linien_minwave[i] = a[linien_argmin]

    fig = plt.figure(i)
    plt.plot(a, b, "-", linewidth=0.5)
    plt.title(filelist[i].rsplit(".", 1)[0] + "_" + linie)

    # Das folgende Intervall [Angström] anpassen, bei verrauschten Linien breit,
    # bei nicht verrauschten, schmalen Linien klein.
    linienminimum_wave = a[abs(a - linien_minwave[i]) <= 0.5]
    linienminimum_flux = b[abs(a - linien_minwave[i]) <= 0.5]

    # Regression, Grad anpassen (2 oder 4 oder 6)
    grad = 4
    model = np.poly1d(np.polyfit(linienminimum_wave, linienminimum_flux, grad))

    polyline = np.linspace(linienminimum_wave[0], linienminimum_wave[-1], 100)
    modflux = model(polyline)

    # Linienminimum subpixelgenau berechnen
    miniflux[i] = modflux.min()
    miniwave[i] = polyline[modflux.argmin()]

    # Plotten
    plt.plot(polyline, modflux, "--")
    plt.plot(miniwave[i], miniflux[i], "o", color="black")
    plt.pause(.1)
    if frage_grafikenspeichern == 'j':
        plt.savefig(filelist[i].rsplit(".", 1)[0] + "_" + linie + ".png")

    # RV ohne baryzentrische Korrektur
    RV[i] = (miniwave[i] - wellenlaenge) / wellenlaenge * 299710

# # # Plot der RV's
# # fig = plt.figure()
# # plt.plot(obs_time, RV, 'bo', markersize=1)

# Abspeichern als ascii-Datei
ascii.write(
    [filelist, RV, miniflux],
    linie + "_RV_Regression_grad" + str(grad) + ".csv",
    overwrite=True,
    names=["Spektrum               ",
           "        RV          ", "         miniFlux", ],
    format="csv",
)

fr = input('Wenn Sie fertig mit der Betrachtung der Grafiken sind und das\
 Programm beenden möchten, drücken sie die Eingabe-Taste.')
plt.close('all')
print('\nProgramm ist beendet')
