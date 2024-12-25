#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Liest Spektrenkatalog ein (Zeitserie von normierten 1d-Spektren im
tab-Format). Festlegung des ungefähren Linienminimums per Interaktion
und bestimmt aus dem per Regression ermittelten Minimum die Radialgeschwindigkeit
RV. Auf Wunsch kann eine zweite Linie behandelt werden.
Plottet alle fittings und gibt ermittelte Daten als ascii-Dateien
(Komma-separiert, als .csv) aus.
Die RV's sind nicht baryzentrisch korrigiert.

Stand 20220408

@author: lothar schanne
"""

import numpy as np
from astropy.io import ascii
import glob
from PyAstronomy import pyasl
import matplotlib.pyplot as plt

# lokaler Modul, muss im gleichen Verzeichnis stehen wie das Skript (oder das
# Verzeichnis des Moduls im Pythonpath, damit Python ihn findet)
import Linienlisten

plt.switch_backend('Qt5Agg')  # Interaktives GrafikBackend einschalten
plt.style.use('seaborn-whitegrid')


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


# Abarbeiten der filelist, Einlesen von flux und Wellenlängen,
# Auswahl des Flux um die Linie:

# Definition von Variablen
RV1 = np.zeros(len(filelist))
RV2 = np.zeros(len(filelist))
miniwave1 = np.zeros(len(filelist))
miniwave2 = np.zeros(len(filelist))
miniflux1 = np.zeros(len(filelist))
miniflux2 = np.zeros(len(filelist))


for i in range(len(filelist)):
    ts = ascii.read(filelist[i], format="tab")

    # **********************************************************
    # Breite des Suchintervalls, Breite in Angström, anpassen.
    suchintervall = 4
    # ************************************************************

    intervall = np.extract(
        abs(ts.columns[0] - wellenlaenge) < suchintervall, ts)
    a = np.zeros(len(intervall))
    b = np.zeros(len(intervall))
    for j in range(len(intervall)):
        a[j], b[j] = intervall[j]

    fig = plt.figure(i)
    plt.plot(a, b, "-", linewidth=0.5)
    plt.title(filelist[i].rsplit(".", 1)[0] + "_" + linie)

    print("\nSpektrum ", filelist[i])
    print("Klicke das Minimum der ersten Linie an: ")
    pts = np.asarray(plt.ginput(n=1, timeout=-1))

    # Das folgende Intervall [Angström] anpassen, bei verrauschten Linien breit,
    # bei nicht verrauschten, schmalen Linien klein.
    linienminimum_wave1 = a[abs(a - pts[0, 0]) <= 0.5]
    linienminimum_flux1 = b[abs(a - pts[0, 0]) <= 0.5]

    ################ Regression, Grad anpassen (2 oder 4 oder 6) ###########
    grad = 4
    ########################################################################
    model = np.poly1d(np.polyfit(
        linienminimum_wave1, linienminimum_flux1, grad))

    polyline = np.linspace(
        linienminimum_wave1[0], linienminimum_wave1[-1], 100)
    modflux = model(polyline)

    # Linienminimum berechnen
    miniflux1[i] = modflux.min()
    miniwave1[i] = polyline[modflux.argmin()]

    # Plotten
    plt.plot(polyline, modflux, "--")
    plt.plot(miniwave1[i], miniflux1[i], "o", color="black")
    # plt.show()
    # plt.savefig(filelist[i].rsplit(".", 1)[0] + "_" + linie + ".png")

    RV1[i] = (miniwave1[i] - wellenlaenge) / wellenlaenge * 299710

    frage = input(
        "Möchten Sie das Minimum einer zweiten Linie anklicken? Dann y eingeben: "
    )
    if frage == "y":
        print("Klicke die zweite Linie an")
        pts = np.asarray(plt.ginput(n=1, timeout=-1))

        # Das folgende Intervall [Angström] anpassen, bei verrauschten Linien breit,
        # bei nicht verrauschten, schmalen Linien klein.
        linienminimum_wave2 = a[abs(a - pts[0, 0]) <= 0.3]
        linienminimum_flux2 = b[abs(a - pts[0, 0]) <= 0.3]

        # Regression
        model = np.poly1d(np.polyfit(
            linienminimum_wave2, linienminimum_flux2, grad))

        polyline = np.linspace(
            linienminimum_wave2[0], linienminimum_wave2[-1], 100)
        modflux = model(polyline)

        # Linienminimum berechnen
        miniflux2[i] = modflux.min()
        miniwave2[i] = polyline[modflux.argmin()]

        # Plotten
        plt.plot(polyline, modflux, "--")
        plt.plot(miniwave2[i], miniflux2[i], "o", color="black")

        RV2[i] = (miniwave2[i] - wellenlaenge) / wellenlaenge * 299710
    else:
        RV2[i] = np.NaN

    plt.pause(.1)
    if frage_grafikenspeichern == 'j':
        plt.savefig(filelist[i].rsplit(".", 1)[0] + "_" + linie + ".png")


# Abspeichern als ascii-Datei
ascii.write(
    [filelist, RV1, RV2],
    linie + "_RV_interaktiv2lines" + ".csv",
    overwrite=True,
    names=["Spektrum", "RV1", "RV2"],
    format="csv",
)

fr = input('Wenn Sie fertig mit der Betrachtung der Grafiken sind und das\
 Programm beenden möchten, drücken sie die Eingabe-Taste.')
plt.close('all')
print('\nProgramm ist beendet')
