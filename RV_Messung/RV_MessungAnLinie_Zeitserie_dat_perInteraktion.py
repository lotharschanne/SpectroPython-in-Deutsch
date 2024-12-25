#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Liest einen Spektrenkatalog ein (Zeitserie von normierten 1d-Spektren im
tab-Format, 2 Spalten WAVE und FLUX). Festlegung des Linienminimums per
Interaktion und Bestimmung der Radialgeschwindigkeit aus dem Minimum.
Plottet alle Spektren zum markieren von bis zu 2 Linien-Minima (im Falle eines
Doppelsterns SB2) und gibt die ermittelten Daten (RV und Apex) als ascii-Dateien
(Komma-separiert, als .csv) aus. Die RV's sind nicht baryzentrisch korrigiert.'

Stand 20230610

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


# Fileliste erstellen. Spektren in einem (Unter)Ordner.
files = input("Pfad und Name der Spektren (nutze wildcards) : ")
filelist = glob.glob(files)

# Alphabetisches Sortieren. Bei richtiger Namensgebung ergibt das eine
# zeitliche Ordnung
filelist.sort()

# Ausdruck der Liste
print("\nSpektrenliste: \n")
print("Anzahl der Spektren: ", len(filelist), "\n")

# Auswahl der Linie bzw. Wellenlänge
linie, wellenlaenge = Linienlisten.linienauswahl()


# Abarbeiten der filelist, Einlesen von flux und Wellenlängen, Auswahl des
# Flux um die Linie:

# Definition von Variablen
RV1 = np.zeros(len(filelist))
RV2 = np.zeros(len(filelist))
apex1 = np.zeros(len(filelist))
apex2 = np.zeros(len(filelist))
linien_minwave1 = np.zeros(len(filelist))
miniwave1 = np.zeros(len(filelist))
miniwave2 = np.zeros(len(filelist))
miniflux1 = np.zeros(len(filelist))
miniflux2 = np.zeros(len(filelist))


for i in range(len(filelist)):
    ts = ascii.read(filelist[i], format="csv")

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
    plt.title(filelist[i].rstrip(".dat") + "_" + linie)

    print("\nSpektrum ", filelist[i])
    print("Klicke das Minimum der ersten Linie an: ")
    pts = np.asarray(plt.ginput(n=1, timeout=-1))
    RV1[i] = (pts[0, 0] - wellenlaenge) / wellenlaenge * 299710
    apex1[i] = pts[0, 1]
    print("Erste Linie, RV =", RV1[i], "Apex =", apex1[i])

    frage = input(
        "Möchten Sie das Minimum einer zweiten Linie anklicken? Dann y eingeben: "
    )
    if frage == "y":
        print("Klicke die zweite Linie an")
        pts = np.asarray(plt.ginput(n=1, timeout=-1))
        RV2[i] = (pts[0, 0] - wellenlaenge) / wellenlaenge * 299710
        apex2[i] = pts[0, 1]
        print("Zweite Linie, RV =", RV2[i], "Apex =", apex2[i])
    else:
        RV2[i] = np.NaN
        apex2[i] = np.NaN

    plt.savefig(filelist[i].rsplit(".", 1)[0] + "_" + linie + ".png")
    plt.close()


# Abspeichern als ascii-Datei
ascii.write(
    [filelist, RV1, apex1, RV2, apex2],
    linie + "_RV_interaktiv" + ".csv",
    overwrite=True,
    names=["Spektrum", "RV1", "Apex1", "RV2", "Apex2"],
    format="csv",
)

fr = input('Wenn Sie fertig mit der Betrachtung der Grafiken sind und das\
 Programm beenden möchten, drücken sie die Eingabe-Taste.')
plt.close('all')
print('\nProgramm ist beendet')
