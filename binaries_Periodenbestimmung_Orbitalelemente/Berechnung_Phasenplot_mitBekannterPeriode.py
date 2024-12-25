#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Das Skript liest eine ascii-Datei ein mit den beiden Spalten (ohne Spalten
überschrift !) JD (Beobachtungszeitpunkte) und RV (Radialgeschwindigkeiten).
Die Zeitpunkte (JD) werden mit der einzugebenden (als bekannt vorausgesetzten)
Periode gefaltet und die zugehörigen RV's mit den Phasen geplottet.
Der Plot wird als pdf gespeichert.
Die Ergebnisse werden in einer csv-Datei abgespeichert mit den Spaltenüberschriften
JD, RV und Phase.

Created on Mon Nov  6 13:58:44 2023

@author: Lothar Schanne
"""

from astropy.io import ascii
import matplotlib.pyplot as plt
from PyAstronomy.pyasl import foldAt
import matplotlib.pylab as plt
import numpy as np


plt.switch_backend('Qt5Agg')  # Interaktives GrafikBackend einschalten

filename = input("Geben Sie den Name des Datenfiles ein: ")
ts = ascii.read(filename)
ts_JD = ts.columns[0]
ts_RV = ts.columns[1]

Periode = float(input('Geben Sie die Periode in Tagen ein: '))

phases = foldAt(ts_JD, Periode, ts_JD[0])

# Sort with respect to phase
# First, get the order of indices ...
sortIndi = np.argsort(phases)
# ... and, second, rearrange the arrays.
phases = phases[sortIndi]
RV = ts_RV[sortIndi]
JD = ts_JD[sortIndi]

fig = plt.figure()
plt.plot(phases*Periode, RV, 'ko', markersize=4)
plt.xlabel('Zeit [Tage]')
plt.ylabel('RV [km/s]')
plt.title('Phasenplot ' + filename + ', mit Periode ' + str(Periode) +
          ' berechnet')
plt.savefig(filename + '_Phasenplot_Periode_' + str(Periode) + '.pdf', format='pdf')

# plt.show(block=True)
plt.pause(.1)

# Abspeichern der Ergebnisse als ascii-Datei (im csv-Format)
ascii.write(
    [JD, RV, phases],
    filename + "_Phasen_Periode" + str(Periode) + "d.csv",
    overwrite=True,
    names=[
        "      JD     ",
        "      RV      ",
        '  Phase  ',
    ],
    format="csv",
)
