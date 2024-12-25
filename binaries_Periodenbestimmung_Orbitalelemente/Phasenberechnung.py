#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Das Skript liest einen csv-File mit den Spalten JD und RV (ohne Spalten-
überschriften) ein. Mit den im Skript vorgegebenen Parametern Periode und T0
werden dann die Beobachtunsgszeitpunkte gefaltet, der Phasenplot angezeigt und
die Ergebnisse in einem ascii-File (Format tab, Spalten 'JD', 'Phase', 'RV') 
gespeichert.
Periode und T0 anpassen (Zeilen 22 und 23).

Created on Tue Nov 14 19:10:56 2023

@author: lothar schanne
"""
from PyAstronomy import pyasl
from astropy.io import ascii
from astropy.table import Table
import matplotlib.pyplot as plt

plt.switch_backend('Qt5Agg')  # Interaktives GrafikBackend einschalten

Periode = 2.867328  # Bitte anpassen, Einheit Tage
T0 = 2441773.49     # Bitte anpassen, Einheit JD (Tage)

filename = input("Geben Sie den Name des Datenfiles ein: ")
ts = ascii.read(filename)

jd = ts.columns[0]
rv = ts.columns[1]


# Folding (Phase)
phases = pyasl.foldAt(jd, Periode, T0)
fig = plt.figure()
plt.plot(phases, rv, 'ko', markersize=3)
plt.xlabel('Phase')
plt.ylabel('RV [km/s]')
plt.title(filename.rsplit('.', 1)[0] + 'Phasendiagramm mit der Periode ' +
          str(Periode)+' d berechnet')
# plt.text(0, 0, 'Periode '+str(Periode.round(5)))
plt.savefig(filename + '_Phasenplot', format='pdf')

# plt.show(block=True)
plt.pause(.1)

data = Table([jd, phases, rv], names=['JD', 'Phase', 'RV'])
ascii.write(data, filename.rsplit('.', 1)[0] + "_Phasen.dat",
            overwrite=True, format="tab")
