#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ansehen eines wählbaren 1d-Spektrums im fits-Format,
Eingabe der Ionen/Elemente, deren Linien im Spektrum markiert werden sollen.
Die Grafik kann auf Wunsch als pdf abgespeichert werden.

Stand 20230905
author: Lothar Schanne
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

# lokaler Modul, muss im gleichen Verzeichnis stehen wie das Skript (oder das
# Verzeichnis des Moduls im Pythonpath, damit Python ihn findet)
# Der Modul befindet sich im Ordner Bausteine
import Linienlisten  # Elementlisten

plt.style.use('seaborn-whitegrid')
plt.switch_backend('Qt5Agg')  # Interaktives GrafikBackend einschalten


# Pfad und Name des Spektrumfiles bitte anpassen
file = input("Pfad und Filebezeichnung eingeben: ")

#   Lesen des Spektrums
sp = fits.open(file)

# Header lesen und in der Konsole ausdrucken
# HD = dict(sp[0].header)
# print("\n\nHeader des Spektrums :\n")
# print(HD)

if 'CRPIX1' not in sp[0].header:
    sp[0].header["CRPIX1"] = 1

# Erzeugen von Arrays mit den Wellenlängen und Fluxes des Spektrums
flux = np.array(sp[0].data)
wave = np.ones(sp[0].header["NAXIS1"], dtype=float)
for i in range(sp[0].header["NAXIS1"]):
    wave[i] = (
        sp[0].header["CRVAL1"]
        + (i + 1 - sp[0].header["CRPIX1"]) * sp[0].header["CDELT1"]
    )
# In der Liste wave sind die Wellenlängen der Pixel enthalten
# In der Liste flux die entsprechenden Intensitäten

#   Schliessen des fits-file
sp.close()

# Plot gesamtes Spektrum
fig = plt.figure(1, figsize=(14, 10))
plt.plot(wave, flux, "b-", linewidth=0.5)
plt.xlabel("Wellenlänge [Angström]", fontsize=14)
plt.ylabel("ADU", fontsize=14)
plt.title("Spektrum " + file, fontsize=18)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.grid(True)
plt.pause(.1)

# Auswahl der Linie oder Wellenlänge
element = Linienlisten.elementauswahl()
for i in element:
    plt.axvline(element[i], color='r')
    plt.text(element[i], 0.3, i)
    plt.pause(.1)

frage = input('Sollen weitere Linien markiert werden, dann y eingeben. ')
while frage == 'y':
    element = Linienlisten.elementauswahl()
    for i in element:
        plt.axvline(element[i])
        plt.text(element[i], 0.3, i)
        plt.pause(.1)
    frage = input('Sollen weitere Linien markiert werden, dann y eingeben. ')

frage2 = input('Soll die Grafik gespeichert werden? Dann y eingeben: ')
if frage2 == 'y':
    grafik = file.rsplit('.fit', 1)[0] + '.pdf'
    fig.savefig(grafik)

plt.close('all')
print('Ende des Programms')
