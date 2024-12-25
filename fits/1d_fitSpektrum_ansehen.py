#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
1d_Spektrum_ansehen.py

Ansehen eines 1d-Spektrums im fits-Format,
Auslesen und Anzeige der Headerdaten.
Plotten des Spektrums.

Version 20230527
@author: lothar schanne
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import glob

plt.switch_backend('Qt5Agg')  # Interaktives GrafikBackend einschalten

# Pfad und Name des Spektrumfiles bitte anpassen
file = input("Pfad und Filebezeichnung eingeben, wildcards benutzbar: ")
file = glob.glob(file)
file = file[0]

#   Lesen des Spektrums
sp = fits.open(file)

# Header lesen und in der Konsole ausdrucken
HD = dict(sp[0].header)
print("\n\nHeader des Spektrums :\n")
for i in HD:
    print(i, ':', HD[i])

try:
    sp[0].header["CRPIX1"]
except:
    sp[0].header["CRPIX1"] = 1

# Erzeugen von Arrays mit den Wellenlängen und Fluxes des Spektrums
flux = np.array(sp[0].data)
wave = np.ones(sp[0].header["NAXIS1"], dtype=float)
for i in range(sp[0].header["NAXIS1"]):
    wave[i] = (
        sp[0].header["CRVAL1"]
        + (i + 1 - sp[0].header["CRPIX1"]) * sp[0].header["CDELT1"]
    )

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

# Zeigen des Plots, wenn das Skript in einer normalen Python-Konsole
# durchgeführt wird. Für die IPython-Konsole in Spyder ist das nicht nötig.
# plt.show(block=True)
