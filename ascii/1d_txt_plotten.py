#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
1d_txt_plotten.py

Das Skript liest ein Spektrum, das in einer ASCII-Datei gespeichert ist. Die
Daten stehen in 2 Spalten mit den Überschriften WAVE und FLUX.
Pfad/Name werden abgefragt.
Das Spektrum wird lediglich geplottet.

Status 20180823
Author = Lothar Schanne
"""

import matplotlib.pyplot as plt
from astropy.io import ascii

spectrum_name = input('Enter the path and name of the text file: ')
spectrum = ascii.read(spectrum_name, guess=True)

# Plotting the spectrum
fig = plt.figure(figsize=(14, 10))
plt.plot(spectrum['WAVE'], spectrum['FLUX'])
plt.xlabel('Wavelength [Angström]')
plt.ylabel('ADU')
plt.title(spectrum_name)
plt.grid(True)

plt.pause(.1)
