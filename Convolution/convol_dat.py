#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Einlesen eines  1d-Spektrums in Form einer ascii-Tabelle mit der extension .dat,
2-spaltig mit den Spaltenüberschriften 'WAVE' und 'FLUX'.
Das Spektrum kann mit einer Standardabweichung convolviert werden
(stddev in der Einheit des step).
Der convolvierte Flux wird zusammen mit der Wellenlänge in einer zweispaltigen
ascii-Tabelle abgespeichert (spacer = tab). Mit den Spaltenüberschriften
'WAVE' und 'FLUX'.

Stand 20180824
@author: Lothar Schanne
"""
# from __future__ import print_function, division
import matplotlib.pyplot as plt
from astropy.convolution import Gaussian1DKernel
from astropy.convolution import convolve
from astropy.io import ascii

spectrum_name = input('Pfad und Filebezeichnung eingeben: ')

std = float(input('Standardabweichung für die Convolution in Vielfachen der\
 Schrittweite: '))

# Spektrum einlesen
spectrum = ascii.read(spectrum_name, guess=True)

#   Convolve mit astropy.convolution
kernel = Gaussian1DKernel(stddev=std)
convoluted = convolve(spectrum['FLUX'], kernel, normalize_kernel=True,
                      boundary='extend')

name = spectrum_name+'_convolved.dat'

# Grafik
fig = plt.figure(1, figsize=(14, 10))
plt.suptitle('Spektrum '+spectrum_name)
plt.subplot(2, 1, 1)
plt.plot(spectrum['WAVE'], spectrum['FLUX'])
plt.xlabel('Wellenlänge in Angström')
plt.ylabel('relative Flux')
plt.subplot(2, 1, 2)
plt.plot(spectrum['WAVE'], convoluted)
plt.xlabel('Wellenlänge in Angström')
plt.ylabel('relative Flux')
plt.title(name)
# speichern
plt.savefig(name+'.png')
plt.savefig(name+'.pdf')
plt.show(block=True)

# Speichern des convolved
ascii.write([spectrum['WAVE'], convoluted], name, overwrite=True,
            names=['WAVE', 'FLUX'], format='tab')
