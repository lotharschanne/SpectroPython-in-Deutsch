#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Das Skript dient dazu ein Spektrum mit zusätzlichem gausschen Rauschen 
zu versehen.

Created on Mon Dec  6 09:55:36 2021

@author: lothar
"""

from PyAstronomy import pyasl
import numpy as np
from astropy.io import fits, ascii
from random import gauss


file = input('Pfad und Name des Spektrums (nutze wildcards) : ')


flux, header = fits.getdata(file, header=True)

# Mit dem sigma s der Gaussfunktion wird die Verteilung des Rauschens eingestellt
s = 0.2
for i in range(len(flux)):
    fehler = gauss(0, s)
    flux[i] = flux[i] + fehler


filename = file + '_zusätzlichesRauschen_' + str(s)+'.fits'
fits.writeto(filename, flux,
             header, overwrite=True)
