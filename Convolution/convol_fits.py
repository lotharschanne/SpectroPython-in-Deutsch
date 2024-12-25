#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Einlesen eines  1d-Spektrums in Form einer fits-Datei.

Das Spektrum kann mit einer Standardabweichung, ausgedrückt durch die
Auflösung R, convolviert werden. Die convolvierte fits-Datei wird abgespeichert.

Stand 20220105
@author: Lothar Schanne
"""
# from __future__ import print_function, division
import matplotlib.pyplot as plt
from astropy.convolution import Gaussian1DKernel
from astropy.convolution import convolve
from astropy.io import ascii, fits

spectrum_name = input('Pfad und Filebezeichnung eingeben: ')

R = float(input('Gewünschtes R: '))

# Spektrum einlesen
sp = fits.open(spectrum_name)
# print('\n\nHeader of the spectrum :\n\n', sp[0].header, '\n\n')
std = sp[0].header['CRVAL1']/R/sp[0].header['CDELT1']

#   Convolve mit astropy.convolution
kernel = Gaussian1DKernel(stddev=std)
sp[0].data = convolve(sp[0].data, kernel, normalize_kernel=True,
                      boundary='extend')

name = spectrum_name.rstrip('.fits')+'_convolved_R'+str(R)+'.fits'


# Speichern des convolved
fits.writeto(name, sp[0].data, sp[0].header,
             overwrite=True, output_verify='silentfix')
