#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
comparison_plot_2spectra_ascii.py

Overplot von zwei 1d-Spektren im ASCII-Format (.dat), mit 2 Spalten
WAVE und FLUX mit Fließkommazahlen. Der Plot wird gespeichert als PDF.

State 20180824
@author: Lothar Schanne
"""

import matplotlib.pyplot as plt
from astropy.io import ascii

# Input of asci-files
template_name = input('Geben Sie den Namen des ersten Spektrums ein: ')
spectrum_name = input('Geben Sie den Namen des zweiten Spektrums ein: ')

template = ascii.read(template_name, format='tab')
spectrum = ascii.read(spectrum_name, format='tab')

fig = plt.figure(figsize=(14, 7))
plt.xlim(6490, 6605)  # please adjust !!!!!
plt.plot(template['WAVE'], template['FLUX'])
plt.plot(spectrum['WAVE'], spectrum['FLUX'])
plt.xlabel('Wavelength [Angstroem]')
plt.ylabel('relative Flux')
plt.title(spectrum_name.rstrip('.dat'))

plt.legend([template_name.rstrip('.dat'), spectrum_name.rstrip('.dat')],
           loc='lower left')

plt.pause(.1)

# fig.savefig(spectrum_name.rstrip('.dat')+'overplot.png')
fig.savefig(spectrum_name.rstrip('.dat')+'overplot.pdf')

print('Ende des Programms')
