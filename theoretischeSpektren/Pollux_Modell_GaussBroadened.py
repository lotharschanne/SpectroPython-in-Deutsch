#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Konvolution mit einer Gaußfunktion auf Spektrographenauflösung,
Auswahl des Wellenlängenbereiches und
Anpassung auf eine gewünschte Schrittbreite (Angström/Pixel).
Speicherung als ascii-tab-Tabelle und als .fit mit dem einzugebenden Filename.
Speicherung des plots als pdf.

Created on Fri Apr  2 13:53:54 2021

@author: lothar
"""
from PyAstronomy import pyasl
import matplotlib.pyplot as plt
from astropy.io import ascii, fits
import numpy as np


plt.switch_backend('Qt5Agg')  # Interaktives GrafikBackend einschalten

file = input("Pfad und Filebezeichnung des Polluxspektrums eingeben: ")
table = ascii.read(file)

specname = input('Filebezeichnung für die Speicherung des erzeugten Spektrums eingeben: ')+'_'

print('\nEingabe des gewünschten Wellenlängenbereichs in Angström:')
a = float(input('Begin: '))
b = float(input('End: '))

newtable = table[table['col1'] >= a]
newtable = newtable[newtable['col1'] <= b]
wave = newtable['col1']
flux = newtable['col3']

print('Schrittweite des Modellspektrums = ', wave[1]-wave[0])

R = float(input('Gewünschte relative Auflösung R eingeben (z.B. 10000): '))
newstep = float(
    input('Geben Sie die gewünschte Schrittweite in Angström des Endergebnisses ein: '))
spectrum_name = specname + '_R' + str(int(R))

convolutedFlux, fwhm = pyasl.instrBroadGaussFast(
    wave, flux, R, edgeHandling='firstlast', fullout=True, maxsig=5.0, equid=True)
print('FWHM des gaußverbreiterten Spektrums = ', fwhm)

neue_wellenlaengen = np.arange(wave[0], wave[-1], newstep)
newwave, newflux = pyasl.equidistantInterpolation(wave, convolutedFlux, neue_wellenlaengen)

name = spectrum_name + '_convolved.dat'
ascii.write([newwave, newflux], name, overwrite=True,
            names=['WAVE', 'FLUX'], format='tab')
hdu = fits.PrimaryHDU()
hdu.header['CRVAL1'] = newwave[0]
hdu.header['CRPIX1'] = 1
hdu.header['NAXIS1'] = len(newwave)
hdu.header['CDELT1'] = newstep
hdu.header['NAXIS'] = 1
name = spectrum_name + '_convolved.fit'
fits.writeto(name, newflux, hdu.header,
              overwrite=True, output_verify='silentfix')


# Grafik
fig, ax = plt.subplots(nrows=2, ncols=1, sharex=True)
fig.suptitle('Spektrum: ' + file, fontsize=10)

ax[0].plot(wave, flux, linewidth=.5)
ax[0].set(ylabel='relativer Flux')
ax[0].set_title('Pollux-Modell', fontsize=8)

ax[1].plot(newwave, newflux, linewidth=.5)
ax[1].set(xlabel='Wellenlänge in Angström', ylabel='relativer Flux')
ax[1].set_title('convolved auf R ' + str(R), fontsize=8)
# ax[1].set_xlim(6500,6700)
plt.pause(1)

plt.savefig(spectrum_name + '_convolved.png', format='png')
# plt.savefig(spectrum_name + '_convolved.pdf', format='pdf')
