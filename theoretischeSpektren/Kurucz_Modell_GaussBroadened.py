#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Laden eines Kuruzc-Sternatmosphärenmodells (1d-Spektrum) nach Auswahl von Teff
und logg (Internetverbindung erforderlich).
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


plt.switch_backend('Qt5Agg')  # Interaktives GrafikBackend einschalten

# Auswahl des Modells mittels Teff, logg und Metallizität
model = pyasl.resBased.spectralLib.SpectralLib(lib='A')
print('Liste der zur Verfügung stehenden Modelle (Parameter und Filenamen)')
model.listInventory()

Teff = float(input('Geben Sie Teff ein: '))
logg = float(input('Geben Sie logg ein: '))
met = 0.0  # momentan gibt es nur Modelle mit Metallicity = 0.0
# met = float(input('Geben Sie die Metallizität ein: '))

Model = model.requestModel(Teff, logg, met)
wave, flux = model.read1dFitsSpec(Model)

specname = input('Filebezeichnung für Speicherung eingeben: ')+'_'

print('Schrittweite des Modellspektrums = ', wave[1]-wave[0])
R = float(input('Gewünschte relative Auflösung R eingeben (z.B. 10000): '))
newstep = float(
    input('Geben Sie die gewünschte Schrittweite des Endergebnisses ein: '))
spectrum_name = specname + str(Teff) + '_' + str(logg) + '_R' + str(int(R))

print('\nEingabe des gewünschten Wellenlängenbereichs in Angström:')
a = float(input('Begin: '))
b = float(input('End: '))

for m in range(len(wave)):
    if wave[m] > a:
        aind = m
        break
for m in range(len(wave)):
    if wave[m] > b:
        bind = m
        break
wave = wave[aind:bind]
flux = flux[aind:bind]


convolutedFlux, fwhm = pyasl.instrBroadGaussFast(
    wave, flux, R, edgeHandling='firstlast', fullout=True, maxsig=5.0)
print('FWHM des gaußverbreiterten Spektrums = ', fwhm)

data, dt = pyasl.binningx0dt(
    wave, convolutedFlux, x0=wave[0], dt=newstep, useBinCenter=True, useMeanX=True)
name = spectrum_name + '_convolved.dat'
ascii.write([data[:, 0], data[:, 1]], name, overwrite=True,
            names=['WAVE', 'FLUX'], format='tab')
hdu = fits.PrimaryHDU()
hdu.header['CRVAL1'] = data[0, 0]
hdu.header['CRPIX1'] = 1
hdu.header['NAXIS1'] = len(data)
hdu.header['CDELT1'] = newstep
hdu.header['NAXIS'] = 1
name = spectrum_name + '_convolved.fit'
fits.writeto(name, data[:, 1], hdu.header,
             overwrite=True, output_verify='silentfix')


# Grafik
fig, ax = plt.subplots(nrows=2, ncols=1)
fig.suptitle('Spektrum: ' + spectrum_name, fontsize=12)

ax[0].plot(wave, flux, linewidth=.3)
ax[0].set(ylabel='relative Flux')
ax[0].set_title('Kurucz-Modell', fontsize=8)

ax[1].plot(data[:, 0], data[:, 1], linewidth=.3)
ax[1].set(xlabel='Wellenlänge in Angström', ylabel='relative Flux')
ax[1].set_title('convolved', fontsize=8)
# ax[1].set_xlim(6500,6700)

plt.savefig(spectrum_name + '_convolved.pdf')
