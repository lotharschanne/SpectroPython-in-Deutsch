#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Die normierten Ordnungen eines Echellespektrums (im fits-Format) werden
aneinander gehängt (mergen). Die Ordnungen müssen so geordnet sein, dass die
Wellenlängen monoton aufsteigend sind.

20241001
@author: lothar
"""

import numpy as np
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import glob


def wavecalc(flux, header):
    """
    Berechnet die Wellenlängen aus den Headerdaten der fits-Datei
    """
    wave = np.zeros(header["NAXIS1"])
    if "CRPIX1" not in header:
        header["CRPIX1"] = 1
    header["CRVAL1"] = header["CRVAL1"] + \
        (1 - header["CRPIX1"]) * header["CDELT1"]
    for i in np.arange(header["NAXIS1"]):
        wave[i] = header["CRVAL1"] + i * header["CDELT1"]
    return wave, flux


plt.switch_backend('Qt5Agg')  # Interaktives GrafikBackend einschalten

# Create file list for objects. All spectra in one (sub)folder.
files = input('Pfad und Name der normierten Ordnungen (use wildcards) : ')
filelist_obj = glob.glob(files)
filelist_obj.sort()

# Print the list
# print('\nList of spectra: \n')
print('Number of spectra: ', len(filelist_obj), '\n')

# Überprüfen auf lückenlosen Anschluss der Ordnungen
crval1_obj = np.zeros(len(filelist_obj))
cdelt1_obj = np.zeros(len(filelist_obj))
for i in range(len(filelist_obj)):
    header = fits.getheader(filelist_obj[i])
    crval1_obj[i] = header['CRVAL1']
    cdelt1_obj[i] = header['CDELT1']

for i in range(len(crval1_obj)-1):
    header1 = fits.getheader(filelist_obj[i])
    header2 = fits.getheader(filelist_obj[i+1])
    if cdelt1_obj[i] == cdelt1_obj[i+1] and (header2['CRVAL1'] - header1['CRVAL1']) % header1['CDELT1'] < 1e-6:
        print(i, 'OK, Ordnungen passen aufeinander')
    else:
        print(filelist_obj[i], ' Ordnungen nicht passend')

print()

# Einlesen der ersten (blauseitigen) Ordnung
flux_merged, header_obj = fits.getdata(filelist_obj[0], header=True)
# Generate the spectrum sections and save them as fit:
wave_merged, flux_merged = wavecalc(flux_merged, header_obj)

print('Spektrum: ', filelist_obj[0])
print(f'Objekt Wellenlänge des ersten Pixels: {wave_merged[0]:.2f}')
print(f'Objekt Wellenlänge des letzten Pixels: {wave_merged[-1]:.2f}')
print('Anzahl der Pixel: ', header_obj['NAXIS1'])
print('Schrittweite', header_obj['CDELT1'])
print()

letztesPixel = wave_merged[-1]

# Anhängen der weiteren Ordnungen
for i in range(1, len(filelist_obj)):
    flux, header_obj = fits.getdata(filelist_obj[i], header=True)
    wave_obj, flux_obj = wavecalc(flux, header_obj)
    print('Spektrum: ', filelist_obj[i])
    print(f'Objekt Wellenlänge des ersten Pixels: {wave_obj[0]:.2f}')
    print(f'Objekt Wellenlänge des letzten Pixels: {wave_obj[-1]:.2f}')
    print('Anzahl der Pixel: ', header_obj['NAXIS1'])
    print('Schrittweite', header_obj['CDELT1'])
    print()

    if letztesPixel >= wave_obj[0]:  # im Falle der Überlappung von Ordnungen
        for k in np.arange(len(wave_merged)-1, -1, -1):
            if (wave_obj[0] >= wave_merged[k]):
                break
        for m in np.arange(0, len(wave_obj)):
            if (wave_obj[m] == wave_merged[-1]) or wave_obj[m] > wave_merged[-1]:
                break

        ueberlapp_flux = np.zeros(m)
        for s in range(m):
            ueberlapp_flux[s] = (flux_merged[k+s] + flux_obj[s]) / 2

        wave_merged_red = wave_merged[0:k]
        flux_merged_red = flux_merged[0:k]
        wave_merged_red = np.hstack((wave_merged_red, wave_obj[0:m]))
        flux_merged_red = np.hstack((flux_merged_red, ueberlapp_flux))
        wave_merged = np.hstack((wave_merged_red, wave_obj[m:]))
        flux_merged = np.hstack((flux_merged_red, flux_obj[m:]))
        letztesPixel = wave_merged[-1]

    else:  # im Falle einer Lücke zwischen Ordnungen (keine Überlappung)
        delta = wave_obj[0] - letztesPixel
        pixellücke = int(delta / header_obj['CDELT1'])
        lückenfüller = np.linspace(letztesPixel + header_obj['CDELT1'],
                                   wave_obj[0],
                                   num=int(delta / header_obj['CDELT1']),
                                   endpoint=True)
        fluxfüller = np.zeros(len(lückenfüller))
        wave_merged = np.hstack((wave_merged, lückenfüller))
        flux_merged = np.hstack((flux_merged, fluxfüller))
        wave_merged = np.hstack((wave_merged, wave_obj))
        flux_merged = np.hstack((flux_merged, flux_obj))
        letztesPixel = wave_merged[-1]


plt.figure(figsize=(10, 8))
plt.plot(wave_merged, flux_merged, linewidth=.5)
plt.grid(True)
plt.xlabel('Wellenlänge [Angström]')
plt.ylabel('auf das Kontinuum normierte Intensität')
plt.savefig('gemerged.pdf')

# Abspeichern des germergten Spektrums
header_obj["CRVAL1"] = wave_merged[0]
header_obj["NAXIS1"] = len(wave_merged)
header_obj["CRPIX1"] = 1
name = "gemerged.fits"

fits.writeto(
    name, flux_merged, header_obj, overwrite=True, output_verify="silentfix"
)

ascii.write([wave_merged, flux_merged], "gemerged.dat", overwrite=True,
            names=['WAVE', 'FLUX'], format='tab')
