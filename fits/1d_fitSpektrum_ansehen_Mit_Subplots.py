#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
1d_Spektrum_ansehen_Mit_Subplots.py

Ansehen eines wellenlängenkalibrierten 1d-Spektrums, Auslesen und Anzeige der 
Headerdaten. Plotten des gesamten Spektrums und einer Aufteilung in eine 
wählbare Zahl von Ausschnitten, die dann in entsprechenden Unterplots 
dargestellt werden. Zusätzlich kann ein beliebiger Wellenlängenbereich 
ausgewählt und geplottet werden.

Created on Sunday Jul 22 2018

@author: lothar schanne
"""


import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

plt.switch_backend('Qt5Agg')  # Interaktives GrafikBackend einschalten

# filelist erstellen
# Pfad und Name der Spektren, files bitte anpassen
file = input('Pfad und Filebezeichnung eingeben: ')

#   Einlesen von Header und Daten
flux, header = fits.getdata(file, header=True)

print('Minimum und Maximum im flux: ', flux.min(), '  ', flux.max())

#   Header-Check Spektrum
header_flag = input('Möchten Sie den header komplett sehen? j/n:')
if header_flag == 'j':
    print('Headerdaten: ')
    print(header)

#   Prüfung auf die nötigen header-Einträge
print('\nAusgabe der zur Wellenlängenberechnung nötigen Headereinträge:')
if 'NAXIS' in header:
    print('Dimension, NAXIS:                        ', header['NAXIS'])
else:
    print('Das ist kein 1d-Spektrum !')

if 'NAXIS1' in header:
    nax = header['NAXIS1']
    print('Anzahl der Werte (Abszisse), NAXIS1:     ', nax)
else:
    print('NAXIS1 fehlt im header !')

if 'CRVAL1' in header:
    crval = header['CRVAL1']
    print('Anfangs-Wellenlänge, CRVAL1:             ', crval)
else:
    print('CRVAL1 fehlt im header !')

if 'CDELT1' in header:
    cdel = header['CDELT1']
    print('Schrittweite der Wellenlänge, CDELT1:    ', cdel)
else:
    print('CDELT1 fehlt im header !')

if 'CRPIX1' not in header:
    header['CRPIX1'] = 1

#   Erzeugen von Arrays mit den Wellenlängen und Fluxes des Spektrums
#

wave = np.ones(nax, dtype=float)
crval = crval + (1 - header['CRPIX1']) * cdel
for i in range(nax):
    wave[i] = crval + i * cdel
#
# In der Liste wave sind die Wellenlängen der Pixel enthalten
# In der Liste flux die entsprechenden Intensitäten


# Plot gesamtes Spektrum
plt.style.use('seaborn-white')

fig = plt.figure(1, figsize=(10, 10))
plt.plot(wave, flux)
plt.xlabel('Wellenlänge [Angström]')
plt.ylabel('ADU')
plt.title('Spektrum '+file)
plt.grid(True)
fig.savefig(file.strip('.fit')+'.pdf')
# plt.show(block=True)
plt.pause(.1)

#   n Subplots
n = int(input('Gib die Anzahl der gewünschten Unterplots an: '))
if n > 1:
    pin = nax/n     # Anzahl der Werte pro Intervall

    fig, ax = plt.subplots(n, figsize=(10, 2*n))
    fig.subplots_adjust(bottom=0.15, left=0.2)
    for i in range(n):
        links = int(i*pin)
        rechts = int((i+1)*pin)
        ax[i].plot(wave[links:rechts], flux[links:rechts])
        ax[i].grid(True, linestyle='-.')
        ax[i].set_ylim(0.95*flux.min(), 1.05*flux.max())
    ax[0].set_title('Spektrum '+file)
    ax[n-1].set_xlabel('Wellenlänge [Angström]')
    fig.savefig(file.strip('.fit')+'_zoom.pdf')
    # plt.show(block=True)
    plt.pause(.1)

# Plot eines speziellen Wellenlängenbereichs
frage = input('Möchten Sie einen speziellen Wellenlängenbereich herausvergrößern?\n\
Anhand des bereits gespeicherten pdf kann der Bereich ausgesucht werden.j/n: ')
if frage == 'j':
    a = int(input('Anfang des Spektrums: '))
    b = int(input('Ende des Spektrums: '))
    aindex = int((a - crval) / cdel)
    bindex = int((b - crval) / cdel)
    fig = plt.figure(3, figsize=(10, 10))
    plt.plot(wave[aindex:bindex], flux[aindex:bindex])
    plt.xlabel('Wellenlänge [Angström]')
    plt.ylabel('ADU')
    plt.title('Spektrum ')
    fig.savefig(file.strip('.fit')+'Ausschnitt_1.pdf')

# plt.show(block=True)
plt.pause(.1)
