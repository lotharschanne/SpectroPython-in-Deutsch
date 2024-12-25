#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Crop_Ausschnitt_Zeitserie_fits.py

Fileliste erstellen für  wellenlängenkalibrierte 1d_Spektren einer Zeitreihe im
fits-Format. Plotten der Spektren und Bestimmung des gemeinsamen
Wellenlängenbereichs.
Beschneiden aller Spektren auf einen wählbaren Wellenlängenbereich und
abspeichern aller beschnittenen als fits mit Wellenlängenbereich im Dateinamen

20220125
@author: lothar
"""

import numpy as np
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import glob

plt.switch_backend('Qt5Agg')  # Interaktives GrafikBackend einschalten


# Fileliste erstellen. Spektren in einem (Unter)Ordner.
files = input("Pfad und Name der dat-Spektren (nutze wildcards) : ")
filelist = glob.glob(files)

# Alphabetisches Sortieren. Bei richtiger Namensgebung ergibt das eine
# zeitliche Ordnung
filelist.sort()

# Ausdruck der Liste
# print("\nSpektrenliste: \n")
# print(filelist)
print("Anzahl der Spektren: ", len(filelist), "\n")
print('Bitte warten. Berechnungen laufen.')


fig = plt.figure(figsize=(14, 20))
# Einlesen von flux und header:

tab = ascii.read(filelist[0], format='tab')
print('Tabellenspalten-Namen: ', tab.colnames)
flux = tab['FLUX']
wave = tab['WAVE']
plt.plot(wave, flux, linewidth=1)
print(filelist[0], wave[0], wave[-1])
print('Bitte warten. Berechnungen laufen.')

plt.grid(True)
plt.xlabel("Wellenlänge [Angström]")

plt.pause(20)  # Während dieser Zeit kann die Grafik aktiv bearbeitet werden.

# Generate the spectrum section, plot and save as fit
min = 0
max = 10000
for i in range(len(filelist)):
    wave = ascii.read(filelist[i])['WAVE']
    minimum = wave[0]
    maximum = wave[-1]
    if minimum > min:
        pass
    else:
        min = minimum
    if maximum < max:
        pass
    else:
        max = maximum

print("\nGemeinsamer Wellenlängenbereich: ",
      int(min), ' bis ', int(max),  ' Angström')

print("\nSpecification of the wavelength range to be transferred ")
a = float(input("Begin: "))
b = float(input("End: "))

plt.close()

for i in range(len(filelist)):
    tab = ascii.read(filelist[i], format='tab')
    print('Tabellenspalten-Namen: ', tab.colnames)
    flux = tab['FLUX']
    wave = tab['WAVE']
    for j in range(len(wave)):
        if wave[j] >= a:
            aindex = j
            break
    for j in range(len(wave)):
        if wave[j] >= b:
            bindex = j
            break

    newflux = flux[aindex:bindex]
    newwave = wave[aindex:bindex]

    ascii.write([newwave, newflux],
                filelist[i].rstrip('.dat')+'_'+str(a)+'_'+str(b)+'.dat',
                overwrite=True,
                names=["WAVE", "FLUX"],
                format="tab",)

    plt.figure()
    plt.plot(newwave, newflux, 'k-', linewidth=1)
    plt.title(filelist[i], fontsize=9)
    plt.xlabel('Wellenlänge [Angström]')
    plt.ylabel('normierter Flux')
    plt.ylim(-.2, .2)  # anpassen !!!!!!!!!!!!!!!!!!
    plt.grid(True)
    plt.savefig(filelist[i].rstrip('.dat')+'_'+str(a)+'_'+str(b)+'.png',
                format='png')
    plt.pause(.5)
    plt.close()


print('Ende des Programms')
plt.close('all')
