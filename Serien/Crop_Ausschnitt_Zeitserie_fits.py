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
from astropy.io import fits
import matplotlib.pyplot as plt
import glob

plt.switch_backend('Qt5Agg')  # Interaktives GrafikBackend einschalten


# Fileliste erstellen. Spektren in einem (Unter)Ordner.
files = input("Pfad und Name der Spektren (nutze wildcards) : ")
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
for i in range(len(filelist)):
    flux, header = fits.getdata(filelist[i], header=True)
    # Grafik mit allen Spektrenauschnitten:
    wave = np.zeros(header["NAXIS1"], dtype=float)
    for k in range(len(wave)):
        wave[k] = header["CRVAL1"] + \
            (k + 1 - header["CRPIX1"]) * header["CDELT1"]
    plt.plot(wave, flux, linewidth=1)
    print(filelist[i], wave[0])
print('Bitte warten. Berechnungen laufen.')

# Ploteigenschaften (Gitter, Beschriftungen etc.) anpassen:
# plt.title('Zeitserie von '+header['OBJECT'])
plt.grid(True)
plt.xlabel("Wellenlänge [Angström]")
# # Legende kann angepasst werden, z.B. Schriftgröße (fontsize)
# plt.legend(
#     bbox_to_anchor=(0.0, 1.02, 1.0, 0.102),
#     loc=3,
#     ncol=4,
#     mode="expand",
#     borderaxespad=0.0,
#     fontsize=6,
# )

# plt.pause(20)  # Während dieser Zeit kann die Grafik aktiv bearbeitet werden

frage = input(
    'Möchten Sie die Grafik als png oder pdf speichern? Dann bitte "png" oder "pdf" eingeben: ')
if frage == 'pdf':
    fig.savefig(files+'_Zeitserie.pdf')
elif frage == 'png':
    fig.savefig(files+'_Zeitserie.png')


# Generate the spectrum section, plot and save as fit
min = 0
max = 10000
for i in range(len(filelist)):
    flux, header = fits.getdata(filelist[i], header=True)
    if 'CRPIX1' not in header:
        header['CRPIX1'] = 1
    if header["CRVAL1"] + \
            (1 - header["CRPIX1"]) * header["CDELT1"] >= min:
        min = header["CRVAL1"] + \
            (1 - header["CRPIX1"]) * header["CDELT1"]
    if header["CRVAL1"] + \
        (1 - header["CRPIX1"]) * header["CDELT1"] + \
            header["CDELT1"] * header["NAXIS1"] <= max:
        max = header["CRVAL1"] + \
            (1 - header["CRPIX1"]) * header["CDELT1"] +\
            header["CDELT1"] * header["NAXIS1"]
print("\nGemeinsamer Wellenlängenbereich: ",
      int(min), ' bis ', int(max),  ' Angström')

print("\nSpecification of the wavelength range to be transferred ")
a = float(input("Begin: "))
b = float(input("End: "))
n = 0

for i in range(len(filelist)):
    flux, header = fits.getdata(filelist[i], header=True)
    if 'CRPIX1' not in header:
        header['CRPIX1'] = 1
    header["CRVAL1"] = header["CRVAL1"] + \
        (1 - header["CRPIX1"]) * header["CDELT1"]
    if header["CRVAL1"] <= a and header["CRVAL1"] +\
            header["CDELT1"] * header["NAXIS1"] >= b:
        aindex = int((a - header["CRVAL1"]) / header["CDELT1"])
        bindex = int((b - header["CRVAL1"]) / header["CDELT1"])
        filename = (
            filelist[i].rsplit(".", 1)[0] + "_" + str(int(a)) +
            "_" + str(int(b)) + ".fits"
        )
        newflux = flux[aindex:bindex]
        header["CRVAL1"] = header["CRVAL1"] + aindex * header["CDELT1"]
        header["CRPIX1"] = 1
        header["NAXIS1"] = bindex - aindex
        print('Verwendet: ', n, i,
              header["CRVAL1"], header["NAXIS1"], ' Pixel')
        fits.writeto(filename, newflux, header,
                     overwrite=True, output_verify="silentfix")
        n += 1

print('Ende des Programms')
plt.close('all')
