#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Erzeugt aus einer Folge von (baryzentrisch korrigierten) fits-1d-Spektren einen
farbigen dynamischen Spektrenplot, wobei die Ordinate die Beobachtungszeitpunkte bilden
(ein Headereintrag namens 'JD' muss im Header jeden Spektrums exstieren !).
Der abgebildete Wellenlängenbereich wird gewählt.

20231212
@author: lothar
"""

import numpy as np
from astropy.io import fits
import glob
import matplotlib.pylab as plt


# Fileliste erstellen. Spektren in einem (Unter)Ordner.
files = input("Pfad und Name der fits-Spektren (nutze wildcards) : ")
filelist = glob.glob(files)

# Alphabetisches Sortieren. Bei richtiger Namensgebung ergibt das eine
# zeitliche Ordnung
filelist.sort()

# Ausdruck der Liste
print("\nSpektrenliste: \n")
print("Anzahl der Spektren: ", len(filelist), "\n")

jd = np.zeros(len(filelist))
fluxmin = np.zeros(len(filelist))
fluxmax = np.zeros(len(filelist))
for i in range(len(filelist)):
    flux, header = fits.getdata(filelist[i], header=True)
    try:
        jd[i] = header['JD']
    except:
        print(filelist[i], ': Kein JD im Header')
        pass
print('JD-Bereich der Spektren reicht von ', jd.min(), ' bis ', jd.max())

# Eingabe des Bereichs der Beobachtungszeitpunkte als JD
JD_anfang = float(input('Geben Sie den Anfang des Bereichs der \
Beobachtungszeitpunkte (JD), den Sie abbilden möchten, ein: '))
JD_ende = float(input('Geben Sie das Ende des Bereichs der \
Beobachtungszeitpunkte (JD) ein: '))

# Eingabe des abgebildeten Wellenlängenbereichs
lambda_anfang = float(input('Geben Sie den Anfang des abzubildenden \
Wellenlängenbereichs an: '))
lambda_ende = float(input('Geben Sie das Ende des abzubildenden \
Wellenlängenbereichs an: '))

obj = input("Please enter the object name: ")

plt.figure(1, figsize=(7, 10))
plt.title("Linienentwicklung " + obj)
plt.grid(True)
plt.xlabel('Wellenlänge [Angström]')
plt.ylabel('JD')

# Bestimmung der minimalen und maximalen Fluxe
flux_bereich_max = np.zeros(len(filelist))
flux_bereich_min = np.zeros(len(filelist))

for i in range(len(filelist)):
    print(i, filelist[i], ' wird untersucht')
    flux, header = fits.getdata(filelist[i], header=True)
    if "CRPIX1" in header:
        refpix = int(header["CRPIX1"])
    else:
        refpix = 1
    step = float(header["CDELT1"])
    lambda0 = float(header["CRVAL1"]) - (refpix - 1) * step
    JD = header['JD']
    JD_float = float(JD)
    wave_bereich = np.array([])
    flux_bereich = np.array([])

    if (JD_float >= JD_anfang) and (JD_float <= JD_ende):
        wave = np.zeros(header['NAXIS1'])
        for k in range(header['NAXIS1']):
            wave[k] = lambda0 + k * step
            if wave[k] >= lambda_anfang and wave[k] <= lambda_ende:
                wave_bereich = np.hstack([wave_bereich, wave[k]])
                flux_bereich = np.hstack([flux_bereich, flux[k]])
        flux_bereich_max[i] = flux_bereich.max()
        flux_bereich_min[i] = flux_bereich.min()

for i in range(len(flux_bereich_max)):
    if flux_bereich_max[i] == False:
        flux_bereich_max[i] = 1.

for i in range(len(flux_bereich_min)):
    if flux_bereich_min[i] == False:
        flux_bereich_min[i] = 1.

# Abarbeiten der Spektrenliste:
for i in range(len(filelist)):
    print(i, filelist[i], ' wird bearbeitet')
    flux, header = fits.getdata(filelist[i], header=True)
    if "CRPIX1" in header:
        refpix = int(header["CRPIX1"])
    else:
        refpix = 1
    step = float(header["CDELT1"])
    lambda0 = float(header["CRVAL1"]) - (refpix - 1) * step
    JD = header['JD']
    JD_float = float(JD)
    wave_bereich = np.array([])
    flux_bereich = np.array([])

    if (JD_float >= JD_anfang) and (JD_float <= JD_ende):
        wave = np.zeros(header['NAXIS1'])
        for k in range(header['NAXIS1']):
            wave[k] = lambda0 + k * step
            if wave[k] >= lambda_anfang and wave[k] <= lambda_ende:
                wave_bereich = np.hstack([wave_bereich, wave[k]])
                flux_bereich = np.hstack([flux_bereich, flux[k]])

    plt.scatter(wave_bereich, np.full(len(flux_bereich), JD),
                c=flux_bereich, s=1., cmap='seismic')  # evtl. die Breite s der
    # Spektrenstreifen anpassen
    plt.clim(flux_bereich_min.min()*0.95, flux_bereich_max.max()*1.05)
    plt.pause(.1)

plt.figure(1)
plt.clim(flux_bereich_min.min()*0.95, flux_bereich_max.max()*1.05)
plt.colorbar()

plt.savefig("Linienentwicklung"+'JD_'+str(JD_anfang)+'_'+str(JD_ende)
            + '_'+str(lambda_anfang)+'_'+str(lambda_ende) + '.pdf', format='pdf')
plt.savefig("Linienentwicklung"+'JD_'+str(JD_anfang)+'_'+str(JD_ende)
            + '_'+str(lambda_anfang)+'_'+str(lambda_ende) + '.png', format='png')

print('Ende des Programms')
