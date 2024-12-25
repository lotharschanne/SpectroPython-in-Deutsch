#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Die 1d-Spektren einer Serie im fits-Format werden eingelesen, in gleiche 
Intervalle geteilt und ein Perzentil des Flux und die mittlere Wellenlänge 
darin berechnet und als Stützpunkt für die Normierung gewertet. 
Dabei können bestimmte Wellenlängenintervalle (breite Balmerlinien) von der 
Bildung von Stützpunkten ausgeschlossen werden. Aus den Stützpunkten wird ein 
Spline berechnet.

Die Einteilung in gleiche Intervalle stellt sicher, dass Stützpunkte im gesamten
Spektrum gebildet werden, allerdings mit dem Nachteil, dass sie auch innerhalb
von Linien liegen können.

Dann werden in einer Schleife mit mehreren Durchgängen Stützpunkte gelöscht, die
zu weit oberhalb oder unterhalb des Splines liegen. Da mehr Punkte unterhalb
des Splines gelöscht werden wie oberhalb, wandert der Spline bei jedem
Schleifendurchgang in Richtung des Quasi-Kontinuums.
Am Schluß wird die Intensitäten (flux) durch den Spline geteilt und damit die 
Normierung auf das Quasi-Kontinuum erreicht.
Die Ergebnisse werden grafisch 0,1 Sekunde lang dargestellt, die Grafiken
gespeichert und die normierten Ordnungen als fits und ascii-Datei abgespeichert.

20241111
@author: Lothar Schanne
"""

import numpy as np
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import glob
from scipy.interpolate import make_smoothing_spline


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
files = input('Path and name of the OBJECT files (use wildcards) : ')
filelist = glob.glob(files)
filelist.sort()

# Print the list
# print('\nList of spectra: \n')
print('Number of spectra: ', len(filelist), '\n')

forderung = input(
    "Möchten Sie Stützpunkte in den breiten Linien löschen? Dann y eingeben: "
)
if forderung == 'y':
    forderung1 = input(
        'Wellenlängenbereiche manuell (m) oder per Liste (l) definieren? m oder l eingeben: ')


# Abarbeiten der filelist
for i in range(len(filelist)):
    print()
    print(filelist[i], ':')

    flux, header = fits.getdata(filelist[i], header=True)
    wave, flux = wavecalc(flux, header)

    inter = int(header['NAXIS1'] / 100)  # Die Zahl der berechneten Intervalle

    print(f'Objekt Wellenlänge des ersten Pixels: {wave[0]:.2f}')
    print('Anzahl der Pixel: ', header['NAXIS1'])
    print('Schrittweite', header['CDELT1'])
    print('Intervallbreite = ', inter)

    wave_points = np.array([])
    flux_points = np.array([])

    wave_points = np.append(wave_points, wave[0])
    flux_points = np.append(flux_points, flux[0:5].mean())

    for j in range(int(len(flux) / inter) + 1):
        index = j*inter
        if index <= len(flux) - inter:
            wave_points = np.append(
                wave_points, wave[index:index+inter].mean())
            # flux_points = np.append(flux_points, flux[index:index+inter].max())
            # flux_points = np.append(flux_points, flux[index:index+inter].mean())
            flux_points = np.append(
                flux_points, np.percentile(flux[index:index+inter], 80))
            # Perzentil evtl. anpassen !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        else:
            wave_points = np.append(wave_points, wave[-1])
            flux_points = np.append(flux_points, flux[-1:-5:-1].mean())

    plt.figure(figsize=(12, 8))
    plt.plot(wave, flux, linewidth=.5, label='Spektrum')
    plt.scatter(wave_points, flux_points, s=12, color='r', label='Stützpunkte')
    plt.title('Lage der primären Stützpunkte ' + filelist[i])
    plt.grid(True)
    plt.legend()
    plt.pause(.1)
    plt.savefig(filelist[i].rsplit('.')[0] + '_Primärstützpunkte.png')
    plt.close()

    print("Anzahl der Stützpunkte-Iteration am Anfang: ", len(wave_points))

    # ++++++++++++++++++ Exclude wavelength ranges:
    loesch1 = np.zeros(20, dtype=float)
    loesch2 = np.zeros(20, dtype=float)

    if forderung1 == 'l':
        # Diese Liste der breiten Absorptionslinien sollte an die Breite der
        # Linien des jeweiligen Spektrums angepasst werden

        # Breite Linien im Spektrum:
        # loeschintervalle = [(3960, 3980), (4080, 4125),
        #                     (4320, 4360), (4700, 4720),
        #                     (4820, 4910),
        #                     (6500, 6620)]

        # Schmale Linien im Spektrum:
        loeschintervalle = [(3960, 3980), (4090, 4115),
                            (4330, 4350), (4700, 4720),
                            (4850, 4870),
                            (6540, 6585)]

        for j in range(len(loeschintervalle)):
            loesch1[j], loesch2[j] = loeschintervalle[j]
            for n in np.arange(len(wave_points)):
                if loesch1[j] <= wave_points[n] and loesch2[j] >= wave_points[n]:
                    flux_points[n] = 0

    if forderung1 == 'm':
        j = 0
        while forderung:
            loesch1[j] = float(input("Initial wavelength of the range? :"))
            loesch2[j] = float(input("End wavelength of the range? :"))
            for n in np.arange(len(wave_points)):
                if loesch1[j] < wave_points[n] and loesch2[j] > wave_points[n]:
                    flux_points[n] = 0
            j += 1
            forderung = input(
                "Do you want to remove further support points? Then enter y: ")

    loesch1 = loesch1[loesch1 != 0]
    loesch2 = loesch2[loesch2 != 0]

    wave_points = wave_points[flux_points != 0]
    flux_points = flux_points[flux_points != 0]

    print("Anzahl der Stützpunkte Iteration, nach Bereichslöschung: ", len(wave_points))

    # Iterierende Löschung von abweichenden Stützpunkten:
    # Festlegung des Abstands von der Normierungsfunktion:
    s1 = .01  # unter der Normierungsfunktion liegende Stützpunkte
    s2 = .02  # über  der Normierungsfunktion liegende Stützpunkte
    # bei Emissionslinienspektren s2= 0.01 wählen, ansonsten 0.1
    K = 3  # Anzahl der Iterationen zum Löschen von abweichenden Stützpunkten
    # evtl. anpassen !!!!!
    smoothing = 100  # Anfangs-Rigizität des Splines, evtl. anpassen !!!!!!!!!!

    normfunc = make_smoothing_spline(
        wave_points[:], flux_points[:], lam=smoothing)

    for k in range(K):
        smoothing = smoothing * (k+1)/1.4  # Der Nenner des Bruchs bestimmt den
        # Zuwachs des smoothing-Parameters während der Iteration, evtl. anpassen
        # if k == K-1:
        #     smoothing = 500
        #     s1 = .01
        #     s2 = .01

        for n in range(1, len(wave_points)-1):
            for m in range(len(wave)):
                if wave[m] >= wave_points[n]:
                    if (flux_points[n] < (1 - s1) * normfunc(wave[m])
                            or flux_points[n] > (1 + s2) * normfunc(wave[m])):
                        flux_points[n] = 0
                        wave_points[n] = 0
                    break

        wave_points = wave_points[flux_points != 0]
        flux_points = flux_points[flux_points != 0]

        flux_points[0] = flux[0:4].mean()
        wave_points[0] = wave[0]
        # flux_points[0] = normfunc(wave[0])
        flux_points[-1] = flux[-1:-4:-1].mean()
        # flux_points[-1] = (normfunc(wave[-1]) + flux[-1:-4:-1].mean()) / 2
        wave_points[-1] = wave[-1]

        print('Iteration', k, ':   Anzahl der Stützpunkte: ', len(wave_points),
              'Smoothing: ', smoothing)

        normfunc = make_smoothing_spline(
            wave_points[:], flux_points[:], lam=smoothing)

    normfunction = normfunc(wave)

    plt.figure(figsize=(10, 8))
    plt.plot(wave, flux, 'g-', linewidth=.5, label='Spektrum')
    plt.plot(wave, normfunction, 'b-', label='Normierungsfunktion')
    plt.scatter(wave_points, flux_points, color='r', s=12, label='Stützpunkte')
    plt.title(filelist[i] + 'Normierungsfunktion und Stützpunkte')
    plt.grid(True)
    plt.legend()
    plt.pause(.1)
    plt.savefig(filelist[i].rsplit('.')[0]+'_Normierungsfunktion.png')
    # plt.show(block=True)
    plt.close()

    flux_normiert = flux / normfunction

    plt.figure(figsize=(10, 8))
    plt.plot(wave, flux_normiert,  'b-',
             label=filelist[i] + ' normiert', linewidth=.6)
    plt.grid(True)
    plt.title(filelist[i] + ' normiert', fontsize=12)  # fontsize anpassen !!!
    plt.xlabel('Wavelength in Angström')
    plt.ylabel('normierter Flux')
    plt.ylim(0, 1.2)  # anpassen, falls Emissionslinien vorhanden sind
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.legend()
    plt.pause(.1)
    plt.savefig(filelist[i].rsplit('.')[0] + '_normiert.png', format='png')
    plt.close()

    fits.writeto(filelist[i].rsplit('.')[0] + '_normiert.fits',
                 flux_normiert, header,
                 overwrite=True, output_verify='silentfix')

    ascii.write([wave, flux_normiert],
                filelist[i].rsplit('.')[0] + "normiert.dat",
                overwrite=True, names=['WAVE', 'FLUX'], format='tab')

# Programmende
plt.close('all')
print("END of program, all done")
