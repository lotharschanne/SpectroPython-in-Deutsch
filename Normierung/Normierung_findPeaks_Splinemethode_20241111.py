#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Die 1d-Spektren einer Serie im fits-Format werden eingelesen.
Dazu werden mit einer peakfinder-Funktion lokale peaks gefunden und als 
primäre Stützpunkte verwendet.
Aus diesen primären Stützpunkten wird ein Spline berechnet, der als Normierungs-
funktion dient.
Danach können bestimmte Wellenlängenintervalle (Balmerlinien) von der Bildung von
Stützpunkten ausgeschlossen werden.
Dann werden in einer Schleife mit mehreren Durchgängen Stützpunkte gelöscht, die
zu weit oberhalb oder unterhalb des immer neu berechneten Splines liegen.
Der letzte berechnete Spline dient als Normierungsfunktion.
Die normierten Spektren werden als fits und als ascii-Datei abgespeichert.

20241111
@author: Lothar Schanne
"""

import numpy as np
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import glob
from scipy.interpolate import make_smoothing_spline
from scipy.signal import find_peaks


def wavecalc(flux, header):
    """Berechnet die Wellenlängen aus den Headerdaten der fits-Datei."""
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
print('Number of spectra: ', len(filelist), '\n')

forderung = input(
    'Möchten Sie Stützpunkte in Wellenlängenbereichen manuell (m) oder per Liste (l) \
löschen? Dann m oder l eingeben: ')

# ********************** Normierung ****************************


# Read headers and fluxes
for i in range(len(filelist)):
    print()
    print(filelist[i], ':')

    # für Object fits-Spektren:
    flux, header = fits.getdata(filelist[i], header=True)
    wave, flux = wavecalc(flux, header)

    flux_mask = flux == 0
    flux[flux_mask] = 1e-10

    print('\nObjektspektrum: ', filelist[i])

    # Auffinden von Peaks:
    # peaks, params = find_peaks(flux, threshold=(.00000001, 100), distance=3)
    peaks, params = find_peaks(flux, distance=5)  # Die Distanz bezeichnet den
    # Mindestabstand in Pixeln, die benachbarte Peaks haben müssen.

    wave_points = np.array([])
    flux_points = np.array([])
    wave_points = np.append(wave_points, wave[0])
    flux_points = np.append(flux_points, flux[0])
    for x in peaks:
        if x > 1 and x < len(wave) - 1:
            f = np.mean(flux[x-1:x+2])
            flux_points = np.append(flux_points, f)
            wave_points = np.append(wave_points, wave[x])

    wave_points = np.append(wave_points, wave[-1])
    flux_points = np.append(flux_points, flux[-1])

    flux_points[0] = np.median(flux[0:5])
    flux_points[-1] = np.median(flux[-1:-5:-1])

    print("Anzahl der Stützpunkte-Iteration am Anfang: ", len(wave_points))

    # ++++++++++++++++++ Exclude wavelength ranges:
    loesch1 = np.zeros(20, dtype=float)
    loesch2 = np.zeros(20, dtype=float)

    if forderung == 'l':
        # Diese Liste der breiten Absorptionslinien sollte an die Breite der
        # Linien des jeweiligen Spektrums angepasst werden.

        # Breite Linien:
        # loeschintervalle = [(3960, 3980), (4080, 4125),
        #                     (4320, 4360), (4700, 4720),
        #                     (4820, 4910),
        #                     (6500, 6620)]

        # Schmale Linien:
        loeschintervalle = [(3960, 3980), (4090, 4115),
                            (4330, 4350), (4700, 4720),
                            (4850, 4870),
                            (6540, 6585)]

        for j in range(len(loeschintervalle)):
            loesch1[j], loesch2[j] = loeschintervalle[j]
            for n in np.arange(len(wave_points)):
                if loesch1[j] <= wave_points[n] and loesch2[j] >= wave_points[n]:
                    flux_points[n] = 0

    if forderung == 'm':
        j = 0
        while forderung:
            loesch1[j] = float(input("Initial wavelength of the range? :"))
            loesch2[j] = float(input("End wavelength of the range? :"))
            for n in np.arange(len(wave_points)):
                if loesch1[j] <= wave_points[n] and loesch2[j] >= wave_points[n]:
                    flux_points[n] = 0
            j += 1
            forderung = input(
                "Do you want to remove further support points? Then enter y: ")

    wave_points = wave_points[flux_points != 0]
    flux_points = flux_points[flux_points != 0]

    flux_points[0] = np.median(flux[0:5])
    wave_points[0] = wave[0]
    flux_points[-1] = np.median(flux[-1:-5:-1])
    wave_points[-1] = wave[-1]

    print("Anzahl der Stützpunkte Iteration, nach Bereichslöschung: ",
          len(wave_points))

    plt.figure(figsize=(10, 8))
    plt.plot(wave, flux, 'g-', linewidth=.5, label='Spektrum')
    plt.scatter(wave_points, flux_points, color='r', s=6, label='Stützpunkte')
    plt.title(filelist[i] + ': Stützpunkte nach Bereichslöschung und Spektrum')
    plt.grid(True)
    plt.legend()
    plt.pause(.1)
    plt.savefig(filelist[i]+'_Anfangsstützpunkte.png')
    plt.close()

    # *****************************************************************
    # Iterierende Löschung von abweichenden Stützpunkten:
    # Festlegung des zur Löschung erforderlichen Abstands von der
    # Normierungsfunktion:
    s1 = .01  # unter der Normierungsfunktion liegende Stützpunkte
    s2 = .02  # über  der Normierungsfunktion liegende Stützpunkte, für sehr
    # linienreiche Spektren 0.2 oder mehr wählen, für Spektren mit
    # Emissionslinien 0.01 wählen, für normale Spektren ist s2=0.05 eine gute Wahl.

    K = 5  # Anzahl der Iterationen zum Löschen von abweichenden Stützpunkten,
    # evtl. anpassen
    smoothing = 10  # Anfangs-Rigizität des Splines, evtl. anpassen !!!!!!!!!!

    for k in range(K):
        smoothing = smoothing * (k+1)/1.4  # Der Nenner des Bruchs bestimmt den
        # Zuwachs des smoothing-Parameters während der Iteration, evtl. anpassen

        # if k == K:
        #     smoothing = 500
        #     s1 = .01
        #     s2 = .01

        normfunc = make_smoothing_spline(
            wave_points[:], flux_points[:], lam=smoothing)

        for n in range(1, len(wave_points)-1):
            for m in range(len(wave)):
                if wave[m] >= wave_points[n]:
                    if (flux_points[n] < (1 - s1) * normfunc(wave[m])
                            or flux_points[n] > (1 + s2) * normfunc(wave[m])):
                        flux_points[n] = 0
                    break

        wave_points = wave_points[flux_points != 0]
        flux_points = flux_points[flux_points != 0]

        flux_points[0] = np.median(flux[0:5])
        wave_points[0] = wave[0]
        flux_points[-1] = np.median(flux[-1:-5:-1])
        wave_points[-1] = wave[-1]

        print('Iteration', k, ':   Anzahl der Stützpunkte: ', len(wave_points),
              'smooth = ', smoothing)
        if len(wave_points) < 20:
            break

    smoothing = 200
    normfunc = make_smoothing_spline(
        wave_points[:], flux_points[:], lam=smoothing)
    normfunction = normfunc(wave)

    flux_normiert = flux / normfunction

    plt.figure(figsize=(10, 8))
    plt.plot(wave, flux, 'g-', linewidth=.5, label='Spektrum')
    plt.plot(wave, normfunction, 'b-', label='Normierungsfunktion')
    plt.scatter(wave_points, flux_points, color='r', s=6, label='Stützpunkte')
    plt.title(filelist[i] + ': Stützpunkte, Normierungsfunktion und Spektrum')
    plt.grid(True)
    plt.legend()
    plt.pause(.1)
    plt.savefig(filelist[i]+'nachIteration.png')
    plt.close()

    plt.figure(figsize=(10, 8))
    plt.plot(wave, flux_normiert,  'b-', label=filelist[i], linewidth=.6)
    plt.grid(True)
    plt.title(filelist[i] + ' normiert', fontsize=12)  # fontsize anpassen !!!
    plt.xlabel('Wavelength in Angström')
    plt.ylabel('normierter Flux')
    plt.ylim(0, 1.5)  # anpassen, falls Emissionslinien vorhanden sind
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.savefig(filelist[i].rsplit('.')[0] + '_normiert.png', format='png')
    plt.pause(.1)
    plt.close()

    fits.writeto(filelist[i].rsplit('.')[0] + '_normiert.fits',
                 flux_normiert, header,
                 overwrite=True, output_verify='silentfix')

    ascii.write([wave, flux_normiert],
                filelist[i].rsplit('.')[0] + "normiert.dat",
                overwrite=True, names=['WAVE', 'FLUX'], format='tab')

# Programmende
# plt.close('all')
print("END of program.")
