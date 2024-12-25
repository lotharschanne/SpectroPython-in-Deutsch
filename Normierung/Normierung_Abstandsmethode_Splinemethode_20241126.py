#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Die 1d-Spektren im fits-Format werden eingelesen.
Im ersten Schritt werden für die Fluxe in kleinen Intervallen eine Regression
erster Ordnung (Gerade) berechnet und die Steigungen in einem Array gespeichert. 
Die Differenzen der Fluxe von Ende und Anfang der Intervalle werden berechnet und
durch die Steigung des Intervalls geteilt (normiert) und in einem Array gespeichert.
Dann werden die Elemente dieses Arrays gesucht, die
ein vorgegebenes Perzentil unterschreiten und als Stützpunkte für die Normierungs-
funktion gewertet werden, die anschließend mittels eines Splines berechnet wird.
Dabei können bestimmte Wellenlängenintervalle (breite Balmerlinien) von der 
Bildung von Stützpunkten ausgeschlossen werden.
Anschließend werden in einer Schleife mit mehreren Durchgängen Stützpunkte 
gelöscht, die zu weit oberhalb oder unterhalb des Splines liegen.
Wenn keine Stützpunkte in einem Spektrum gefunden werden, bricht das Programm
ab. Beispielsweise, wenn im Spektrum der Flux für viele Pixel Null oder sonst
einen konstanten Wert beträgt.

20241112
@author: Lothar Schanne
"""

import numpy as np
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import glob
from scipy.interpolate import make_smoothing_spline
from scipy.optimize import curve_fit


# Definition von Funktionen:

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


def gerade(x, a, b):
    y = a + b * x
    return y


plt.switch_backend('Qt5Agg')  # Interaktives GrafikBackend einschalten

# Create file list for object spectra. All spectra in one (sub)folder.
files = input('Path and name of the OBJECT files (use wildcards) : ')
filelist = glob.glob(files)
filelist.sort()

print('Number of spectra: ', len(filelist), '\n')

forderung = input(
    "\nDo you want to remove interpolation points in certain wavelength ranges?\
 Then enter y: "
)
if forderung == 'y':
    forderung1 = input(
        'Wellenlängenbereiche manuell oder per Liste definieren? m oder l eingeben: ')
else:
    forderung1 = ''

# kann anders gewählt werden = Halber (Abstand - 1) der Pixel, deren
abstand = 2
# Differenzen berechnet werden.

# Wahl des Perzentils der Differenzen, die als primäre Stützpunkte gelten
perc = 10  # in %, kann angepasst werden, so klein wie möglich wählen


for i in range(len(filelist)):
    flux, header = fits.getdata(filelist[i], header=True)
    wave, flux = wavecalc(flux, header)
    print()
    print('Spektrum: ', filelist[i])
    print(f'Wellenlänge des ersten Pixels: {wave[0]:.2f}')
    print(f'Wellenlänge des letzten Pixels:  {wave[-1]:.2f}')
    print('Anzahl der Pixel: ', header['NAXIS1'])
    print('Schrittweite', header['CDELT1'])

    wave_points = np.zeros(header['NAXIS1']-2*abstand)
    delta = np.zeros(header['NAXIS1']-2*abstand)
    steigungen = []

    # Berechnung der Steigungen in den Intervallen
    for j in range(abstand, header['NAXIS1']-2*abstand):
        fl = flux[j-abstand:j+abstand+1]
        wa = wave[j-abstand:j+abstand+1]
        popt, pcov = curve_fit(gerade, wa, fl)
        delta[j] = abs(flux[j+abstand] - flux[j-abstand]) / abs(popt[1])
        steigungen = np.append(steigungen, popt[1])
        wave_points[j] = wave[j]

    plt.figure()
    plt.plot(wave[abstand:-2*abstand], steigungen, linewidth=.5)
    plt.grid(True)
    plt.title(filelist[i].rsplit('.')[0] + '  Steigungen der Intervalle')
    # plt.ylim(-3,3)
    plt.pause(.1)
    plt.savefig(filelist[i].rsplit('.')[0]+'_Steigungen.png')
    plt.close()

    # plt.figure()
    # plt.hist(delta, bins=100)
    # plt.title(filelist[i].rsplit('.')[0]
    #           + '   Histogramm der normierten Differenzen benachbarter Pixel')
    # plt.grid(True)
    # plt.pause(.1)
    # plt.savefig(filelist[i].rsplit('.')[0]+'_Histogramm.png')
    # plt.close()

    indices = delta < np.percentile(delta, perc)
    flux_points = flux[abstand:-abstand][indices]
    wave_points = wave[abstand:-abstand][indices]

    if len(wave_points) == 0:
        print(filelist[i] + ': keine Stützpunkte gefunden')
        break

    plt.figure(figsize=(12, 8))
    plt.plot(wave, flux)
    plt.scatter(wave_points, flux_points, s=10, color='r')
    plt.title(filelist[i].rsplit('.')[0] + 'Lage der primären Stützpunkte')
    plt.grid(True)
    plt.pause(.1)
    plt.savefig(filelist[i].rsplit('.')[0] + '_Primärstützpunkte.png')
    plt.close()
    print(str(perc)+'%-Perzentil der Abstände: ', np.percentile(delta, perc))
    print("Anzahl der Stützpunkte-Iteration am Anfang: ", len(wave_points))

   # ++++++++++++++++++ Exclude wavelength ranges

    loesch1 = np.zeros(100, dtype=float)
    loesch2 = np.zeros(100, dtype=float)

    if forderung1 == 'l':
        # Diese Liste der breiten Absorptionslinien sollte an die Breite der
        # Linien des jeweiligen Spektrums angepasst werden

        # Breite Linien im Spektrum:
        loeschintervalle = [(3960, 3980), (4080, 4125),
                            (4320, 4360), (4700, 4720),
                            (4820, 4910),
                            (6500, 6620)]

        # Schmale Linien im Spektrum:
        # loeschintervalle = [(3960, 3980), (4090, 4115),
        #                     (4330, 4350), (4700, 4720),
        #                     (4850, 4870),
        #                     (6540, 6585)]

        for j in range(len(loeschintervalle)):
            loesch1[j], loesch2[j] = loeschintervalle[j]
            for n in np.arange(len(wave_points)):
                if loesch1[j] < wave_points[n] and loesch2[j] > wave_points[n]:
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

    # Creation of a cubic spline from the points (= normfunction)
    # Rigizität des Splines:
    smoothing = 200  # anpassen !!!!!!!!!!!!!!!!!!!!!!!!!!!

    normfunc = make_smoothing_spline(
        wave_points[:], flux_points[:], lam=smoothing)

    # Iterierende Löschung von abweichenden Stützpunkten:
    # Festlegung des Abstands von der Normierungsfunktion:
    s1 = .01  # unter der Normierunngsfunktion liegende Stützpunkte
    s2 = .05  # über  der Normierunngsfunktion liegende Stützpunkte

    for k in range(3):  # Wenn nicht ausreichend Stützpunkte übrig bleiben
        # range verkleinern
        for n in range(len(wave_points)):
            for m in range(len(wave)):
                if wave[m] == wave_points[n]:
                    if (flux_points[n] < (1 - s1) * normfunc(wave[m])
                            or flux_points[n] > (1 + s2) * normfunc(wave[m])):
                        flux_points[n] = 0
                        wave_points[n] = 0
                    break

        if len(flux_points[flux_points != 0]) <= 5:
            print('Iteration', k, ':   Anzahl der Stützpunkte: ',
                  len(flux_points[flux_points != 0]))
            print('Abbruch der Iteration')
            break

        flux_points = flux_points[flux_points != 0]
        wave_points = wave_points[wave_points != 0]

        # flux_points[0] = normfunc(wave[0])
        flux_points[0] = np.mean(flux[0:10])
        wave_points[0] = wave[0]
        # flux_points[-1] = normfunc(wave[-1])
        flux_points[-1] = np.mean(flux[-1:-10:-1])
        wave_points[-1] = wave[-1]

        print('Iteration', k, ':   Anzahl der Stützpunkte: ', len(wave_points))

        normfunc = make_smoothing_spline(
            wave_points[:], flux_points[:], lam=smoothing)

    normfunction = normfunc(wave)

    plt.figure(figsize=(10, 8))
    plt.plot(wave, flux, 'g-', linewidth=.5)
    plt.plot(wave, normfunction, 'b-')
    plt.scatter(wave_points, flux_points, color='r', s=8)
    plt.title(filelist[i].rsplit('.')[0] +
              ' letzliche Stützpunkte und Normierungsspline')
    plt.grid(True)
    plt.pause(.1)
    plt.savefig(filelist[i].rsplit('.')[0] +
                '_letzliche Stützpunkte und Normierungsspline.png')
    plt.close()

    # Abschließende Normierung auf das Kontinuum (= normfunction)
    normierter_flux = flux / normfunction

    plt.figure(figsize=(10, 8))
    plt.plot(wave, normierter_flux, 'k', linewidth=.5)
    plt.title(filelist[i].rsplit('.')[0] + ' normiert')
    plt.pause(.1)
    plt.grid(True)
    plt.savefig(filelist[i].rsplit('.')[0] + '_normiert.png')
    plt.close()

    header["CRVAL1"] = wave[0]
    header["NAXIS1"] = len(wave)
    header["CRPIX1"] = 1
    filename = filelist[i].rsplit('.')[0]
    name = filename + "_normalized.fits"

    fits.writeto(
        name, normierter_flux, header, overwrite=True, output_verify="silentfix"
    )

    ascii.write([wave, normierter_flux], filename +
                "_normalized"+'.dat', overwrite=True,
                names=['WAVE', 'FLUX'], format='tab')


# Programmende
plt.close('all')
print("END of program.")
