#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Berechnet für eine Zeitserie AUF DAS KONTINUUM NORMIERTER SPEKTREN im
fits-Format die Äquivalentweite einer Linie, die Linientiefe, die FWHM und die
Wellenlänge des Linienminimums.
Eingabe der Integrationsgrenzen manuell oder grafisch.
Normierungsfehler werden durch eine lineare Renormierungroutine im
Integrationsbereich kompensiert. Dabei wird der Flux an den
Integrationsgrenzen für die Kontinuumsberechnung genommen und dazwischen linear
interpoliert und mit der so gewonnenen linearen Kontinuumsfunktion erneut
normiert. Daher ist es wichtig, dass die Integrationsgrenzen auch wirklich
auf dem Kontinuum liegen.

Die Berechnungen setzen voraus, dass für alle Spektren der Serie das gleiche
Wellenlängenintervall für die Berechnung des Integrals verwendet werden kann,
also keine wesentlichen RV-Änderungen stattfinden, oder wenn doch, dass die
Linie isoliert ist (also der Flux = 1 im Umfeld ist). Deshalb am besten
baryzentrisch korrigierte Spektren benutzen.

Das erste Spektrum wird geplottet und für 20 Sekunden angezeigt. In diesem
Zeitraum ist das Grafikfenster interaktiv geschaltet, so dass man das Spektrum
vergrößern kann und den Integrations-Wellenlängenbereich für die EW-Berechnung
optisch aussuchen kann. Dieses Reaktionszeitfenster von 20 Sekunden kann in
Zeile 91 geändert werden. Falls die grafische/interaktive Eingabe der Integrations-
grenzen gewählt wurde, ist nach der Pause mit zwei Mausklicks im Grafikfenster
der Integrationsbereich zu wählen, ansonsten sind die Integrationsgrenzen als
Zahlen einzugeben.

Die gewählte Linie wird für alle Spektren als schwarze Linie geplottet,
dazu die renormierte in blauer Farbe und die interpolierte in rot.
Auf Nachfrage können die plots auch als pdf abgespeichert werden.

Version 2023-05-08

@author: lothar schanne
"""

import numpy as np
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import glob
from scipy.interpolate import Rbf

plt.switch_backend('Qt5Agg')  # Interaktives GrafikBackend einschalten


# Fileliste erstellen. Spektren in einem (Unter)Ordner.
files = input("Pfad und Name der Spektren (nutze wildcards) : ")
filelist = glob.glob(files)

# Frage ob die Grafiken gespeichert werden sollen
frage = input('Sollen die Grafiken der Spektrenausschnitte als PDF gespeichert\
 werden?. Wenn ja bitte "j" eingeben, ansonsten einen anderen \
Buchstaben oder return: ')

# Alphabetisches Sortieren. Bei richtiger Namensgebung ergibt das eine
# zeitliche Ordnung
filelist.sort()

# Ausdruck der Liste
print("\nSpektrenliste: \n")
print(filelist)
print("Anzahl der Spektren: ", len(filelist), "\n")


# Grafik erstes Spektrum
sp = fits.open(filelist[0], ignore_missing_end=True)
# print('\n\nHeader of the spectrum :\n\n', sp[0].header, '\n\n')

flux = np.array(sp[0].data)
wave = np.ones(sp[0].header["NAXIS1"], dtype=float)

for i in np.arange(sp[0].header["NAXIS1"]):
    wave[i] = (
        sp[0].header["CRVAL1"]
        + (i - sp[0].header["CRPIX1"] + 1) * sp[0].header["CDELT1"]
    )
    # The list wave contains the wavelengths of the pixels.
# Close the fits-file:
sp.close()

frage1 = input(
    "Möchten Sie die Integrationsgrenzen zahlenmäßig eingeben oder \
per Mausklick (grafisch)? Geben Sie m oder g ein: "
)

# Plot the spectrum
fig = plt.figure()
plt.plot(wave, flux)
plt.xlabel("Wavelength [Angström]")
plt.ylabel("ADU")
plt.title("Spectrum " + filelist[0])
plt.grid(True)
print('\nBitte Plot auf den Integrationsbereich innerhalb 15 sek interaktiv einengen')
plt.pause(15)
print('Jetzt per Mausklick die Integrationsgrenzen setzen')


if frage1 == "g":
    # Interaktives Festlegen der Integrationsgrenzen
    # Falls die Darstellung des Spektrums interaktiv vergrößert wird die Eingabepunkte
    # jeweils mit der rechten Maustaste löschen, bis dann wirklich die linke
    # Seite der zu messenden Linie angeklickt wird.
    pts = []
    pts = np.asarray(plt.ginput(n=2, timeout=-1))
    plt.plot(pts[:, 0], pts[:, 1], "o", markersize=3)
    begin = pts[0, 0]
    end = pts[1, 0]
    print("Gewählter Integrationsbereich:", begin, " bis", end)
    print()

if frage1 == "m":
    # Eingabe der Integrationsgrenzen
    begin = float(
        input("Geben Sie die kurzwellige Wellenlänge-Integrationsgrenze ein: ")
    )
    end = float(
        input("Geben Sie die langwellige Wellenlänge-Integrationsgrenze ein: "))
    print()

# Initiieren von Arrays
EW = np.zeros(len(filelist))
JD = np.zeros(len(filelist))
linienminimum_flux = np.zeros(len(filelist))
linienminimum_wave = np.zeros(len(filelist))
fwhm = np.zeros(len(filelist))


# Abarbeiten der filelist
for i in range(len(filelist)):
    flux, header = fits.getdata(filelist[i], header=True)

    # Header Überprüfung
    if "JD-OBS" in header:
        JD[i] = float(header["JD-OBS"])
    elif 'JD' in header:
        JD[i] = float(header["JD"])
    elif "JD_OBS" in header:
        JD[i] = float(header["JD_OBS"])
    elif "BAS_MJD" in header:
        JD[i] = float(header["BAS_MJD"])
    else:
        print("Kein Beobachtungszeitpunkt im Header")
        break

    try:
        header["CRPIX1"]
    except:
        header["CRPIX1"] = 1

    # Generation of array with the wavelengths of the spectrum
    wave = np.ones(len(flux), dtype=float)
    for k in range(header["NAXIS1"]):
        wave[k] = (
            header["CRVAL1"]
            + (k - header["CRPIX1"] + 1) * header["CDELT1"])

    # Isolieren des Integrationsbereichs:
    ind = []
    for n in range(len(wave)):
        if wave[n] >= begin and wave[n] <= end:
            ind.append(n)
    intervallwave = wave[ind]
    intervallflux = flux[ind]

    fig = plt.figure()
    plt.plot(intervallwave, intervallflux, '-k',
             label='Originale Messung')

    intervallflux_begin = intervallflux[0]
    intervallflux_end = intervallflux[-1]
    # Steigungsbereechnung zur linearen Interpolation des Kontinuums (faktor)
    faktor = (intervallflux_end - intervallflux_begin) / len(intervallflux)

    # Renormierung des Integrationsbereichs auf das (lineare) Pseudokontinuum
    # gebildet aus begin und end (den Integrationsgrenzwellenlängen)
    ew = 0
    for p in range(len(intervallflux)):
        intervallflux[p] = intervallflux[p] / (intervallflux_begin + p * faktor)
        # Berechnung der EW
        dif = 1 - intervallflux[p]
        ew += dif * header["CDELT1"]
    EW[i] = ew

    # Erhöhung der Auflösung auf das Zehnfache durch Interpolation
    # Radial basis function (RBF) über das Spektrum im Intervall
    # smooth anpassen, 0. = Funktion geht durch alle Punkte, >0. = Ausgleich
    rbf = Rbf(intervallwave, intervallflux, smooth=1)  # smooth evtl. anpassen
    newwaveintervall = np.arange(
        intervallwave[0], intervallwave[-1], header['CDELT1'] / 10)
    newfluxinterpolated = rbf(newwaveintervall)
    linienminimum_flux[i] = newfluxinterpolated.min()
    linienminimum_wave[i] = newwaveintervall[newfluxinterpolated.argmin()]
    print(filelist[i] + ':')
    print('EW:', EW[i])
    print('Wellenlänge Linienminimum = {:.2f}'.format(
        linienminimum_wave[i]), 'Angström')
    print('Flux Linienminimum = {:.3f}'.format(linienminimum_flux[i]))

    # FWHM ermitteln
    hm = (1 + linienminimum_flux[i]) / 2
    for o in range(newfluxinterpolated.argmin(), 0, -1):
        if newfluxinterpolated[o] >= hm:
            links = newwaveintervall[o]
            break
    for o in range(newfluxinterpolated.argmin(), len(newfluxinterpolated), 1):
        if newfluxinterpolated[o] >= hm:
            rechts = newwaveintervall[o]
            break
    fwhm[i] = rechts - links
    print('FWHM = {:.2f}'.format(fwhm[i]), ' Angström\n')

    # Plotten der Linien
    plt.plot(intervallwave, intervallflux, color='b',
             label='renormiertes Spektrum',
             linewidth=1.)
    plt.plot(linienminimum_wave[i],
             linienminimum_flux[i], "or",
             label='Minimum')
    plt.title('Verwendeter Spektrumbereich von ' + filelist[i])
    plt.plot(newwaveintervall, newfluxinterpolated, color='r', linewidth=.6,
             label='interpoliertes Spektrum')
    plt.legend(loc='best')
    plt.grid(True)
    plt.pause(.1)

    if frage == 'j':
        plt.savefig(filelist[i] + '{:.2f}_'.format(begin) + '{:.2f}'.format(end)
                    + '.pdf')

    plt.close('all')

# Abspeichern der ermittelten Daten
ascii.write(
    [filelist, JD, EW, linienminimum_wave, linienminimum_flux, fwhm],
    'Wellenlängenbereich' +
    '{:.2f}_'.format(begin) + '{:.2f}'.format(end) + '.csv',
    names=["Spektrum ", 'JD', "EW", 'Wellenlaenge Minimum',
           'Flux Linienminimum', 'FWHM'],
    overwrite=True,
    format="csv",
)


# # Plotten der EWs
# fig = plt.figure()
# plt.stem(JD, -EW)
# plt.xlabel("JD")
# plt.ylabel("-EW in Angstroem")
# plt.title('EW ' + str(begin) + ' bis ' + str(end) + ' Angstroem')
# plt.grid(True)
# plt.savefig('EW_' + str(begin) + '_' + str(end) + '.pdf')

print('Ende des Programms')
