#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Berechnet für eine Zeitserie auf das Kontinuum normierter Spektren im
fits-Format die Äquivalentweiten mehrerer zusammenhängender Bereiche einer Linie.
Eingabe der Integrationsgrenzen grafisch-interaktiv.

Das erste Spektrum wird geplottet und für 15 Sekunden angezeigt. In diesem
Zeitraum ist das Grafikfenster interaktiv geschaltet, so dass man das Spektrum
vergrößern kann und den Integrations-Wellenlängenbereich für die EW-Berechnung
optisch aussuchen kann. Dieses Reaktionszeitfenster von 15 Sekunden kann in
Zeile 73 geändert werden.

Nach Auswahl des darzustellenden Wellenlängenbereichs für alle Spektren der Serie
wird nach der Anzahl der Linienbereiche, für die die EW bestimmt werden soll,
gefragt. Danach werden die einzelnen Spektren der Serie grafisch dargestellt und
es wird verlangt, dass die Grenzen der getrennt auszuwertenden Bereiche per 
Mausklick auf die Grenzen eingegeben werden sollen. Dadurch werden die 
Wellenlängen der Integrationsgrenzen individuell für jedes Spektrum der Serie 
ermittelt.
Am Ende werden die ermittelten EW's in eine ascii-Datei geschrieben.

Version 20230524

@author: lothar schanne
"""

import numpy as np
from astropy.io import fits, ascii
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

# Plot the spectrum
fig = plt.figure()
plt.plot(wave, flux)
plt.xlabel("Wavelength [Angström]")
plt.ylabel("ADU")
plt.title("Spectrum " + filelist[0])
plt.grid(True)

print('Bitte  das Grafikfenster bearbeiten (die Linie heraus \
vergrössern) und dann warten')
plt.pause(15)
print('Jetzt mit 2 Mausklicks den Wellenlängenbereich der Linie auswählen für \
alle Spektren gemeinsam')

pts = []
pts = np.asarray(plt.ginput(n=2, timeout=-1))
plt.plot(pts[:, 0], pts[:, 1], "o", markersize=3)
bereichsanfang = pts[0, 0]
bereichsende = pts[-1, 0]

anzahl = 1 + int(input('Geben Sie die gewünschte Anzahl der getrennt zu \
integrierenden Linienbereiche ein: '))


EW = np.zeros([len(filelist), anzahl-1])
JD = np.zeros(len(filelist))

# Abarbeiten der filelist
for k in np.arange(len(filelist)):
    sp = fits.open(filelist[k], ignore_missing_end=True)
    # Header Überprüfung
    if "JD-OBS" in sp[0].header:
        JD[k] = float(sp[0].header["JD-OBS"])
    elif 'JD' in sp[0].header:
        JD[k] = float(sp[0].header["JD"])
    elif "JD_OBS" in sp[0].header:
        JD[k] = float(sp[0].header["JD_OBS"])
    elif "BAS_MJD" in sp[0].header:
        JD[k] = float(sp[0].header["BAS_MJD"])
    else:
        print("Kein Beobachtungszeitpunkt im Header")
        JD[k] = 0

    try:
        sp[0].header["CRPIX1"]
    except:
        sp[0].header["CRPIX1"] = 1

    # Generation of arrays with the wavelengths and fluxes of the spectrum
    flux = np.array(sp[0].data)
    wave = np.ones(sp[0].header["NAXIS1"], dtype=float)

    for i in np.arange(sp[0].header["NAXIS1"]):
        wave[i] = (
            sp[0].header["CRVAL1"]
            + (i - sp[0].header["CRPIX1"] + 1) * sp[0].header["CDELT1"]
        )

    # Close the fits-file
    sp.close()

    WAVE = []
    FLUX = []

    for n in np.arange(len(wave)):
        if wave[n] >= bereichsanfang and wave[n] <= bereichsende:
            WAVE.append(wave[n])
            FLUX.append(flux[n])

    # Plot the spectrum
    fig = plt.figure()
    plt.plot(WAVE, FLUX)
    plt.xlabel("Wavelength [Angström]")
    plt.ylabel("ADU")
    plt.title("Spectrum " + filelist[k])
    plt.grid(True)

    # Interaktives Festlegen der Integrationsgrenzen
    # durch anklicken des Spektrums von links nach rechts
    pts = []
    pts = np.asarray(plt.ginput(n=anzahl, timeout=-1))
    plt.plot(pts[:, 0], pts[:, 1], "o", markersize=3)
    begin = pts[0, 0]
    end = pts[-1, 0]
    print("Gewählter Integrationsbereich:", begin, " bis", end)

    for m in range(len(pts)-1):
        for n in np.arange(len(WAVE)):
            if WAVE[n] <= pts[m, 0] and WAVE[n + 1] > pts[m, 0]:
                begin_n = n
                break
        for n in np.arange(len(WAVE)):
            if WAVE[n] <= pts[m+1, 0] and WAVE[n + 1] > pts[m+1, 0]:
                end_n = n
                break
        ew = 0
        for o in np.arange(end_n - begin_n):
            dif = 1 - FLUX[begin_n + o]
            ew += dif
        EW[k, m] = ew * sp[0].header["CDELT1"]

    print(filelist[k], "   EW =", EW[k])
    print()


ascii.write(
    [filelist, JD, EW],
    'EW_' + str(int(begin)) + "_" + str(int(end)) + '.ecsv',
    names=["Spektrum", 'JD', "EW"],
    overwrite=True,
    format="ecsv",
)


print('Zum beenden des Programms in das zuletzt geöffnete Diagramm klicken')
plt.waitforbuttonpress(-1)
plt.close('all')
