#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Das Skript dient der Gauss-Modellierung eines Emissionslinienprofils, das wie bei
beta Lyrae aus einer gaussförmigen Emissionslinie mit einer darin enthaltenen
zentralen Absorption besteht.

Das Skript lädt eine Serie auf das Kontinuum normierter Spektren im fits-Format
ein. Wenn in einer Serie das Linienprofil dopplerverschoben wird oder ganz anders
ist wie im ersten Spektrum kann es sein, dass das fitting nicht funktioniert.
Solche Spektren müssen dann einzeln behandelt werden.

Dann wird das erste Spektrum geplottet und man kann die zu modellierende Linie
für eine in Zeile 89 definierte und änderbare Zeit entweder per Mausklick
vergrößern und auswählen.

Bei der grafischen Auswahl der Linie muss man vier mal mit der Maus in den plot
klicken, zuerst links der Linie (mit möglichst viel Kontinuum), dann rechts
(mit möglichst viel Kontinuum) und als drittes das ungefähre Linienmaximum und
als viertes auf das das Minimum der eingelagerten Absorption.
Dadurch werden die betreffenden Wellenlängen und Flüsse ausgewählt.

Die Ergebnisse des Fittings werden in einer Excel-Datei gespeichert.

Die Grafiken können auf Wunsch gespeichert werden.

20240326

@author: lothar schanne
"""

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from lmfit.models import LinearModel, GaussianModel, SplitLorentzianModel, VoigtModel
import glob
import pandas as pd


plt.switch_backend('Qt5Agg')  # Interaktives GrafikBackend einschalten
plt.ion()

# Fileliste erstellen. Spektren in einem (Unter)Ordner.
files = input("Pfad und Name der Spektren (nutze wildcards) : ")
filelist = glob.glob(files)

# Alphabetisches Sortieren. Bei richtiger Namensgebung ergibt das eine
# zeitliche Ordnung
filelist.sort()

# Ausdruck der Liste
print("\nSpektrenliste:")
print(filelist)
print("\nAnzahl der Spektren: ", len(filelist), "\n")

# Grafik erstes Spektrum
sp = fits.open(filelist[0], ignore_missing_end=True)
# print('\n\nHeader of the spectrum :\n\n', sp[0].header, '\n\n')

flux = np.array(sp[0].data)
wave = np.ones(sp[0].header["NAXIS1"], dtype=float)

if 'CRPIX1' not in sp[0].header:
    sp[0].header['CRPIX1'] = 1

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
print('Jetzt die Linie aus dem Spektrum heranzoomen, dann warten')
# Die Wartezeit, in der man das Spektrum vergrößern kann, evtl. anpassen:
plt.pause(10)

print('Jetzt können die Punkte in der Grafik angeklickt werden.')

# Interaktives Festlegen der Wellenlängengrenzen
# linke Seite mit möglichst viel Kontinuum, dann rechte Seite mit möglichst
# viel Kontinuum, dann das Maximum der Linie anklicken,
# dann das Absorptionsminimum anklicken.
pts = []
pts = np.asarray(plt.ginput(n=4, timeout=-1))
plt.plot(pts[:, 0], pts[:, 1], "o", markersize=3)
begin = pts[0, 0]
end = pts[1, 0]
extremum1 = pts[2, 0]
extremum2 = pts[3, 0]
print("Gewählter Bereich:", begin, " bis", end, ', Maximum =', extremum1,
          'Wellenlänge der Absorption =', extremum2)
print()

plt.ioff()

graphspeichern = input('Möchten Sie die Grafiken speichern? Dann y eingeben:')

# Abarbeiten der filelist
for k in np.arange(len(filelist)):
    print()
    print(filelist[k])
    sp = fits.open(filelist[k], ignore_missing_end=True)
    # Auslesen des Beobachtungszeitpunktes als JD, Header Überprüfung
    if "JD-OBS" in sp[0].header:
        JD = float(sp[0].header["JD-OBS"])
    elif 'JD' in sp[0].header:
        JD = float(sp[0].header["JD"])
    elif "JD_OBS" in sp[0].header:
        JD = float(sp[0].header["JD_OBS"])
    elif "BAS_MJD" in sp[0].header:
        JD = float(sp[0].header["BAS_MJD"])
    else:
        print("Kein Beobachtungszeitpunkt im Header")
        JD = 0
        break

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

    # Ausschneiden des Linienbereichs:
    for n in np.arange(len(wave)):
        if wave[n] <= begin and wave[n + 1] > begin:
            begin_n = n
            break
    for n in np.arange(len(wave)):
        if wave[n] <= end and wave[n + 1] > end:
            end_n = n
            break

    x = wave[begin_n:end_n]
    y = flux[begin_n:end_n]


# Modelldefinitionen:

    # Modell für das Kontinuum:
    lin_mod = LinearModel(prefix='lin_')
    pars = lin_mod.make_params(intercept=1., slope=0.)

    # Doppelgauss:
    gauss1 = GaussianModel(prefix='g1_')
    pars.update(gauss1.make_params(center=dict(value=extremum1),
                                       sigma=dict(value=5., min=0.),
                                       amplitude=dict(value=pts[2,1])))
    gauss2 = GaussianModel(prefix='g2_')
    pars.update(gauss2.make_params(center=dict(value=extremum2),
                                       sigma=dict(value=1., min=0.),
                                       amplitude=dict(value=-(pts[2,1]-pts[3,1]), max=0.)))
    mod = lin_mod + gauss1 + gauss2



    # Fitten des Modells
    init = mod.eval(pars, x=x)
    out = mod.fit(y, pars, x=x)

    print(out.fit_report(correl_mode='list'))

    if k == 0:
        ergebnisse = pd.Series(out.values)
        ergebnisse = pd.concat([ergebnisse, pd.Series({'JD': JD})], axis=0)
        ergebnis = pd.DataFrame({filelist[k]: ergebnisse})

    if k > 0:
        ergebnisse = pd.Series(out.values)
        ergebnisse = pd.concat([ergebnisse, pd.Series({'JD': JD})], axis=0)
        ergebnis = pd.concat(
            [ergebnis, pd.DataFrame({filelist[k]:ergebnisse})], axis=1)

    # Plotten
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))

    axes[0].plot(x, y)
    axes[0].plot(x, init, '--', label='initial fit')
    axes[0].plot(x, out.best_fit, '-', label='best fit')
    axes[0].legend()
    axes[0].set_xlabel("Wellenlänge in Angström", fontsize=10)

    comps = out.eval_components(x=x)

    axes[1].plot(x, y)
    axes[1].plot(x, comps['lin_'], '--', label='Kontinuum Komponente')
    axes[1].plot(x, comps['g1_'], '--', label='Gauss Komponente 1')
    axes[1].plot(x, comps['g2_'], '--', label='Gauss Komponente 2')
    axes[1].legend()
    axes[1].set_xlabel("Wellenlänge in Angström", fontsize=10)
    plt.pause(.2)

    if graphspeichern == 'y':
        plt.savefig('Linefit_Gauss_'+filelist[k] + '.pdf')

plt.show()


# Gefittete Modellparameter als Excel-Datei im Arbeitsverzeichnis abspeichern
ergebnis.to_excel('Fitergebnisse.xlsx', sheet_name='Fitergebnisse')


fr = input('Sind Sie fertig mit dem Betrachten der Grafiken und möchten das\
Programm beenden, dann drücken Sie die Taste e: ')
if fr == 'e':
    print('Ende des Programms')
    plt.close('all')
