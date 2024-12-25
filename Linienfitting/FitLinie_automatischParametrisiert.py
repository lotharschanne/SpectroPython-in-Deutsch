#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Das Skript dient der Modellierung eines Absorptions- oder Emissionslinienprofils.
Das Profil kann eines oder zwei Maxima oder Minima haben.

Verschiedene Modelle (Gauss, Voigt, SplitLorentz für asymmetrische Profile)
können ausgewählt werden.

Das Skript lädt eine Serie auf das Kontinuum normierter Spektren im fits-Format
ein. Wenn in einer Serie das Linienprofil dopplerverschoben ist oder ganz anders
ist wie im ersten Spektrum kann es sein, dass das fitting nicht funktioniert.
Solche Spektren müssen dann einzeln behandelt werden.

Das erste Spektrum wird geplottet und man kann die zu modellierende Linie
für eine in Zeile 89 definierte und änderbare Zeit entweder per Mausklick
vergrößern und auswählen oder alternativ manuell den Wellenlängenbereich und die
Wellenlänge des Linienextremums eingeben.

Bei der grafischen Auswahl der Linie muss man vier mal mit der Maus in den plot
klicken, zuerst links der Linie (mit möglichst viel Kontinuum), dann rechts
(mit möglichst viel Kontinuum) und als drittes das ungefähre stärkere
Linienminimum/maximum und als viertes das schwächere Linienminimum/maximum
(oder wieder das Linienminimum/maximum).
Dadurch werden die betreffenden Wellenlängen und Flüsse ausgewählt.

Nach welchem Modell gefittet werden soll, ist wählbar. Zur Auswahl stehen
Gauss, Doppelgauss, Lorentz, Voigt, Doppelvoigt.
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
# Die Wartezeit, in der man das Spektrum vergrößern kann, evtl. anpassen:
plt.pause(10)


frage1 = input(
    "Möchten Sie die Spektrumgrenzen und die Linienminimamaximazahlenmäßig \
eingeben (m) oder per Mausklick (grafisch, g)? Geben Sie m oder g ein: "
)

if frage1 == "g":
    # Interaktives Festlegen der Wellenlängengrenzen
    # linke Seite mit möglichst viel Kontinuum, dann rechte Seite mit möglichst
    # viel Kontinuum, dann das stärkere Minimum/Maximum der Linie anklicken,
    # dann das schwächere zweite Minimum/Maximum anklicken.
    pts = []
    pts = np.asarray(plt.ginput(n=4, timeout=-1))
    plt.plot(pts[:, 0], pts[:, 1], "o", markersize=3)
    begin = pts[0, 0]
    end = pts[1, 0]
    extremum1 = pts[2, 0]
    extremum2 = pts[3, 0]
    print("Gewählter Bereich:", begin, " bis", end, ', Extremum1 =', extremum1,
          'Extremum2 =', extremum2)
    print()

if frage1 == "m":
    # Eingabe der Wellenlängengrenzen
    begin = float(
        input("Geben Sie die kurzwellige Wellenlängengrenze ein: ")
    )
    end = float(
        input("Geben Sie die langwellige Wellenlängengrenze ein: "))
    print()
    extremum1 = float(
        input("Geben Sie die ungefähre Wellenlänge des stärkeren Linienminimums/maximums ein: "))
    extremum2 = float(
        input("Geben Sie die ungefähre Wellenlänge des schwächeren Linienminimums/maximums ein: "))

plt.ioff()

modell = input('Welches Modell möchten Sie anwenden?\
 Gauss (1), Doppel-Gauss (2), Split-Lorentzian (3), Voigt (4), Doppel-Voigt (5)\
 Bitte angegebene Zahl wählen: ')

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
    pars = lin_mod.make_params(intercept=1, slope=0)

    # Modelle für die Linie:
    if modell == '1':
        # Gauss:
        gauss1 = GaussianModel(prefix='g1_')
        if pts[2,1] > 1:
            pars.update(gauss1.make_params(center=dict(value=extremum1),
                                       sigma=dict(value=5, min=0),
                                       amplitude=dict(value=pts[2,1], min=0.)))
        if pts[2,1] < 1:
            pars.update(gauss1.make_params(center=dict(value=extremum1),
                                       sigma=dict(value=5, min=0),
                                       amplitude=dict(value=pts[2,1]-1., max=0.)))

        mod = lin_mod + gauss1

    if modell == '2':
        # Doppelgauss, falls zweite Gausslinie überlappend vorhanden ist:
        gauss1 = GaussianModel(prefix='g1_')
        if pts[2,1] > 1:
            pars.update(gauss1.make_params(center=dict(value=extremum1),
                                       sigma=dict(value=5, min=0),
                                       amplitude=dict(value=pts[2,1], min=0.)))
        if pts[2,1] < 1:
            pars.update(gauss1.make_params(center=dict(value=extremum1),
                                       sigma=dict(value=5, min=0),
                                       amplitude=dict(value=pts[2,1]-1., max=0.)))
        gauss2 = GaussianModel(prefix='g2_')
        if pts[3,1] > 1:
           pars.update(gauss2.make_params(center=dict(value=extremum1),
                                      sigma=dict(value=2, min=0),
                                      amplitude=dict(value=pts[3,1], min=0.)))
        if pts[3,1] < 1:
           pars.update(gauss2.make_params(center=dict(value=extremum1),
                                      sigma=dict(value=2, min=0),
                                      amplitude=dict(value=pts[3,1]-1., max=0.)))
        mod = lin_mod + gauss1 + gauss2

    if modell == '3':
        # SplitLorentzianModell, für asymmetrische Linien gedacht
        lorentz_mod = SplitLorentzianModel(prefix='lor_')
        if pts[2,1] > 1:
            pars.update(lorentz_mod.make_params(center=dict(value=extremum1),
                                       sigma=dict(value=5, min=0),
                                       sigma_r=dict(value=10, min=0),
                                       amplitude=dict(value=pts[2,1], min=0.)))
        if pts[2,1] < 1:
            pars.update(lorentz_mod.make_params(center=dict(value=extremum1),
                                       sigma=dict(value=5, min=0),
                                       sigma_r=dict(value=10, min=0),
                                       amplitude=dict(value=pts[2,1]-1., max=0.)))
        mod = lin_mod + lorentz_mod

    if modell == '4':
        # Voigt-Modell
        voigt_mod = VoigtModel(prefix='voi_')
        if pts[2,1] > 1:
            pars.update(voigt_mod.make_params(center=dict(value=extremum1),
                                       sigma=dict(value=5, min=0),
                                       amplitude=dict(value=pts[2,1], min=0.)))
        if pts[2,1] < 1:
            pars.update(voigt_mod.make_params(center=dict(value=extremum1),
                                       sigma=dict(value=5, min=0),
                                       amplitude=dict(value=pts[2,1]-1., max=0.)))
        mod = lin_mod + voigt_mod

    if modell == '5':
        # Doppeltes Voigt-Modell
        voigt_mod1 = VoigtModel(prefix='voi_1')
        if pts[2,1] > 1:
            pars.update(voigt_mod1.make_params(center=dict(value=extremum1),
                                       sigma=dict(value=3, min=0),
                                       amplitude=dict(value=pts[2,1], min=0.)))
        if pts[2,1] < 1:
            pars.update(voigt_mod1.make_params(center=dict(value=extremum1),
                                       sigma=dict(value=3, min=0),
                                       amplitude=dict(value=pts[2,1]-1., max=0.)))
        voigt_mod2 = VoigtModel(prefix='voi_2')
        if pts[3,1] > 1:
            pars.update(voigt_mod2.make_params(center=dict(value=extremum2),
                                       sigma=dict(value=1, min=0),
                                       amplitude=dict(value=pts[3,1], min=0.)))
        if pts[3,1] < 1:
            pars.update(voigt_mod2.make_params(center=dict(value=extremum2),
                                       sigma=dict(value=1, min=0),
                                       amplitude=dict(value=pts[3,1]-1., max=0.)))
        mod = lin_mod + voigt_mod1 + voigt_mod2

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
    if modell == '1':
        axes[1].plot(x, comps['g1_'], '--', label='Gauss Komponente 1')
    if modell == '2':
        axes[1].plot(x, comps['g1_'], '--', label='Gauss Komponente 1')
        axes[1].plot(x, comps['g2_'], '--', label='Gauss Komponente 2')
    if modell == '3':
        axes[1].plot(x, comps['lor_'], '--', label='Lorentz Komponente')
    if modell == '4':
        axes[1].plot(x, comps['voi_'], '--', label='Voigt Komponente')
    if modell == '5':
        axes[1].plot(x, comps['voi_1'], '--', label='Voigt Komponente 1')
        axes[1].plot(x, comps['voi_2'], '--', label='Voigt Komponente 2')
    axes[1].legend()
    axes[1].set_xlabel("Wellenlänge in Angström", fontsize=10)
    plt.pause(.2)

    if graphspeichern == 'y':
        plt.savefig('Linefit_'+filelist[k] + '.pdf')

plt.show()


# Gefittete Modellparameter als Excel-Datei im Arbeitsverzeichnis abspeichern
ergebnis.to_excel('Fitergebnisse.xlsx', sheet_name='Fitergebnisse')


fr = input('Sind Sie fertig mit dem Betrachten der Grafiken und möchten das\
Programm beenden, dann drücken Sie die Taste e: ')
if fr == 'e':
    print('Ende des Programms')
    plt.close('all')
