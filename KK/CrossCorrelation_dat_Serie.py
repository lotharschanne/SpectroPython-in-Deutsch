#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Kreuzkorrelation
abgeleitet von einem Beispiel in PyAstronomy
https://www.hs.uni-hamburg.de/DE/Ins/Per/Czesla/PyA/PyA/pyaslDoc/aslDoc/
crosscorr.html

Es wird eine Kreuzkorrelation einer Serie von target-Spektren bzgl. eines
template-Spektrums durchgeführt. Beide liegen als ascii-Datei, Format .dat,
(2 Spalten, erste Wellenlänge, zweite flux, mit oder ohne Überschrift) vor.

Die berechneten RV's werden ausgedruckt und in einer Datei gespeichert. Die 
Kreuzkorrelationsfunktionen weden geplottet und die Grafiken werden 
abgespeichert.

Stand 20220408

@author: Lothar Schanne
"""


from PyAstronomy import pyasl
import numpy as np
import matplotlib.pylab as plt
from astropy.io import ascii
from astropy.table import Table
import glob
from linetools.spectra.xspectrum1d import XSpectrum1D


obj = input("Geben Sie die Überschrift für die Grafiken ein (ohne blanks): ")

# Template auswählen
# Pfad und Name des templates
template_name = input("Pfad und Filebezeichnung des template eingeben: ")
template = np.loadtxt(template_name, skiprows=1)
print("Minimum und Maximum im Template [ADU]: ",
      template.min(), "  ", template.max())

w = template[:, 0]
f = template[:, 1]

dt = (w[-1] - w[0]) / len(w)

template_newdata, dt_template = pyasl.binningx0dt(
    w, f, dt=dt, x0=w[0], useBinCenter=True
)

# print("Template Beginn und Ende: ", template_newdata[0, 0], template_newdata[-1, 0])
print("Template Beginn und Ende: ", w[0], w[-1])

# Zu korrelierendes Spektrum (target) auswählen
# Pfad und Name des target
# Create file list. Spectra in a (sub)folder.
files = input("Path and name of the target spectra (use wildcards) : ")
filelist = glob.glob(files)

# Sort alphabetically. If the spectrum files are named correctly, this results
# in a temporal order.
filelist.sort()

# Printout of the list for control purposes.
print("\nList of target spectra:")
print(filelist)
print("\nNumber of target spectra: ", len(filelist), "\n")
print("Bitte warten. Berechnung läuft")

RV = np.zeros(len(filelist))
jd = np.zeros(len(filelist))

for i in range(len(filelist)):
    target = np.loadtxt(filelist[i], skiprows=1)

    #   Erzeugen je eines numpy-Arrays mit den Wellenlängen und Flux des target
    tw = target[:, 0]
    tf = target[:, 1]
    binnumber = int((tw[-1] - tw[0]) / dt_template)

    target_newdata, dt_target = pyasl.binningx0dt(
        tw, tf, x0=tw[0], nbins=binnumber, useBinCenter=True
    )

    print(filelist[i])
    print("Target Beginn und Ende:",
          target_newdata[0, 0], target_newdata[-1, 0])

    # Cross-Correlation ausführen.
    # Das RV-range ist (Parameter 1) - bis (Parameter 2) km/s in
    # Schritten von (Parameter 3) km/s.
    # Die ersten und letzten (Parameter 4) Datenpunkte sind unberücksichtigt.
    # Muss angepasst werden, damit die Spektren noch überlappen können bei der
    # maximalen Verschiebung (Parameter 2)
    rv, cc = pyasl.crosscorrRV(
        target_newdata[:, 0],
        target_newdata[:, 1],
        w,
        f,
        -200.0,
        200.0,
        0.1,
        mode="doppler",
        skipedge=100,
    )

    # Maximum der cross-correlation function finden
    maxind = np.argmax(cc)
    RV[i] = rv[maxind]

    # Plot CroCo-Funktion
    fig = plt.figure()
    plt.plot(rv, cc, "b-")
    plt.plot(rv[maxind], cc[maxind], "ro")
    plt.title(obj + " : Kreuzkorrelationsfunktion " + filelist[i])
    plt.grid(True)
    plt.xlabel("km/s")
    plt.pause(.1)
    plt.savefig(filelist[i]+'_CroCo.png')


print("Spektrum         ", 'JD', "RV [km/s]")
for i in range(len(filelist)):
    print(filelist[i], jd[i], RV[i])

data = Table([filelist, jd, RV], names=["Spektrum", 'JD', "RV"])
ascii.write(data, "RV_Liste_" + obj + ".dat", overwrite=True, format="tab")

print('Zum beenden des Programms in das zuletzt geöffnete Diagramm klicken')
plt.waitforbuttonpress(-1)
plt.close('all')

print("Ende des Programs")
