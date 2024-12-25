#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
newbinning_mittleresSpektrum.py

Liest einen Spektrenkatalog ein (fits), rebinnt die Spektren und erzeugt
tab-Spektren (tab-Tabelle mit den Spaltenbenennungen WAVE und FLUX) sowie
fits-Dateien, alle mit gleicher wählbarer Schrittweite und gleichem
Wellenlängenbereich. Die Art der Interpolation (linear oder per kubischem
Spline) ist durch auskommentieren wählbar. Außerdem wird ein gemitteltes
Spektrum berechnet und als fits sowie als tab-Tabelle abgespeichert mit den Spaltennamen
WAVE und FLUX. Das gemittelte Spektrum wird auch geplottet.

@author: lothar schanne
24.01.2022
"""

import numpy as np
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import glob
from specutils import Spectrum1D
from specutils.manipulation import (
    LinearInterpolatedResampler,
    SplineInterpolatedResampler,
)
from astropy import units as u
from astropy.visualization import quantity_support
quantity_support()

# Fileliste erstellen. Spektren in einem (Unter)Ordner.
files = input("Pfad und Name der Spektren (nutze wildcards) : ")
filelist = glob.glob(files)

# Alphabetisches Sortieren. Bei richtiger Namensgebung ergibt das eine
# zeitliche Ordnung
filelist.sort()

# Ausdruck der Liste
print("\nSpektrenliste: \n")
print(filelist, end='\n')
print("\nAnzahl der Spektren: ", len(filelist), "\n")

# Berechnung und Ausgabe des gemeinsamen Wellenlängenbereichs und der Schrittweiten
well_min = 0
well_max = 10000
step_min = 100
step_max = 0
for i in range(len(filelist)):
    sp = fits.open(filelist[i])
    crval = sp[0].header["CRVAL1"]
    if "CRPIX1" not in sp[0].header:
        sp[0].header["CRPIX1"] = 1
    crpix = sp[0].header["CRPIX1"]
    cdel = sp[0].header["CDELT1"]
    wave_erstesPixel = crval - cdel * (crpix - 1)
    if wave_erstesPixel > well_min:
        well_min = wave_erstesPixel
    if wave_erstesPixel + cdel * sp[0].header["NAXIS1"] < well_max:
        well_max = wave_erstesPixel + cdel * sp[0].header["NAXIS1"]
    if cdel < step_min:
        step_min = cdel
    if cdel > step_max:
        step_max = cdel
    sp.close()

print("\nGemeinsamer Wellenlängenbereich: ", well_min, well_max)
print("Minimale und maximale Schrittweite: ", step_min, step_max)

# zu wählende Schrittweite und Wellenlängenbereich
newstep = float(input("Geben Sie die gewünschte Schrittweite ein: "))
print("\nSpecification of the wavelength range to be transferred ")
a = float(input("Wellenlänge Begin: "))
b = float(input("Wellenlänge End: "))

pixzahl = int((b - a) / newstep)
DATA = np.zeros((len(filelist), pixzahl+1))

# Berechnung der neuen Spektren der Serie
for i in range(len(filelist)):
    sp = fits.open(filelist[i])
    crval = sp[0].header["CRVAL1"]
    if "CRPIX1" not in sp[0].header:
        sp[0].header["CRPIX1"] = 1
    crpix = sp[0].header["CRPIX1"]
    cdel = sp[0].header["CDELT1"]
    wave_erstesPixel = crval - cdel * (crpix - 1)
    aindex = int((a - wave_erstesPixel) / cdel)
    bindex = int((b - wave_erstesPixel) / cdel)
    newflux = sp[0].data[aindex:bindex]
    newflux = newflux * u.dimensionless_unscaled
    newwave = np.zeros(bindex - aindex)
    for k in range(len(newwave)):
        newwave[k] = wave_erstesPixel + (aindex + k) * cdel
    newwave = newwave * u.AA

    # new binning:
    input_spec = Spectrum1D(spectral_axis=newwave, flux=newflux)
    new_disp_grid = np.arange(a, b, newstep) * u.AA

    # Interpolation
    # lineare Interpolation:
    linear = LinearInterpolatedResampler()
    new_spec = linear(input_spec, new_disp_grid)
    # Interpolation per Spline:
    # spline = SplineInterpolatedResampler()
    # new_spec = spline(input_spec, new_disp_grid)

    filename = (
        filelist[i].rsplit(".")[0] + "_" + str(int(a)) +
        "_" + str(int(b)) + ".dat"
    )
    # ascii-Datei speichern
    ascii.write(
        [new_spec.spectral_axis, new_spec.flux],
        filename,
        overwrite=True,
        names=["WAVE", "FLUX"],
        format="tab",
    )
    # fits-Datei berechnen und speichern
    sp[0].header["CRVAL1"] = float(new_spec.spectral_axis[0].value)
    sp[0].header["CRPIX1"] = 1
    sp[0].header["NAXIS1"] = len(new_spec.spectral_axis)
    sp[0].header["CDELT1"] = newstep
    newfile = (
        filelist[i].rsplit(".")[0] + "_" + str(int(a)) +
        "_" + str(int(b)) + ".fit"
    )
    fits.writeto(
        newfile,
        new_spec.flux.value,
        sp[0].header,
        overwrite=True,
        output_verify="silentfix",
    )

    for j in range(len(new_spec.flux)):
        DATA[i, j] = new_spec.flux[j].value
    sp.close()

# Berechnung und Abspeichern des mittleren Spektrums
obj = input("Geben Sie den Objektnamen ohne blanks ein: ")
mean = np.zeros(len(DATA[0]) - 1)
for m in range(len(mean)):
    mean[m] = DATA[:, m].mean()
file = obj + "_" + str(int(a)) + "_" + str(int(b))
filename = file + "_mean.dat"
ascii.write(
    [new_spec.spectral_axis / u.AA, mean[:]],
    filename,
    overwrite=True,
    names=["WAVE", "FLUX"],
    format="tab",
)
filename = file + "_mean.fits"
fits.writeto(filename, mean[:], sp[0].header,
             overwrite=True, output_verify="silentfix")

# Plotten des mittleren Spektrums
fig = plt.figure(figsize=(14, 14))
plt.xlabel("Wavelength [Angstroem]")
plt.ylabel("Flux")
plt.title("Mittleres Spektrum " + filename)
plt.plot(new_spec.spectral_axis, mean[:])
plt.pause(1)
fig.savefig(filename.rstrip(".")[0] + "_mittleresSpektrum.png")


print('Zum beenden des Programms in das zuletzt geöffnete Diagramm klicken.')
plt.waitforbuttonpress(-1)
plt.close('all')
