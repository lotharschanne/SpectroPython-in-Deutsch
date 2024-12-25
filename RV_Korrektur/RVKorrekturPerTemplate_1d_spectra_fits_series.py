#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RVKorrektur_1d_spectra_fits_series.py

The script corrects a series of 1d-spectra in fits format of an object
by an per KK calculated RV [km/s].
Writes 1 fit-file of the RV-corrected spectra into the
working directory. The respective RV is noted in the the header of the
generated fit.

@author: Lothar Schanne
Stand 20211225
"""


from astropy.io import fits, ascii
import glob
from PyAstronomy import pyasl
import numpy as np
import matplotlib.pyplot as plt


# Create file list. Spectra in a (sub)folder.
files = input("Path and name of the spectra (use wildcards) : ")
filelist = glob.glob(files)

# Alphabetical sorting. If the name is correct (e.g. 20180922-xyz.fit),
# this results in a temporal order.
filelist.sort()

# Print the list
print("\nList of spectra: \n")
print("Number of spectra: ", len(filelist), "\n")

# Template auswählen
# Pfad und Name des templates
tfile = input('Pfad und Filebezeichnung des template eingeben: ')
#   Einlesen von Header und Daten (Flux vom template in tf,
#   Header des template in theader gespeichert)
tf, theader = fits.getdata(tfile, header=True)
print(
    'Flux-Minimum und -Maximum im Template [ADU]: ', tf.min(), '  ', tf.max())
#   Prüfung auf die nötigen header-Einträge
print('Ausgabe der zur Wellenlängenberechnung nötigen Headereinträge im template:')
if 'NAXIS' in theader:
    print('Dimension, NAXIS:                        ', theader['NAXIS'])
else:
    print('Das ist kein 1d-Spektrum !')
if 'NAXIS1' in theader:
    tnax = theader['NAXIS1']
    print('Anzahl der Werte (Abszisse), NAXIS1:     ', tnax)
else:
    print('NAXIS1 fehlt im header !')
if 'CRVAL1' in theader:
    tcrval = theader['CRVAL1']
    print('Anfangs-Wellenlänge, CRVAL1:             ', tcrval)
else:
    print('CRVAL1 fehlt im header !')
if 'CDELT1' in theader:
    tcdel = theader['CDELT1']
    print('Schrittweite der Wellenlänge, CDELT1:    ', tcdel)
else:
    print('CDELT1 fehlt im header !')
if 'CRPIX1' in theader:
    print('Referenzpixel, CRPIX1: ', theader['CRPIX1'])
else:
    theader['CRPIX1'] = 1
#   Erzeugen eines numpy-Arrays mit den Wellenlängen des template
tw = np.ones(tnax, dtype=float)
tcrval = tcrval + (1 - theader['CRPIX1']) * tcdel
for i in range(tnax):
    tw[i] = tcrval + i * tcdel
print('Wellenlängenbereich des template in Angström: ', tw[0], tw[-1])
print('\n')

RV = np.zeros(len(filelist))


# Processing the filelist
for i in range(len(filelist)):
    #   Einlesen von Header und Daten
    f, header = fits.getdata(filelist[i], header=True)
    print('\nFlux-Minimum und -Maximum im target ' +
          filelist[i], f.min(), '  ', f.max())
    #   Prüfung auf die nötigen header-Einträge
    if 'NAXIS' in header:
        print('Dimension, NAXIS:                        ', header['NAXIS'])
    else:
        print('Das ist kein 1d-Spektrum !')
    if 'NAXIS1' in header:
        nax = header['NAXIS1']
        print('Anzahl der Werte (Abszisse), NAXIS1:     ', nax)
    else:
        print('NAXIS1 fehlt im header !')
    if 'CRVAL1' in header:
        crval = header['CRVAL1']
        print('Anfangs-Wellenlänge, CRVAL1:             ', crval)
    else:
        print('CRVAL1 fehlt im header !')
    if 'CDELT1' in header:
        cdel = header['CDELT1']
        print('Schrittweite der Wellenlänge, CDELT1:    ', cdel)
    else:
        print('CDELT1 fehlt im header !')
    if 'CRPIX1' in header:
        print('Referenzpixel, CRPIX1: ', header['CRPIX1'])
    else:
        header['CRPIX1'] = 1
    #   Erzeugen eines numpy-Arrays mit den Wellenlängen des target
    w = np.ones(nax, dtype=float)
    for j in range(nax):
        crval = crval + (1 - header['CRPIX1']) * cdel
        w[j] = crval + j * cdel
    print('Wellenlängenbereich in Angström: ', w[0], w[-1])

    # Cross-Correlation ausführen
    # Das RV-range ist (Parameter 1) - bis (Parameter 2) km/s in
    # Schritten von (Parameter 3) km/s.
    # Die ersten und letzten (Parameter 4) Datenpunkte sind unberücksichtigt.
    # muss angepasst werden, damit die Spektren noch überlappen können bei der
    # maximalen Verschiebung (Parameter 2)
    rv, cc = pyasl.crosscorrRV(
        w, f, tw, tf, -200., 200., 0.1, mode='doppler', skipedge=10)
    # Maximum der cross-correlation function
    maxind = np.argmax(cc)
    print("Die Cross-correlation function ist maximal bei dRV = ", rv[maxind],
          " km/s.")
    if rv[maxind] > 0.0:
        print("Rotverschiebung des target gegenüber dem template.")
    else:
        print("Blauverschiebung des target gegenüber dem template")

    # Shift that spectrum
    flux_rv, wave_rv = pyasl.dopplerShift(
        w, f, -rv[maxind], edgeHandling="firstlast")

    header["CRVAL1"] = w[0]
    header['CRPIX1'] = 1

    newfile = filelist[i].rsplit(".", 1)[0] + "_RVcorrected_perTemplate.fit"
    header["RVCORR"] = (rv[maxind], "km/s, corrected")
    fits.writeto(newfile, flux_rv, header, overwrite=True,
                 output_verify="silentfix")

    RV[i] = rv[maxind]

    plt.figure(figsize=(10, 10))
    plt.title(filelist[i], size=10)
    plt.xlabel('Wellenlänge [Angström]')
    plt.ylabel('Flux auf das Kontinuum normiert')
    plt.plot(tw, tf, '-k', label='Template', linewidth=1.5)
    plt.plot(w, f, '-b', label='Original', linewidth=.5)
    plt.plot(w, flux_rv, '-r', label='an Template angepasst', linewidth=1.5)
    plt.legend(fontsize='x-small')
    plt.pause(1)
    plt.savefig(filelist[i]+'RVcorr.pdf', format='pdf')
    plt.close()


# Abspeichern der RV's als ascii-Datei
filename = filelist[i].rsplit(".", 1)[0]
ascii.write(
    [filelist, RV],
    filename + "_RVs_perTemplate.dat",
    overwrite=True,
    names=[
        "Spektrum",
        "RV"
    ],
    format="tab",
)
