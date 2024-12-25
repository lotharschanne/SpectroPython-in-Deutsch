#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Das Skript liest eine Serie von fits-1d-Spektren ein und überprüft, ob im Header
ein Eintrag für das Beobachtungsdatum (in JD) vorhanden ist.
Die Befunde werden je Spektrum ausgedruckt.

20231122

@author: lothar
"""

import glob
from astropy.io import fits
from astropy.time import Time


# filelist erstellen
# Pfad und Name der Spektren eingeben
files = input("Geben Sie Pfad und Namen zu den Spektren ein: ")
filelist = glob.glob(files)

# Sortierung
filelist.sort()

# Ausdruck der Liste
print("\nSpektrenliste: \n", filelist, "\n")
print("Anzahl der Spektren: ", len(filelist), "\n")


# Abarbeitung der Spektrenliste
for i in range(len(filelist)):
    g = fits.open(filelist[i])
    print()
    print(filelist[i])

    if 'JD' in g[0].header:
        t = g[0].header['JD']
        dat = Time(t, format='jd')
        print('JD gefunden: ')
        print('Aufnahmedatum: ', t, ' = ', dat.strftime('%Y %b %d %H:%M:%S'))

    if 'JD' not in g[0].header:
        if 'OBSDATE' in g[0].header:
            print('OBSDATE gefunden:')
            t = [0].header['OBSDATE']
            dat = Time(t, format='jd')
            print('Aufnahmedatum: ', t, ' = ',
                  dat.strftime('%Y %b %d %H:%M:%S'))

        if 'JD-OBS' in g[0].header:
            print('JD-OBS gefunden:')
            t = g[0].header['JD-OBS']
            dat = Time(t, format='jd')
            print('Aufnahmedatum: ', t, ' = ',
                  dat.strftime('%Y %b %d %H:%M:%S'))

        if 'JD-MID' in g[0].header:
            print('JD-MID gefunden:')
            t = g[0].header['JD-MID']
            dat = Time(t, format='jd')
            print('Aufnahmedatum: ', t, ' = ',
                  dat.strftime('%Y %b %d %H:%M:%S'))

        if 'MID-HJD' in g[0].header:
            print('MID-HJD gefunden:')
            t = g[0].header['MID-HJD']
            dat = Time(t, format='jd')
            print('Aufnahmedatum: ', t, ' = ',
                  dat.strftime('%Y %b %d %H:%M:%S'))

        if 'DATE-OBS' in g[0].header:
            print('DATE-OBS gefunden:')
            t = g[0].header['DATE-OBS']
            try:
                dat = Time(t, format='isot')
                print('Aufnahmedatum: ', t, ' = ',
                      dat.strftime('%Y %b %d %H:%M:%S'))
            except:
                print('Aufnahmedatum: ', t)

        if 'DATE-BEG' in g[0].header:
            print('DATE-BEG gefunden:')
            t = g[0].header['DATE-BEG']
            pdat = Time(t, format='jd')
            print('Aufnahmedatum: ', t, ' = ',
                  dat.strftime('%Y %b %d %H:%M:%S'))

        if "MJD-OBS" in g[0].header:
            print('MJD-OBS gefunden:')
            t = g[0].header["MJD-OBS"]
            dat = Time(t, format='jd')
            print('Aufnahmedatum: ', t, ' = ',
                  dat.strftime('%Y %b %d %H:%M:%S'))

        if "BAS-MJD" in g[0].header:
            print('BAS-MJD gefunden:')
            t = g[0].header["BAS-MJD"]
            dat = Time(t, format='jd')
            print('Aufnahmedatum: ', t, ' = ',
                  dat.strftime('%Y %b %d %H:%M:%S'))
        else:
            print('Kein Beobachtungszeitpunkt im Header gefunden !')

    # dat = Time(t, format='jd')
    # print('Aufnahmedatum: ', t, ' = ', dat.strftime('%Y %b %d %H:%M:%S'))
    g.close()
print('Ende des Programms')
