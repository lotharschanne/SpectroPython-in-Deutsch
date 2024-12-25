#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Erzeugung einer Spektrenliste der 1d-Spektren im fits-Format in einem Ordner,
Korrektur/Ergänzung des Headers. Es können beliebig viele Headereinträge 
geändert oder um neue Headereinträge ergänzt werden.

Created on Fr Jul 19 14:48:03 2018
@author: lothar
"""

import glob
from astropy.io import fits


# filelist erstellen
# Pfad und Name der Spektren eingeben
files = input("Geben Sie Pfad und Namen zu den Spektren ein: ")
filelist = glob.glob(files)

# Sortierung
filelist.sort()

# Ausdruck der Liste
print("\nSpektrenliste: \n", filelist, "\n")
print("Anzahl der Spektren: ", len(filelist), "\n")

#   Header-Check erstes Spektrum
hdul = fits.open(filelist[0])
print("\nInfo zum ersten Spektrum der Liste:")
hdul.info()
print("\nKompletter Header des ersten Spektrums der Liste:")
print("\n", repr(fits.getheader(filelist[0], 0)))

Eingabe = input(
    "Möchten Sie für alle Spektren der Serie Headerdaten ändern oder ergänzen? \
Wenn nein 0 eingeben, wenn der Headereintrag eine Zahl ist: 1, wenn \
der Headereintrag ein String ist: 2 : "
)

k = 0
KW = []
Wert = []
while Eingabe == '1' or Eingabe == '2':
    KW.append(input("Eingabe des Header-Keywords: "))
    if Eingabe == "1":
        Wert.append(float(input("Eingabe des Wertes dazu: ")))
    elif Eingabe == "2":
        Wert.append(input("Eingabe des Wertes dazu: "))
    else:
        print("Falsche Eingabe")
        break
    Eingabe = input(
        "Weitere Headereinträge ändern? Wenn nein Eingabe von 0, Wenn ja, \
Eingabe von 1 oder 2: "
    )
    k += 1


for i in range(len(filelist)):
    g = fits.open(filelist[i])
    for j in range(len(KW)):
        g[0].header[KW[j]] = Wert[j]
    # fits.update(filelist[i], g[0].data, g[0].header)
    fits.writeto(filelist[i], g[0].data, g[0].header,
                 overwrite=True, output_verify='silentfix')
    g.close()
