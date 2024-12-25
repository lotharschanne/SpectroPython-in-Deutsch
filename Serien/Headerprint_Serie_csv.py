#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Liest für eine 1d-Spektrenserie im fits-Format die Headerdaten ein und
schreibt sie in getrennte ascii-Dateien (.csv).

Created on Fr Jul 19 14:48:03 2018
@author: Lothar Schanne
"""

import glob
from astropy.io import fits
import pandas as pd

# filelist erstellen
# Pfad und Name der Spektren eingeben
files = input("Geben Sie Pfad und Namen zu den Spektren ein: ")
filelist = glob.glob(files)

# Sortierung
filelist.sort()

# Ausdruck der Liste
print("\nSpektrenliste: \n", filelist, "\n")
print("Anzahl der Spektren: ", len(filelist), "\n")

# a = fits.open(filelist[0])
# header = a[0].header
# headerdict = dict(header)
# file = filelist[0].rsplit(".")[0]
# pd_ser = pd.Series(headerdict.values(), index=headerdict.keys())
# pd_fra = pd.DataFrame({file: pd_ser})

# f = open(filelist[0].rstrip(".fit") + "_Header.csv", "w")
# for i in headerdict:
#     f.write(str(i) + "," + str(headerdict[i]) + "\n")
# f.close()


for f in range(0, len(filelist)):
    print('Spektrum: ', filelist[f])
    with fits.open(filelist[f]) as hdul:
        hdul.verify('fix')
        header = hdul[0].header
    # g = fits.open(filelist[f], ignore_missing_simple=True)
    # header = g[0].header
    headerdict = dict(header)
    fi = open(filelist[f].rsplit(".")[0] + "_Header.csv", "w")
    for i in headerdict:
        fi.write(str(i) + "," + str(headerdict[i]) + "\n")
    fi.close()
    # Falls alle Header den gleichen Keywords-Satz haben, kann ein Pandas_Frame
    # erzeugt werden. Dann die folgende Zeile auskommentieren.
    # pd_fra[filelist[f].rstrip(".fit")] = headerdict.values()
