#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Einlesen eines synthetischen Spektrums in Form einer Tabelle (wie als .spec
in der Pollux-Datenbank erhältlich).
In dem begleitenden Pollux-file Spektrum.txt sind Startwellenlänge und step
angegeben.
Berechnung eines wählbaren Ausschnitts (Pandas Dataframe mit 'bereich'
bezeichnet).
Geplottet wird die Spalte NFLUX = normierter Flux für das gesamte Spektrum
und für bereich.
Der Bereich kann mit einer Standardabweichung gefaltet werden (stddev in der
Einheit des step).
Der convolvierte Flux wird zusammen mit der Wellenlänge in einer zweispaltigen
ascii-Tabelle abgespeichert (spacer = tab, Spaltennamen WAVE und FLUX).

Stand 20180823
@author: lothar
"""
# from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.convolution import Gaussian1DKernel
from astropy.convolution import convolve
from astropy.io import ascii

file = input('Pfad und Filebezeichnung eingeben: ')
lambda_min = float(input('Eingabe der Startwellenlänge im Spektrum: '))
deltalambda = float(input('Eingabe der Schrittweite (Pixel) im Spektrum: '))
std = float(input('Standardabweichung für die Convolution in Vielfachen der\
 Schrittweite: '))

table = pd.read_fwf(file, names=['WAVE', 'AFLUX', 'NFLUX'], header=0)

wl_li = int(input('Wellenlängenbereich links: '))
wl_re = int(input('Wellenlängenbereich rechts: '))
# bereich =  table.query('(wl_li<WAVE) & (WAVE<wl_re)')
ind_li = int((wl_li-lambda_min) / deltalambda)
ind_re = int((wl_re-lambda_min) / deltalambda)
bereich = table[ind_li:ind_re]

#   Convolve mit astropy.convolution
kernel = Gaussian1DKernel(stddev=std)
con = np.array(bereich['NFLUX'])
convoluted = convolve(con, kernel, normalize_kernel=True,
                      boundary='extend')

# Grafik
plt.style.use('seaborn-whitegrid')
fig, ax = plt.subplots(3)
fig.set_size_inches(15, 15)
plt.xlabel('Wellenlänge [Angström]')
fig.suptitle('Spektrum '+file)
plt.grid(True)
ax[0].plot(table['WAVE'], table['NFLUX'])
ax[1].plot(bereich['WAVE'], bereich['NFLUX'])
ax[2].plot(bereich['WAVE'], convoluted)
plt.savefig(file+'.png')
plt.savefig(file+'.pdf')
plt.show()

# Speichern des convolved, NFLUX jetzt als FLUX
convol_file = pd.DataFrame(bereich['WAVE'], convoluted,
                           columns=['WAVE', 'NFLUX'])
name = file.rstrip('.spec')+'_convolved.dat'
ascii.write([bereich['WAVE'], convoluted], name, overwrite=True,
            names=['WAVE', 'FLUX'], format='tab')
