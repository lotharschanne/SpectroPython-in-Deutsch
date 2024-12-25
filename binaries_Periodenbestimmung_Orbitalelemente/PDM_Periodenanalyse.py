#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Liest eine Datei (im csv- oder tab-Format) ein mit den beiden Spalten JD
und zugehöriger Wert (z.B. RV's), aber ohne Spaltenüberschriften.
Macht dann mit beiden Spalten eine Periodenanalyse (PDM-Analyse),
zeigt grafisch das Periodogramm
und gibt die Periode (in der Einheit Tage) aus.

Der Periodenbereich, für den die Berechnung durchgeführt werden soll, wird 
abgefragt.

Angelehnt an:
    https://pyastronomy.readthedocs.io/en/latest/pyTimingDoc/pyPDMDoc/pdm.html

Created on Thu Sep  7 16:00:55 2023

@author: lothar
"""


from astropy.io import ascii
import matplotlib.pylab as plt
from PyAstronomy.pyTiming import pyPDM
from astropy.timeseries import TimeSeries
from astropy import time
from astropy import units as u

plt.switch_backend('Qt5Agg')  # Interaktives GrafikBackend einschalten

filename = input("Geben Sie den Name des Datenfiles ein: ")
ts = ascii.read(filename)
ts_JD = ts.columns[0]
ts_RV = ts.columns[1]

print('Geben Sie nachfolgend den Beginn und das Ende des Periodenbereichs ein,\
der überprüft werden soll.')
pa = float(input('Geben Sie das minimale P in der Einheit Tage ein: '))
pe = float(input('Geben Sie das maximale P in der Einheit Tage ein: '))

# Get a ``scanner'', which defines the frequency interval to be checked.
# Alternatively, also periods could be used instead of frequency.
# Bitte anpassen !!!!!!!!!!!!!!!
S = pyPDM.Scanner(minVal=pa, maxVal=pe, dVal=0.001, mode="period")

# Carry out PDM analysis. Get frequency array
# (f, note that it is period, because the scanner's
# mode is ``period'') and associated Theta statistic (t).
P = pyPDM.PyPDM(ts_JD, ts_RV)
# Use 10 phase bins and 3 covers (= phase-shifted set of bins).
# Den ersten Parameter <= Datenpunkte/2 wählen, den zweiten auch variieren:
f1, t1 = P.pdmEquiBinCover(3, 3, S)

# Show the result
plt.figure(facecolor='white')
plt.title("Result of PDM analysis")
plt.xlabel("Frequency")
plt.ylabel("Theta")
plt.plot(f1, t1, 'bp-')
plt.pause(.1)


Periode = f1[t1.argmin()] * u.d
print('Periode = ', Periode.round(6))


# Time Series Objekt erstellen
ts_JD = time.Time(ts.columns[0], format='jd')
serie = TimeSeries(time=ts_JD)
# RV's einfügen
serie.add_column(ts_RV)

# Folding (Phasenplot)
a, b = P.phase(ts_JD, Periode)
fig = plt.figure()
plt.plot(a, ts_RV, 'ko', markersize=3)
plt.xlabel('Time (days)')
plt.ylabel('RV [km/s]')
plt.title(filename + '\nPDM-Phasendiagramm mit der Periode ' +
          str(Periode.round(6))+' berechnet')
plt.pause(.1)
plt.savefig(filename + '_Phasenplot_PDM', format='pdf')
print('Zum beenden des Programms in das zuletzt geöffnete Diagramm (Figure 2) klicken')
plt.waitforbuttonpress(-1)
plt.close('all')
