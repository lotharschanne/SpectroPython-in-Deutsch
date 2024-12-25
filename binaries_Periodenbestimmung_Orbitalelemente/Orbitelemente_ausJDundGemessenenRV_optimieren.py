#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Das Skript liest eine ascii-Datei mit zwei Spalten ohne Überschrift mit den
Zeitpunkten (JD, erste Spalte) und den RV's einer binary-Komponente
(RV, zweite Spalte) ein.
Das Skript berechnet mit dem als bekannt vorausgesetzten (unabhängig bestimmten)
Orbitalparameter P den Phasenplot. Für die übrigen Orbitalelemente reichen im
allgemeinen die Schätzwerte der Zeilen 35 bis 38 aus.
Dann werden die restlichen Orbitalparameter K1, e, w1 und gamma mittels einer
Kurvenfitting-Routine optimiert, damit die theoretischen RV's berechnet und im
Phasenplot zusätzlich zum Vergleich mit den gemessenen RV's dargestellt.
Die optimierten Orbitalparameter K1, e, w1, und gamma und ihre
Standardwabweichungen werden ausgedruckt.

20231007
@author: lothar schanne
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Die RV-Tabelle im ascii-Format hat KEINE Spaltenüberschriften
table = input('Geben Sie den Dateinamen mit den RV-Daten ein: ')
data = pd.read_table(table, names=['JD', 'RV'])

P = float(input('Geben Sie die zu verwendende Periode [in Tagen] ein: '))
# # Anfangs-Orbitalparameter des Binärsystems (Schätzwerte)
# omega (w1) in Grad,
# T in JD, gamma (Systemgeschwindigkeit in km/s) sowie K1 in km/s,
gamma = 0
K1 = 50
e = .5
w1 = 50
T = data.JD[0]


# Umrechnung der Zeitpunkte relativ zu To und P (für Phasenplot):
JD = data.JD
zp = (JD - T) / P - (JD - T) // P  # zp = Phase

# Grafik
plt.plot(zp, data.RV, 'or', markersize=3., label='gemessene RV')


def pipe(zp, K1, e, w1, gamma):
    # Ermittlung der exzentrischen Anomalien E durch Lösung per Newtonverfahren
    def newton(E):
        return E - (E - args1['e'] * np.sin(E) - args1['M']) / (1.0 - args2['e'] * np.cos(E))

    E = np.ones(len(zp))
    v = np.ones(len(zp))
    vr1 = np.ones(len(zp))

    for n in range(len(zp)):
        rv = data.RV[n]
        while n >= 0:  # Iterative Berechnung der exzentrischen Anomalie E
            args1 = {'e': e, 'M': zp[n] * 2 * np.pi}
            args2 = {'e': e, 'M': zp[n] * 2 * np.pi}
            alt = E[n]
            E[n] = newton(alt)
            if abs(alt - E[n]) < 10e-12:
                break

        # Berechnung der wahren Anomalien v:
        v[n] = 2 * np.arctan(np.tan(E[n] / 2) *
                             np.sqrt(1 + e) / np.sqrt(1 - e))
        # Berechnung der beiden Radialgeschwindigkeiten:
        vr1[n] = gamma + K1 * (e * np.cos(w1 * (np.pi/180))
                               + np.cos(w1 * (np.pi/180) + v[n]))
    return vr1


# Optimierung der Orbitalparameter K1, e, w1 und gamma, bounds bitte anpassen,
# wenn sinnvoll
popt, pcov = curve_fit(pipe, zp.to_numpy(),
                       data.RV.to_numpy(),
                       bounds=([0., 0., 0., -50.], [150., 1., 360., 50.]))

print('Optimierte Orbitalparameter K1, e, w1 und gamma:\n', popt)
asini = (popt[0] * np.sqrt(1 - popt[1]**2) * P * 24*60*60) / 2 * np.pi
print('asini =', asini, 'km')
stdabw = np.sqrt(np.diag(pcov))
print('Standardabweichungen der optimierten Parameter K1, e, w1 und gamma:\n',
      stdabw)


def pipe2(zp, popt):  # Berechnung der RV's mit den optimierten Orbitalparametern

    def newton(E):
        return E - (E - args1['e'] * np.sin(E) - args1['M']) / (1.0 - args2['e'] * np.cos(E))

    E = np.ones(len(zp))
    v = np.ones(len(zp))
    vr1 = np.ones(len(zp))

    for n in range(len(zp)):
        while n >= 0:  # Iterative Berechnung der exzentrischen Anomalie E
            args1 = {'e': popt[1], 'M': zp[n] * 2 * np.pi}
            args2 = {'e': popt[1], 'M': zp[n] * 2 * np.pi}
            alt = E[n]
            E[n] = newton(alt)
            if abs(alt - E[n]) < 10e-12:
                break

        # Berechnung der wahren Anomalien v:
        v[n] = 2 * np.arctan(np.tan(E[n] / 2) *
                             np.sqrt(1 + popt[1]) / np.sqrt(1 - popt[1]))

        # Berechnung der Radialgeschwindigkeit:
        vr1[n] = popt[3] + popt[0] * (popt[1] * np.cos(popt[2] * (np.pi/180))
                                      + np.cos(popt[2] * (np.pi/180) + v[n]))

    # Grafik:
    plt.plot(zp, vr1, 'b-',
             label='RV mit den optimierten Orbitalparametern berechnet')
    plt.legend(fontsize='x-small')
    plt.xlabel('Phase')
    plt.ylabel('RV [km/s]')
    plt.title(table + ', P=' + str(P) + ' d' +
              '\nK=' + str(popt[0].round(1)) + '+-' +
              str(stdabw[0].round(1)) + ' km/s'
              ', e=' + str(popt[1].round(5)) + '+-' + str(stdabw[1].round(5))
              + ', omega=' + str(popt[2].round(1)) + '+-' +
              str(stdabw[2].round(1)) + '°'
              + ', gamma=' + str(popt[3].round(1)) + '+-' +
              str(stdabw[3].round(1)) + ' km/s',
              fontsize='small')
    plt.pause(1)
    grsave = input(
        'Möchten Sie die Grafik speichern ? Dann antworten Sie mit y :')
    if grsave == 'y':
        plt.savefig(table+'_Phasenplot.pdf')
    plt.close('all')


z = np.linspace(0., 1., 1001)
pipe2(z, popt)
