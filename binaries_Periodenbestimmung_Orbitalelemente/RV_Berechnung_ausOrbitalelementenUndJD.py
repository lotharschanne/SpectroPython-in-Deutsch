#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Das Skript berechnet aus den als bekannt vorausgesetzten Orbitalparametern und 
einer Reihe von Zeitpunkten (in JD, werden aus einer ascii-Tabelle eingelesen) 
die theoretischen Radialgeschwindigkeiten aus.
Schreibt die Ergebnisse in ein csv-file namens RVs.csv und plottet das Phasen-
diagramm.

20231007
@author: lothar schanne
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii

datenherkunft = input('Möchten Sie die Zeitpunkte aus einer\
ascii-Datei einlesen? Wenn ja dann y eingeben: ')
if datenherkunft == 'y':
    t = input('Geben Sie den Dateinamen der Zeitpunkte ein: ')
    t = np.loadtxt(t)
else:
    # Beispiel-Zeitpunkte der Beobachtungen (JD),
    t = np.array([
        2.459971531251342967e+06,
        2.459972730320548173e+06,
        2.459977545927382074e+06,
        2.459990486079244874e+06,
        2.459994471415292937e+06,
        2.459995491419071797e+06])

# # Orbitalparameter des Binärsystems,
# omegas (w) in Grad),
# P und T in JD, gamma (Systemgescheindigkeit) und K1 und K2 in km/s,
# bitte anpassen an den konkreten Fall #############################
OP = {
    'gamma': 0.82407167,
    'K1': 31.8565906,
    'K2': 29.8379765,
    'T': 2459685.699575991,
    'e': 0.00212,
    'w1': 0.5168984*180/np.pi,
    'w2': (0.5168984 + np.pi)*180/np.pi,
    'P': 71.656017}

# Umrechnung der Zeitpunkte relativ zu To und P (in Phase):
zp = (t - OP['T']) / OP['P'] -\
    (t - OP['T']) // OP['P']

# Mittlere Anomalien an den Zeitpunkten zp:
M = zp * 2 * np.pi

# Ermittlung der exzentrischen Anomalien E durch Lösung per Newtonverfahren


def newton(E):
    return E - (E - args1['e'] * np.sin(E) - args1['M']) / (1.0 - args2['e'] * np.cos(E))


E = np.ones(len(M))
v = np.ones(len(M))
vr1 = np.ones(len(M))
vr2 = np.ones(len(M))


for n in range(len(M)):
    iter = 0
    while n >= 0:  # Iterative Berechnung der exzentrischen Anomalie E
        args1 = {'e': OP['e'], 'M': M[n]}
        args2 = {'e': OP['e'], 'M': M[n]}
        alt = E[n]
        E[n] = newton(alt)
        iter += 1
        # print(n, 'Iter.: ', iter, 'E :', E[n])
        if abs(alt - E[n]) < 10e-10:
            break

    # Berechnung der wahren Anomalien v:
    v[n] = 2 * np.arctan(np.tan(E[n] / 2) * np.sqrt(1 + OP['e'])
                         / np.sqrt(1 - OP['e']))
    # Berechnung der beiden Radialgeschwindigkeiten:
    vr1[n] = OP['gamma'] + OP['K1'] *\
        (OP['e']*np.cos(OP['w1']*(np.pi/180))
         + np.cos(OP['w1']*(np.pi/180) + v[n]))

    vr2[n] = OP['gamma'] + OP['K2']\
        * (OP['e']*np.cos(OP['w2']*(np.pi/180))
           + np.cos(OP['w2']*(np.pi/180) + v[n]))
#    Ergebisausdruck
for i in range(len(vr1)):
    print('JD:', t[i], 'rv1:', vr1[i].round(2), 'rv2:', vr2[i].round(2))

ascii.write([t, vr1, vr2], 'RVs.csv', overwrite=True, names=['JD', 'RV1', 'RV2'],
            format="csv")

# Grafik:
plt.plot(zp, vr1, 'or', label='vr1')
plt.plot(zp, vr2, 'ob', label='vr2')
plt.legend()
plt.pause(1)
grsave = input('Möchten Sie die Grafik speichern ? Dann antworten Sie mit y :')
if grsave == 'y':
    plt.savefig('Phasenplot.pdf')
plt.close('all')
