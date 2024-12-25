#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Übernimmt die Lichtkurvendaten aus einer csv-Datei (welche z.B. von der 
AAVSO-Datenbank importiert wurde). Man gibt den auszuwertenden
Zeitraum ein (JD Anfang und Ende) und das Programm berechnet die Tagesmittelwerte
der mag-Werte. Diese werden in eine csv-Datei geschrieben und grafisch dargestellt.

Created on 20230219
@author: lothar
"""

import csv
import numpy as np
from astropy.io import ascii
import matplotlib.pyplot as plt


reader = csv.reader(open('JD_mag.csv'))

data = []
magmittel = []

for row in reader:
    data.append(row)

JD_anfang = int(input('Geben Sie den Beginn des Zeitraums ein (in JD): '))
JD_ende = int(input('Geben Sie das Enden des Zeitraums ein (in JD): '))

for i in range(JD_anfang, JD_ende):
    sum = 0.
    zaehler = 0
    for j in range(len(data)):
        if i == round(float(data[j][0])):
            sum += float(data[j][1])
            zaehler += 1
    magmittel.append(sum / zaehler)

diff = JD_ende - JD_anfang
MittelwertMag = np.zeros((diff, 2))
for k in range(diff):
    MittelwertMag[k][0] = k + JD_anfang
    MittelwertMag[k][1] = magmittel[k]


# Abspeichern als ascii-Datei
ascii.write(
    MittelwertMag,
    'Tagesmittelwerte.csv',
    overwrite=True,
    names=["JD", "mag"],
    format="csv",
)

plt.plot(MittelwertMag[:, 0], -MittelwertMag[:, 1])
plt.xlabel('JD')
plt.ylabel('Helligkeit [mag]')
plt.title('Lichtkurve Nova Cas 2021 nach AAVSO')
plt.grid(True)
plt.savefig('Lichtkurve.png')
plt.savefig('Lichtkurve.pdf')
