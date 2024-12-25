#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ansehen eines wählbaren 1d-Spektrums im fits-Format,
Eingabe der Linien, die im Spektrum markiert werden sollen. Die wählbaren
Linien sind ab Zeile 21 in der Liste Linien aufgeführt. Diese Liste kann
beliebig erweitert werden.
Grafik kann auf Wunsch als pdf abgespeichert werden.

Stand 20230905
author: Lothar Schanne
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

plt.switch_backend('Qt5Agg')  # Interaktives GrafikBackend einschalten

# Liste der wählbaren Linien:
Linien = {
    "O2 6883": 6883.818,
    "CIV 7726": 7726.2,
    "CIII 7037": 7037.25,
    "HeI 6678": 6678.149,
    "FeI 6678": 6677.9865,
    "FeI 6634": 6633.7492,
    "FeI 6609": 6609.1098,
    "H alpha": 6562.817,
    "FeI 6546": 6546.24,
    "FeII 6516": 6516.0783,
    "FeI 6463": 6462.725,
    "FeII 6456": 6456.3805,
    "FeI 6417": 6416.9386,
    "FeI 6412": 6411.6592,
    "FeI 6408": 6408.0272,
    "FeI 6400": 6400.0008,
    "FeI 6394": 6393.6009,
    "SiII 6371": 6371.36,
    "SiII 6347": 6347.10,
    "FeI 6265": 6265.1335,
    "FeI 6256": 6256.3611,
    "FeI 6255": 6254.581,
    "FeI 6253": 6252.555,
    "FeII 6248": 6247.559,
    "FeI 6233": 6232.6408,
    "FeI 6231": 6230.7226,
    "FeI 6213": 6213.299,
    "FeI 6200": 6200.3125,
    "FeI 6192": 6191.558,
    "FeI 6180": 6180.2038,
    "NiI 6177": 6176.81,
    "NiI 6175": 6175.367,
    "FeI 6173": 6173.3352,
    "FeII 6170": 6169.816,
    "FeI 6170": 6169.597,
    "FeI 6164": 6163.5441,
    "CaI 6162": 6162.17,
    "FeII 6149": 6149.231,
    "FeII 6148": 6147.734,
    "HgII 6142": 6141.773,
    "FeI 6142": 6141.7316,
    "FeI 6137": 6137.286,
    "CaI 6122": 6122.22,
    "NiI 6108": 6108.12,
    "CaI 6103": 6102.72,
    "FeI 6065": 6065.482,
    "FeI 6056": 6056.0043,
    "FeI 6027": 6027.0505,
    "FeI 6024": 6024.0576,
    "FeI 6020": 6020.1688,
    "CuII 6013": 6013.411,
    "CII 5920": 5919.6,
    "FeIII 5920": 5920.0,
    "FeIII 5919": 5918.960,
    "Na D1": 5895.924,
    "Na D2": 5889.951,
    "HeI_5876": 5875.989,
    "FeI 5860": 5859.608,
    "FeI 5862": 5862.357,
    "HI 5861": 5861.35,
    "CIV 5812": 5812.140,
    "CIV 5802": 5801.510,
    "CIII 5696": 5696.000,
    "OIII 5593": 5592.37,
    "FeI 5576": 5576.0884,
    "FeI 5573": 5572.842,
    "FeI 5570": 5569.6177,
    "FeI 5567": 5567.3907,
    "FeI 5566": 5565.7036,
    "FeI 5560": 5560.2112,
    "FeI 5558": 5557.9818,
    "FeI 5555": 5554.8947,
    "FeI 5526": 5525.5439,
    "FeI 5498": 5497.5157,
    "TiII 5491": 5490.7,
    "FeI 5447": 5446.8743,
    "FeI 5430": 5429.6964,
    "TiII 5419": 5418.8,
    "FeI 5415": 5415.1989,
    "FeII 5412": 5411.970,
    "HeII 5411": 5411.524,
    "FeI 5406": 5405.7749,
    "TiI 5404": 5404.11,
    "FeI 5383": 5383.3688,
    "TiII 5381": 5381.03,
    "FeI 5367": 5367.466,
    "FeII 5363": 5362.9698,
    "FeI 5307": 5307.36,
    "FeI 5302": 5302.299,
    "TiII 5262": 5262.14,
    "FeI 5233": 5232.94,
    "MgI 5183": 5183.6042,
    "MgI 5172": 5172.6843,
    "MgI 5167": 5167.3216,
    "FeI 5162": 5162.2725,
    "FeII 5154": 5154.242,
    "FeI 5097": 5096.9977,
    "FeI 5075": 5074.748,
    "TiII 5072": 5072.25,
    "FeI 5065": 5065.0181,
    "VI 5062": 5061.79,
    "FeI 5018": 5018.4354,
    "HeI 5016": 5015.6783,
    "FeI 5007": 5007.275,
    "FeI 5002": 5001.8633,
    "HeI 4922": 4921.929,
    "H beta": 4861.332,
    "HeI 4713": 4713.143,
    "HgII 4687": 4686.563,
    "HeI 4686": 4685.682,
    "CIII 4650": 4650.16,
    "CIII 4647": 4647.400,
    "MgI 4571": 4571.0956,
    "FeII 4542": 4541.985,
    "HeI 4542": 4541.59,
    "He 4452": 4541.59,
    "HeI 4471": 4471.688,
    "HeI 4388": 4387.928,
    "CIII 4388": 4388.24,
    "H gamma": 4340.468,
    "NIII 4200": 4200.020,
    "HI 4102": 4101.737,
    "SiIV 4089": 4088.863,
    "HeI 4026": 4026.191,
}
print("Linienauswahl: ", Linien.keys())


# Pfad und Name des Spektrumfiles bitte anpassen
file = input("Pfad und Filebezeichnung eingeben: ")

#   Lesen des Spektrums
sp = fits.open(file)

# Header lesen und in der Konsole ausdrucken
# HD = dict(sp[0].header)
# print("\n\nHeader des Spektrums :\n")
# print(HD)

if 'CRPIX1' not in sp[0].header:
    sp[0].header["CRPIX1"] = 1

# Erzeugen von Arrays mit den Wellenlängen und Fluxes des Spektrums
flux = np.array(sp[0].data)
wave = np.ones(sp[0].header["NAXIS1"], dtype=float)
for i in range(sp[0].header["NAXIS1"]):
    wave[i] = (
        sp[0].header["CRVAL1"]
        + (i + 1 - sp[0].header["CRPIX1"]) * sp[0].header["CDELT1"]
    )
# In der Liste wave sind die Wellenlängen der Pixel enthalten
# In der Liste flux die entsprechenden Intensitäten

#   Schliessen des fits-file
sp.close()

# Plot gesamtes Spektrum
fig = plt.figure(1, figsize=(14, 10))
plt.plot(wave, flux, "b-", linewidth=0.5)
plt.xlabel("Wellenlänge [Angström]", fontsize=14)
plt.ylabel("ADU", fontsize=14)
plt.title("Spektrum " + file, fontsize=18)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.grid(True)
plt.pause(.1)

linie = input('Geben Sie den Namen der zu  markierenden Linie ein: ')
line = Linien[linie]
plt.axvline(line, color='r')
plt.text(line, 1.05, linie)
plt.pause(.1)

frage = input('Sollen weitere Linien markiert werden, dann y eingeben. ')
while frage == 'y':
    linie = input('Geben Sie den Namen der zu  markierenden Linie ein: ')
    line = Linien[linie]
    plt.axvline(line, color='r')
    plt.text(line, 1.05, linie)
    frage = input('Sollen weitere Linien markiert werden, dann y eingeben.')

plt.pause(.1)

frage2 = input('Soll die Grafik gespeichert werden? Dann y eingeben: ')
if frage2 == 'y':
    grafik = file.rsplit('.fit', 1)[0] + '.pdf'
    fig.savefig(grafik)

# Zeigen des Plots, wenn das Skript in einer normalen Python-Konsole
# durchgeführt wird. Für die IPython-Konsole in Spyder ist das nicht nötig.
# plt.show(block=True)

print('Ende des Programms')
