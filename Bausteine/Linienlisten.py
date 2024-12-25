#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Dictionarys der Linien

Das Skript erzeugt Linienlaborwellenlängen der Elemente/Ionen in Form von 
Dictionarys.

Das Skript ist gedacht, um in anderen Skripten Linienwellenlängen zur Verfügung
zu stellen. Ist also dort an geeigneter Stelle zu importieren.

Created on Wed May 31 13:10:24 2023

@author: lothar
"""


C = {
    "CIV 7726": 7726.2,
    "CIII 7037": 7037.25,
    "CII 5920": 5919.6,
    "CIV 5812": 5812.140,
    "CIV 5802": 5801.510,
    "CIII 5696": 5696.000,
    "CIII 4650": 4650.16,
    "CIII 4647": 4647.400,
    "CIII 4388": 4388.24,
    "CII 5920": 5919.6,
    "CII 4267": 4267.258,
}

Ca = {
    "CaI 6162": 6162.17,
    "CaI 6122": 6122.22,
    "CaI 6103": 6102.72,
    "CaI 5857": 5857.45,
    "CaI 5042": 5041.62,
    'CaI 4227': 4226.728,

}

Cu = {
    "CuII 6013": 6013.411,
}

FeI = {
    "FeI 6678": 6677.9865,
    "FeI 6634": 6633.7492,
    "FeI 6609": 6609.1098,
    "FeI 6546": 6546.24,
    "FeI 6463": 6462.725,
    "FeI 6417": 6416.9386,
    "FeI 6412": 6411.6592,
    "FeI 6408": 6408.0272,
    "FeI 6400": 6400.0008,
    "FeI 6394": 6393.6009,
    "FeI 6265": 6265.1335,
    "FeI 6256": 6256.3611,
    "FeI 6255": 6254.581,
    "FeI 6253": 6252.555,
    "FeI 6233": 6232.6408,
    "FeI 6231": 6230.7226,
    "FeI 6213": 6213.299,
    "FeI 6200": 6200.3125,
    "FeI 6192": 6191.558,
    "FeI 6180": 6180.2038,
    "FeI 6173": 6173.3352,
    "FeI 6170": 6169.597,
    "FeI 6164": 6163.5441,
    "FeI 6142": 6141.7316,
    "FeI 6137": 6137.286,
    "FeI 6065": 6065.482,
    "FeI 6056": 6056.0043,
    "FeI 6027": 6027.0505,
    "FeI 6024": 6024.0576,
    "FeI 6020": 6020.1688,
    "FeI 5860": 5859.608,
    "FeI 5862": 5862.357,
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
    "FeI 5447": 5446.8743,
    "FeI 5430": 5429.6964,
    "FeI 5415": 5415.1989,
    "FeI 5406": 5405.7749,
    "FeI 5383": 5383.3688,
    "FeI 5367": 5367.466,
    "FeI 5307": 5307.36,
    "FeI 5302": 5302.299,
    "FeI 5233": 5232.94,
    "FeI 5162": 5162.2725,
    "FeI 5097": 5096.9977,
    "FeI 5075": 5074.748,
    "FeI 5065": 5065.0181,
    "FeI 5018": 5018.4354,
    "FeI 5007": 5007.275,
    "FeI 5002": 5001.8633,
    "FeI 4383": 4383.5447,
    "FeI 4046": 4045.8122,
}

FeII = {
    "FeII 6516": 6516.0783,
    "FeII 6456": 6456.3805,
    "FeII 6248": 6247.559,
    "FeII 6170": 6169.816,
    "FeII 6149": 6149.231,
    "FeII 6148": 6147.734,
    "FeII 5412": 5411.970,
    "FeII 5363": 5362.9698,
    "FeII 5154": 5154.242,
    "FeII 4924": 4923.9216,
    "FeII 4629": 4629.3317,
    "FeII 4584": 4583.8292,
    "FeII 4542": 4541.985,
    "FeII 4233": 4233.1627,
}

FeIII = {"FeIII 5920": 5920.0,
         "FeIII 5919": 5918.960,
         }

Hg = {
    "HgII 6142": 6141.773,
    "HgII 4687": 4686.563,
}

HI = {
    "H alpha": 6562.817,
    "H beta": 4861.332,
    "H gamma": 4340.468,
    "H delta": 4101.737,
    'H epsilon': 3970.07,
}

HeI = {
    "HeI 6678": 6678.149,
    "HeI 5876": 5875.989,
    "HeI 5016": 5015.6783,
    "HeI 4922": 4921.929,
    "HeI 4713": 4713.143,
    "HeI 4686": 4685.682,
    "HeI 4542": 4541.59,
    "HeI 4471": 4471.688,
    "HeI 4388": 4387.928,
    "HeI 4026": 4026.191,
}

HeII = {
    "HeII 5411": 5411.524,
    'HeII 4200': 4199.83,
}

Mn = {
    "MnI 4031": 4030.76,
}

Mg = {

    "MgI 5183": 5183.6042,
    "MgI 5172": 5172.6843,
    "MgI 5167": 5167.3216,
    "MgI 4571": 4571.0956,
    "MgI 4481": 4481.13,
}

N = {
    "NIII 4200": 4200.020,
}

NaI = {
    "Na D1": 5895.924,
    "Na D2": 5889.951,
}

Ni = {
    "NiI 6177": 6176.81,
    "NiI 6175": 6175.367,
    "NiI 6108": 6108.12,
}

O = {
    "OIII 5592": 5592.370,
}

Si = {
    "SiII 6371": 6371.36,
    "SiII 6347": 6347.10,
    "SiIV 4089": 4088.863,
}

terrestrische = {
    "O2 6883": 6883.818,
    'H2O 6515': 6517.727,
    'H2O 6519': 6519.452,
    'H2O 6524': 6523.843,
    'H2O 6532': 6532.359,
    'H2O 6534': 6534.000,
    'H2O 6537': 6536.720,
    'H2O 6542': 6542.313,
    'H2O 6544': 6543.907,
    'H2O 6548': 6547.705,
    'H2O 6549': 6548.622,
    'H2O 6553': 6552.629,
    'H2O 6554': 6553.785,
    'H2O 6557': 6557.171,
    'H2O 6558': 6558.149,
    'H2O 6561': 6560.555,
    'H2O 6564': 6564.206,
    'H2O 6569': 6568.806,
    'H2O 6572': 6572.086,
    'H2O 6575': 6574.852,
}

Ti = {"TiII 5491": 5490.7,
      "TiII 5419": 5418.8,
      "TiI 5404": 5404.11,
      "TiII 5381": 5381.03,
      "TiII 5262": 5262.14,
      "TiII 5072": 5072.25,
      }

V = {
    "VI 5062": 5061.79,
}

Kunstspektrum = {
    '5017': 5016.5,
    '5019': 5018.7,
    '5051': 5050.5,
    '5061': 5061.0,
    '5062': 5061.5,
    '5090': 5090.0,
}

"""
Falls neue Elemente/Ionen als Dictionarys eingefügt werden, sollte auch die
folgende Liste (Elementliste) ergänzt werden.
"""
Elementliste = ['C', 'Ca', 'Cu', 'FeI', 'FeII', 'FeIII', 'HI', 'Hg', 'HeI', 'HeII',
                'Mg', 'N', 'NaI', 'Ni', 'O', 'Si', 'terrestrische', 'Ti', 'V',
                'Kunstspektrum']


# Funktion zur Linienauswahl:
def linienauswahl():
    print('Zur Verfügung stehende Elemente/Ionen: ')
    print(Elementliste)
    element = eval(
        input('Von welchem Element/Ion möchten Sie eine Linie auswählen?: '))
    print(element)
    linie = input('Geben Sie den Bezeichner der gewünschten Linie ein: ')
    return linie, element[linie]  # Wellenlänge der gewählten Linie


# Funktion zur Elementauswahl:
def elementauswahl():
    print('\nZur Verfügung stehende Elemente/Ionen: ')
    print(Elementliste)
    element = input(
        'Von welchem Element/Ion möchten Sie die Linien auswählen?: ')
    element = eval(element)
    print(element)
    return element
