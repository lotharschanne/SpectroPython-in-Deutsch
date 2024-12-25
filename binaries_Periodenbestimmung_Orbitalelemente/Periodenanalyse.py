# -*- coding: utf-8 -*-
"""
Liest eine Datei (im csv- oder tab-Format) mit Spaltenüberschriften ein. Die
Spaltenüberschriften für die unabhängige Variable (z.B. 'JD') und die abhängige 
Variable (z.B. RV [km/s]) sind anzugeben.
Macht dann mit beiden Spalten eine Periodenanalyse (Lomb-Scargle-Methode und
Stringlängenmethode),
zeigt grafisch die Wertepaare, das Periodogramm und den Phasenplot
und gibt die Periode (in der Einheit der unabhängigen Variablen) aus.


Lothar
20231229
"""


from PyAstronomy import pyTiming as pyt
import numpy as np
from astropy.timeseries import LombScargle, TimeSeries
from astropy.io import ascii
from astropy import units as u
from astropy import time
import matplotlib.pyplot as plt

plt.switch_backend('Qt5Agg')  # Interaktives GrafikBackend einschalten

# Mittlerer Fehler der RV-Messungen, bitte anpassen
fehler = 3.  # km/s

filename = input("Geben Sie den Name des Datenfiles ein: ")
x = input('Geben Sie den Namen der unabhängigen Variablen ein, z.B. "JD": ')
y = input('Geben Sie den Namen der abhängigen Variablen ein, z.B. "RV": ')

ts = ascii.read(filename)
ts_JD = time.Time(ts.columns[x], format='jd')
ts_RV = ts.columns[y]

print('Geben Sie nachfolgend den Beginn und das Ende des Periodenbereichs ein,\
der überprüft werden soll.')
pa = float(input('Geben Sie das minimale P in der Einheit Tage ein: '))
pe = float(input('Geben Sie das maximale P in der Einheit Tage ein: '))


# Time Series Objekt erstellen
serie = TimeSeries(time=ts_JD)

# RV's einfügen
serie.add_column(ts_RV)

fig = plt.figure()
plt.plot(serie.time.value, serie[y], "o")
# plt.savefig('Werteplot', format='pdf')

# ************** Lomb-Scargle aus astropy.timeseries **************
frequency, power = LombScargle(serie.time.value, serie[y],
                               fehler).autopower(
    minimum_frequency=1/pe, maximum_frequency=1/pa)

fig = plt.figure()
plt.plot(frequency, power)
plt.title("Lomb-Scargle")
# plt.savefig('Lomb_Scargle', format='pdf')


power_max = power.max()
frequency_max = frequency[power.argmax()]
Periode = 1 / frequency_max * u.d
print("LombScargle-Methode aus astropy ergibt Periode = ", Periode.round(6))


# Folding (Phase)
serie_folded = serie.fold(period=Periode, epoch_time=ts_JD[0])
fig = plt.figure()
plt.plot(serie_folded.time.jd, serie_folded[y], 'ko', markersize=3)
plt.xlabel('Time (days)')
plt.ylabel('RV [km/s]')
plt.title(filename + '\nLomb-Scargle-Phasendiagramm mit der Periode ' +
          str(Periode.round(5))+' berechnet')
# plt.text(0, 0, 'Periode '+str(Periode.round(5)))
plt.savefig(filename + '_Phasenplot_LombScargle', format='pdf')

# plt.show(block=True)
plt.pause(.1)


# ****************** weitere Methoden zur Periodenbestimmung *************
jd = ts.columns[x]
rv = ts.columns[y]

# ************** Stringlängenmethode *****************
# Periodenbereich und Schrittezahl bitte anpassen
tps = (pa, pe, 10000)  # getesteter Periodenbereich [min, max, Schritteanzahl]
# Calculate string length
p, sl = pyt.stringlength_dat(jd, rv, tps)
per2 = p[np.argmin(sl)] * u.d
print('Periode nach der Stringlängenmethode = ', per2.round(6))
# Show the string length.
fig = plt.figure()
plt.plot(p, sl, 'b.-')
plt.ylabel("String length")
plt.xlabel("Trial period")
plt.title('Stringlängenmethode')
plt.pause(.1)

# Folding (Phase)
serie_folded2 = serie.fold(period=per2, epoch_time=ts_JD[0])
fig = plt.figure()
plt.plot(serie_folded2.time.jd, serie_folded2[y], 'ko', markersize=3)
plt.xlabel('Time (days)')
plt.ylabel('RV [km/s]')
plt.title(filename + '\nStringlängen-Phasendiagramm mit der Periode ' +
          str(per2.round(5))+' berechnet')
# plt.text(0, 0, 'Periode '+str(per2.round(5)))
plt.pause(.1)
plt.savefig(filename + '_Phasenplot_Stringlaenge', format='pdf')

print('Zum beenden des Programms in das zuletzt geöffnete Diagramm (Figure 5) klicken')
plt.waitforbuttonpress(-1)
plt.close('all')
