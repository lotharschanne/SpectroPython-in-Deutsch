#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Das Skript dient der möglichst exakten und reproduzierbaren Normierung von
1d-Spektren im Umfeld einzelner Linien in Serien von Spektren im fits-Format.
In Zeile 66 muß die zu normierende Linie durch ihre Wellenlänge definiert
werden. Außerdem eine Liste von Wellenlängen objektspezifischer Stützpunkte 
(ab Zeile 70), durch die im Verlauf der Berechnungen ein Spline (oder durch
Auskommentieren in Zeile 188 ein Polynom) gelegt wird, der oder das die 
Normierungsfunktion defininert, durch die dann diee Intensitäten (flux) geteilt
wird. Alternativ können die Stützpunkte auch interaktiv grafisch gewählt werden.

Die Grafiken werden gespeichert und die normierten Ordnungen als fits und 
ascii-Datei abgespeichert.

Eine interaktiv angelegte Liste der Stützpunkte wird zudem als ascii-File namens
Stützpunktliste.txt zur Wiederverwendung abgespeichert. Diese kann für weitere 
gleichartige Normierungen per copy/paste in die Stützpunkt-Liste in Zeile 69 ff 
eingetragen werden.

20241126
@author: Lothar Schanne
"""

import numpy as np
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import glob
from scipy.interpolate import make_smoothing_spline


def wavecalc(flux, header):
    """
    Berechnet die Wellenlängen aus den Headerdaten der fits-Datei
    """
    wave = np.zeros(header["NAXIS1"])
    if "CRPIX1" not in header:
        header["CRPIX1"] = 1
    header["CRVAL1"] = header["CRVAL1"] + \
        (1 - header["CRPIX1"]) * header["CDELT1"]
    for i in np.arange(header["NAXIS1"]):
        wave[i] = header["CRVAL1"] + i * header["CDELT1"]

    return wave, flux


plt.switch_backend('Qt5Agg')  # Interaktives GrafikBackend einschalten

# Create file list for objects. All spectra in one (sub)folder.
files = input('Path and name of the OBJECT files (use wildcards) : ')
filelist = glob.glob(files)
filelist.sort()

# Print the list
# print('\nList of spectra: \n')
print('Number of spectra: ', len(filelist), '\n')

########################################################################

# Hier anschließend die gewünschten Linie und die Stützpunktwellenlängen
# eintragen. Die Liste der Stützpunkte und obige Konstante namens ueberhang
# definieren den Ausschnitt des Spektrums, der normiert wird. Die Liste muß
# mindestens 5 Einträge enthalten

# H alpha:
wellenlaenge = 6562.817   # H alpha Laborwellenlänge, hier die gewünschte
# Wellenlänge eintragen !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# Liste der Stützpunkte
# stuetzpunkte = np.asarray(
#     [6498, 6511, 6525, 6588, 6596, 6608, 6616])  # für sehr breite Ha-Linien
stuetzpunkte = np.asarray([6500, 6506, 6510, 6520.7, 6529, 6539, 6591, 6604,
                           6608, 6615, 6620])  # für schmale Ha_Linien
# stuetzpunkte = np.asarray(
#     [6498, 6511, 6525, 6551, 6576, 6588, 6596, 6608, 6616])  # Ha für bet Cep
# stuetzpunkte = np.asarray(
#     [6498, 6511, 6525, 6541, 6588, 6596, 6608, 6616])  # Ha für bet Lyr


if input('Möchten Sie die Stützpunkte interaktiv grafisch per Maus im ersten\
 Spektrum festlegen? Dann y eingeben: ') == 'y':
    flux, header = fits.getdata(filelist[0], header=True)
    wave, flux = wavecalc(flux, header)

    plt.figure(figsize=(15, 8))
    plt.plot(wave, flux, linewidth=.5, label='Spektrum')
    plt.xlabel("Wavelength [Angström]")
    plt.ylabel("ADU")
    plt.title("Spektrum " + filelist[0])
    plt.grid(True)
    print('Bitte  das Grafikfenster bearbeiten (die Linie heraus \
vergrössern) und dann auf nächste Anweisung warten')
    plt.pause(10)
    print('Jetzt mit Mausklicks mindestens 6 Stützpunkte von links nach rechts \
auswählen und mit der Taste Esc beenden.')

    pts = []
    pts = np.asarray(plt.ginput(n=-1, timeout=-1, show_clicks=True))
    plt.plot(pts[:, 0], pts[:, 1], "o", markersize=3)

    stuetzpunkte = pts[:, 0]

    fileobj = open('Stützpunkteliste.txt', 'w')
    fileobj.writelines(str(list(stuetzpunkte)))
    fileobj.close()

    plt.close()
else:
    # Hier die gewünschten Linie und die Stützpunktwellenlängen
    # eintragen. Die Liste der Stützpunkte und obige Konstante namens ueberhang
    # definieren den Ausschnitt des Spektrums, der normiert wird. Die Liste muß
    # mindestens 6 Einträge enthalten

    # H alpha:
    wellenlaenge = 6562.817   # H alpha Laborwellenlänge, hier die gewünschte
    # Wellenlänge eintragen !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    # Liste der Stützpunkte
    # stuetzpunkte = np.asarray(
    #     [6498, 6511, 6525, 6588, 6596, 6608, 6616])  # für sehr breite Ha-Linien
    stuetzpunkte = np.asarray([6500, 6506, 6510, 6520.7, 6529, 6539, 6591, 6604,
                               6608, 6615, 6620])  # für schmale Ha_Linien
    # stuetzpunkte = np.asarray(
    #     [6498, 6511, 6525, 6551, 6576, 6588, 6596, 6608, 6616])  # für bet Cep
    # stuetzpunkte = np.asarray(
    #     [6498, 6511, 6525, 6541, 6588, 6596, 6608, 6616])  # für bet Lyr


############################################################################

# Normierungsprogramm:

# Festlegung der Rigizität des Splines, hohe Werte entsprechen großer Steifheit,
# geringere einer geringen Steifheit. Sollte angepaßt werden.
smoothing = 1000

# Halbe Breite (in Pixel) des Intervalls, in dem um einen Stützpunkt der Flux
# gemittelt wird. Evtl. anpassen.
intervall = 0

ueberhang = 0  # Pixel
polynomorder = 5  # Ordnung des gefitteten Polynoms

# Abarbeiten der filelist
for i in range(len(filelist)):
    print()
    print(filelist[i], ':')

    flux, header = fits.getdata(filelist[i], header=True)
    wave, flux = wavecalc(flux, header)

    step = float(header["CDELT1"])
    lambda0 = wave[0]

    # index_wellenlaenge = int((wellenlaenge - lambda0) / step)
    index_start = int((stuetzpunkte[0] - lambda0) / step) - ueberhang
    index_stop = int((stuetzpunkte[-1] - lambda0) / step) + ueberhang
    indices_stuetzpunkte = (stuetzpunkte - lambda0) / step
    indices_stuetzpunkte = np.round(indices_stuetzpunkte).astype(int)

    croppedflux = flux[index_start:index_stop]
    croppedwave = wave[index_start:index_stop]

    flux_stuetzpunkte = np.zeros(len(indices_stuetzpunkte))

    for k in range(len(indices_stuetzpunkte)):
        flux_stuetzpunkte[k] = np.mean(
            flux[indices_stuetzpunkte[k]-intervall:indices_stuetzpunkte[k]+intervall+1])

    plt.figure(figsize=(15, 8))
    plt.plot(wave, flux, linewidth=.3, label='Spektrum')
    plt.plot(croppedwave, croppedflux, linewidth=.5, label='Linienausschnitt')
    plt.scatter(wave[indices_stuetzpunkte],
                flux_stuetzpunkte, c='k', s=15, label='Stützpunkte')
    plt.title(filelist[i] + ', Linienausschnitt und Stützpunkte')
    plt.grid(True)
    plt.legend()
    plt.pause(.1)
    plt.savefig(filelist[i]+'_Ausschnitt.png')
    plt.close()

    # Normierungsfunktion als Spline berechnen
    normfunc = make_smoothing_spline(
        wave[indices_stuetzpunkte],
        flux_stuetzpunkte, lam=smoothing)
    normfunction = normfunc(croppedwave)

    # Normierungsfunktion als Polynom berechnen
    # polynomcoeff = np.polyfit(wave[indices_stuetzpunkte],
    #                           flux_stuetzpunkte, polynomorder)
    # polynomfunc = np.poly1d(polynomcoeff)
    # normfunction = polynomfunc(croppedwave)

    plt.figure(figsize=(15, 8))
    plt.plot(croppedwave, croppedflux, 'g-', linewidth=.3, label='Spektrum')
    plt.plot(croppedwave, normfunction, 'b-', label='Normierungsfunktion')
    plt.scatter(wave[indices_stuetzpunkte],
                flux_stuetzpunkte, c='k', s=15, label='Stützpunkte')
    plt.title(filelist[i])
    plt.grid(True)
    plt.legend()
    plt.pause(.1)
    plt.savefig(filelist[i] + '_Normierungsfunktion.png')
    plt.savefig(filelist[i] +
                '_Normierungsfunktion.pdf', format='pdf')
    plt.close()

    flux_normiert = croppedflux / normfunction

    plt.figure(figsize=(15, 8))
    plt.plot(croppedwave, flux_normiert,  'b-',
             label=filelist[i], linewidth=.6)
    plt.grid(True)
    plt.title(filelist[i] + ' normiert', fontsize=12)  # fontsize anpassen !!!
    plt.xlabel('Wavelength in Angström')
    plt.ylabel('normierter Flux')
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.savefig(filelist[i] + '_normiert.png', format='png')
    plt.savefig(filelist[i] + '_normiert.pdf', format='pdf')
    plt.pause(.1)
    plt.close()

    header['CRVAL1'] = croppedwave[0]
    header['CRPIX1'] = 1
    fits.writeto(filelist[i].rsplit('.', 1)[0] + '_normiert.fits',
                 flux_normiert, header,
                 overwrite=True, output_verify='silentfix')

    ascii.write([croppedwave, flux_normiert],
                filelist[i].rsplit('.', 1)[0] + "normiert.dat",
                overwrite=True, names=['WAVE', 'FLUX'], format='tab')

# Programmende
plt.close('all')
print("END of program, all done")
