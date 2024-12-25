#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The script calculates a theoretical phase-dependent sum spectrum of the 2
or 3 components of a 1-orbit or 2-orbit star system. This is subtracted from the
measured spectra, so that a difference spectrum is obtained which contains the effects
effects (e.g. emissions) that are not included in the theoretical spectra of the components.
The script asks for the file in which the fluxes of the two (or 3) components
that originate from the Phoebe II program.
Then a function 'anteile' is defined, which is calculated from the passed phase
of the inner (and outer) orbit and returns the proportions (weights) of the
spectra of the 2 (3) stars and returns them.
The program then asks for the measured fits spectra for which the difference spectra
are to be calculated. They must all have the same step size and the same
wavelength range and have an header entry 'JD'. From them
the observation date JD is read
The JD is used to calculate the radial velocities of the 2 (3) stars
from the orbital elements and combined with the system velocity to form a total RV.
by which the respective theoretical spectra of the 2 (3) stars are shifted.
The spectra are then used to calculate the theoretical sum spectrum
The difference between the measured flux and the theoretical flux then forms the
sought difference spectrum.
In addition to the difference spectrum, all spectra can also be saved by
commenting out the lines of code.


Created on Sat Dec  2 09:09:28 2023

@author: lothar
"""


from scipy.interpolate import interp1d
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii, fits
import glob
from PyAstronomy import pyasl

plt.switch_backend('Qt5Agg')  # Switch on interactive graphics backend


def newton(E):
    return E - (E - args1['e'] * np.sin(E) - args1['M']) / (1.0 - args2['e'] * np.cos(E))


#########################
# Calculation of the weights of the theoretical spectra of Algol A, B and C
# to form the sum spectrum for the respective phase of the inner orbit
filename = input(
    "Enter the name of the data file with the Phoebe fluxes of the spectra: ")
ts = ascii.read(filename, format='csv')
print(ts.colnames)

# Adjust !!!!!!!!!!!!!!!!!!!!!!!!!!!
phase = ts.columns['col4']
anteil_A = ts.columns['col6']
anteil_B = ts.columns['col8']
# anteil_C = ts.columns['col10']
summe = ts.columns['col12']
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# adjust A, B, C  !!!!!!!!!
A = input('Enter the filename of the theoretical spectrum of star A: ')
B = input('Enter the filename of the theoretical spectrum of star B: ')
# C = input('Enter the filename of the theoretical spectrum of star C: ')

flux_A, header = fits.getdata(A, header=True)
integ_A = flux_A.mean()
print('Mean flux star A: ', str(integ_A))

flux_B, header = fits.getdata(B, header=True)
integ_B = flux_B.mean()
print('Mean flux star B', str(integ_B))

# ADJUST !!!!!!!
# flux_C, header = fits.getdata(C, header=True)
# integ_C = flux_C.mean()
# print('Mittlerer Flux Stern C: ', str(integ_C))


def anteile(PhaseInnererOrbit):
    ant_A = interp1d(phase, anteil_A)(PhaseInnererOrbit) / integ_A
    ant_B = interp1d(phase, anteil_B)(PhaseInnererOrbit) / integ_B
    # ant_C = interp1d(phase, anteil_C)(PhaseInnererOrbit) / integ_C  # Anpassen !!!!!!!
    # Sum = ant_A + ant_B + ant_C
    # return ant_A, ant_B, ant_C, Sum
    Sum = ant_A + ant_B
    return ant_A, ant_B, Sum


########################
# Read in the bc-corrected spectra, all with the same step size and
# same wavelength range, extraction of the JD (observation times)

files = input("Path and name of the fits-spectra (use wildcards) : ")
filelist = glob.glob(files)

# Sort alphabetically. If the spectrum files are named correctly, this results # in a temporal order.
filelist.sort()

# Printout of the list for control purposes.
print("\nList of spectra:")
print(filelist)
print("\nNumber of spectra: ", len(filelist), "\n")

t = np.zeros(len(filelist))


# Processing of the individual barycentrically corrected spectra
# Calculation of the JD from the spectrum headers ( = t)
for i in range(len(filelist)):
    flux, header = fits.getdata(filelist[i], header=True)
    t[i] = float(header["JD"])
    if 'CRPIX1' not in header:
        header["CRPIX1"] = 1

    wave = np.ones(header["NAXIS1"], dtype=float)

    for j in np.arange(header["NAXIS1"]):
        wave[j] = (
            header["CRVAL1"]
            + (j - header["CRPIX1"] + 1) * header["CDELT1"]
        )
        # The list wave contains the wavelengths of the pixels.

    # Star A   ###########################
    # # Orbital parameters of the binary system, omega (w) in degrees),
    # P and T in JD, gamma (system velocity) and K1 for Algol A in km/s,
    # adjust !!!!!!!!!!!!!!!!!!!!
    OP = {
        'gamma': 0,
        'K1': 44.1,
        'T': 2441773.49,
        'e': 0.,
        'w1': 90,
        'P': 2.867328}

    # Conversion of the time relative to To and P to determine the phase zp
    zp = (t[i] - OP['T']) / OP['P'] - (t[i] - OP['T']) // OP['P']

    # Calculation of the mean anomaly at the phase zp:
    M = zp * 2 * np.pi  # M in rad

    # Determine the eccentric anomaly E by solving using Newton's method:
    m = 0
    E = 1
    while m >= 0:  # Iterative calculation of the eccentric anomaly E
        args1 = {'e': OP['e'], 'M': M}
        args2 = {'e': OP['e'], 'M': M}
        alt = E
        E = newton(alt)
        m += 1
        if abs(alt - E) < 10e-10:
            break

    # Calculation of the true anomaly v:
    v = 2 * np.arctan(np.tan(E / 2) * np.sqrt(1 + OP['e'])
                      / np.sqrt(1 - OP['e']))
    # Calculation of the radial velocity of star A, inner orbit:
    vr_A_innen = OP['gamma'] + OP['K1'] *\
        (OP['e']*np.cos(OP['w1']*(np.pi/180)) + np.cos(OP['w1']*(np.pi/180) + v))

    # # Orbital parameters of the outer binary system, omega (w) in degrees),
    # # P and T in JD, gamma (system velocity) and K1 of the center of gravity of the
    # # inner orbit in km/s,
    # OP = {
    #     'gamma': 0,
    #     'K1': 11.9,
    #     'T': 2454433.2,
    #     'e': 0.227,
    #     'w1': 138.1,
    #     'P': 680.168}

    # # Umrechnung des Zeitpunktes relativ zu To und P zur Ermittlung der Phase
    # zp = (t[i] - OP['T']) / OP['P'] - (t[i] - OP['T']) // OP['P']

    # # Berechnung der Mittleren Anomalie an der Phase zp:
    # M = zp * 2 * np.pi  # M in rad

    # # Ermittlung der exzentrischen Anomalie E durch Lösung per Newtonverfahren:
    # m = 0
    # E = 1
    # while m >= 0:  # Iterative Berechnung der exzentrischen Anomalie E
    #     args1 = {'e': OP['e'], 'M': M}
    #     args2 = {'e': OP['e'], 'M': M}
    #     alt = E
    #     E = newton(alt)
    #     m += 1
    #     if abs(alt - E) < 10e-10:
    #         break

    # # Berechnung der wahren Anomalie v:
    # v = 2 * np.arctan(np.tan(E / 2) * np.sqrt(1 + OP['e'])
    #                   / np.sqrt(1 - OP['e']))
    # # Berechnung der Radialgeschwindigkeit von Algol A, äusserer Orbit:
    # vr_A_aussen = OP['gamma'] + OP['K1'] *\
    #     (OP['e']*np.cos(OP['w1']*(np.pi/180)) + np.cos(OP['w1']*(np.pi/180) + v))

    # Total speed star A incl. system speed:
    vr_A = vr_A_innen + OP['gamma']  # + vr_A_aussen


######################
    # Star B
    # # Orbitalparameter des inneren Binärsystems
    # omega (w) in Grad), P und T in JD,
    # gamma (Systemgeschwindigkeit) und K1 von AlgolB in km/s,
    # anpassen !!!!!!!!!!!!!!!!
    OP = {
        'gamma': 0,
        'K1': 194.2,
        'T': 2441773.49,
        'e': 0.,
        'w1': 90+180,
        'P': 2.867328}

    # Umrechnung des Zeitpunktes relativ zu To und P zur Ermittlung der Phase
    zp = (t[i] - OP['T']) / OP['P'] - (t[i] - OP['T']) // OP['P']

    # Berechnung der Mittleren Anomalie an der Phase zp:
    M = zp * 2 * np.pi  # M in rad

    # Ermittlung der exzentrischen Anomalie E durch Lösung per Newtonverfahren:
    m = 0
    E = 1
    while m >= 0:  # Iterative Berechnung der exzentrischen Anomalie E
        args1 = {'e': OP['e'], 'M': M}
        args2 = {'e': OP['e'], 'M': M}
        alt = E
        E = newton(alt)
        m += 1
        if abs(alt - E) < 10e-10:
            break

    # Berechnung der wahren Anomalie v:
    v = 2 * np.arctan(np.tan(E / 2) * np.sqrt(1 + OP['e'])
                      / np.sqrt(1 - OP['e']))
    # Berechnung der Radialgeschwindigkeit von Algol B, innerer Orbit:
    vr_B_innen = OP['gamma'] + OP['K1'] *\
        (OP['e']*np.cos(OP['w1']*(np.pi/180)) + np.cos(OP['w1']*(np.pi/180) + v))

    # # Algol B, RV äusseres Orbitsystem ist gleich Stern A, nämlich gleich der RV
    # # des Schwerpunkts des inneren Orbits
    # vr_B_aussen = vr_A_aussen

    # Total speed Algol B incl. system speed:
    vr_B = vr_B_innen + OP['gamma']  # + vr_B_aussen


###########################
    # # Star C
    # # Orbitalparameter des äusseren Binärsystems
    # # omega (w) in Grad),
    # # P und T in JD, gamma (Systemgeschwindigkeit) und K1 von Stern C in km/s,
    # OP = {
    #     'gamma': 0,
    #     'K1': 32.9,
    #     'T': 2454433.2,
    #     'e': 0.227,
    #     'w1': 138.1+180,
    #     'P': 680.168}  # Kolbas, 2015, Table 1

    # # Umrechnung des Zeitpunktes relativ zu To und P zur Ermittlung der Phase
    # zp = (t[i] - OP['T']) / OP['P'] - (t[i] - OP['T']) // OP['P']

    # # Berechnung der Mittleren Anomalie an der Phase zp:
    # M = zp * 2 * np.pi  # M in rad

    # # Ermittlung der exzentrischen Anomalie E durch Lösung per Newtonverfahren:
    # m = 0
    # E = 1
    # while m >= 0:  # Iterative Berechnung der exzentrischen Anomalie E
    #     args1 = {'e': OP['e'], 'M': M}
    #     args2 = {'e': OP['e'], 'M': M}
    #     alt = E
    #     E = newton(alt)
    #     m += 1
    #     if abs(alt - E) < 10e-10:
    #         break

    # # Berechnung der wahren Anomalie v:
    # v = 2 * np.arctan(np.tan(E / 2) * np.sqrt(1 + OP['e'])
    #                   / np.sqrt(1 - OP['e']))
    # # Berechnung der Radialgeschwindigkeit von Stern C, äusserer Orbit:
    # vr_C_aussen = OP['gamma'] + OP['K1'] *\
    #     (OP['e']*np.cos(OP['w1']*(np.pi/180)) + np.cos(OP['w1']*(np.pi/180) + v))

    # # Gesamtgeschwindigkeit Stern C inkl. Systemgeschwindigkeit:
    # vr_C = vr_C_aussen + 3.7


#########################################
# Shift of the theoretical spectra of Algol A, B and C by the total RVs

# Star A
    flux_A, header = fits.getdata(A, header=True)

    if "NAXIS1" in header:
        nax = header["NAXIS1"]
    else:
        print("NAXIS1 fehlt im header !")
    if "CRVAL1" in header:
        crval = header["CRVAL1"]
    else:
        print("CRVAL1 fehlt im header !")
    if "CRPIX1" in header:
        crpix = header["CRPIX1"]
    else:
        crpix = 1
    if "CDELT1" in header:
        cdel = header["CDELT1"]
    else:
        print("CDELT1 fehlt im header !")

    #   Creation of numpy arrays with the wavelengths and fluxes of the spectrum
    wave_A = np.ones(nax, dtype=float)
    for k in range(nax):
        wave_A[k] = crval + (k - crpix + 1) * cdel

    # Shift that spectrum
    flux_rv_A, wave_rv_A = pyasl.dopplerShift(
        wave_A, flux_A, vr_A, edgeHandling="firstlast")

    # # Writing the RV-corrected spectrum to dat-file
    # ascii.write(
    #     [wave_A, flux_rv_A],
    #     'AlgolA_' + "_OrbitRVcorrected" + filelist[i].rsplit('.')[0] + '.dat',
    #     overwrite=True,
    #     names=["WAVE", "FLUX"],
    #     format="tab",
    # )

    # # Writing the RV-corrected spectrum to fits-file
    # header["CRVAL1"] = wave_A[0]
    # header["CRPIX1"] = 1
    # header["NAXIS1"] = len(wave_A)
    # newfile = 'AlgolA' + "_OrbitRVcorrected" + filelist[i]
    # header["RV_COR"] = (vr_A, "km/s, corrected")
    # fits.writeto(newfile, flux_rv_A, header, overwrite=True,
    #              output_verify="silentfix")


# Star B
    flux_B, header = fits.getdata(B, header=True)

    if "NAXIS1" in header:
        nax = header["NAXIS1"]
    else:
        print("NAXIS1 fehlt im header !")
    if "CRVAL1" in header:
        crval = header["CRVAL1"]
    else:
        print("CRVAL1 fehlt im header !")
    if "CRPIX1" in header:
        crpix = header["CRPIX1"]
    else:
        print("CRPIX1 fehlt im header !")
        crpix = 1
    if "CDELT1" in header:
        cdel = header["CDELT1"]
    else:
        print("CDELT1 fehlt im header !")

    #   Creation of numpy arrays with the wavelengths and fluxes of the spectrum
    wave_B = np.ones(nax, dtype=float)
    for k in range(nax):
        wave_B[k] = crval + (k - crpix + 1) * cdel

    # Shift that spectrum
    flux_rv_B, wave_rv_B = pyasl.dopplerShift(
        wave_B, flux_B, vr_B, edgeHandling="firstlast")

    # # Writing the RV-corrected spectrum to dat-file
    # ascii.write(
    #     [wave_B, flux_rv_B],
    #     'AlgolB' + "_OrbitRVcorrected" + filelist[i].rsplit('.')[0] + '.dat',
    #     overwrite=True,
    #     names=["WAVE", "FLUX"],
    #     format="tab",
    # )

    # # Writing the RV-corrected spectrum to fits-file
    # header["CRVAL1"] = wave_B[0]
    # header["CRPIX1"] = 1
    # header["NAXIS1"] = len(wave_B)
    # newfile = 'AlgolB' + "_OrbitRVcorrected" + filelist[i]
    # header["RV_COR"] = (vr_B, "km/s, corrected")
    # fits.writeto(newfile, flux_rv_B, header, overwrite=True,
    #              output_verify="silentfix")


# # Star C
#     flux_C, header = fits.getdata(C, header=True)

#     if "NAXIS1" in header:
#         nax = header["NAXIS1"]
#     else:
#         print("NAXIS1 fehlt im header !")
#     if "CRVAL1" in header:
#         crval = header["CRVAL1"]
#     else:
#         print("CRVAL1 fehlt im header !")
#     if "CRPIX1" in header:
#         crpix = header["CRPIX1"]
#     else:
#         print("CRPIX1 fehlt im header !")
#         crpix = 1
#     if "CDELT1" in header:
#         cdel = header["CDELT1"]
#     else:
#         print("CDELT1 fehlt im header !")

#     #  Creation of numpy arrays with the wavelengths and fluxes of the spectrum
#     wave_C = np.ones(nax, dtype=float)
#     for k in range(nax):
#         wave_C[k] = crval + (k - crpix + 1) * cdel

#     # Shift that spectrum
#     flux_rv_C, wave_rv_C = pyasl.dopplerShift(
#         wave_C, flux_C, vr_C, edgeHandling="firstlast")

#     # # SWriting the RV-corrected spectrum to dat-file
#     # ascii.write(
#     #     [wave_C, flux_rv_C],
#     #     'AlgolC' +
#     #     "_OrbitRVcorrected" + filelist[i].rsplit('.')[0] + '.dat',
#     #     overwrite=True,
#     #     names=["WAVE", "FLUX"],
#     #     format="tab",
#     # )

#     # # Schreiben des RV-korrigierten Spektrums in fits-file
#     # header["CRVAL1"] = wave_C[0]
#     # header["CRPIX1"] = 1
#     # header["NAXIS1"] = len(wave_C)
#     # newfile = 'AlgolC' + "_OrbitRVcorrected" + filelist[i]
#     # header["RV_COR"] = (vr_B, "km/s, corrected")
#     # fits.writeto(newfile, flux_rv_C, header, overwrite=True,
#     #              output_verify="silentfix")

    ###########################################
    # Adding the weighted shifted and narrowed spectra to a
    # total spectrum

    ant_A, ant_B, Sum = anteile(zp)
    # For 2 orbits in the 3 star system:
    # ant_A, ant_B, ant_C, Sum = anteile(zp)  # do not forget to ADJUST function anteile

    step = wave[1] - wave[0]
    for w in range(len(wave_A) - 1):
        if abs(wave_A[w] - wave[0]) <= step and (wave_A[w+1] - wave[0]) <= step:
            wave_A = wave_A[w:w+len(wave)]
            flux_rv_A = flux_rv_A[w:w+len(wave)]
            wave_B = wave_B[w:w+len(wave)]
            flux_rv_B = flux_rv_B[w:w+len(wave)]
            # wave_C = wave_C[w:w+len(wave)]  # für 3 Sternesystem
            # flux_rv_C = flux_rv_C[w:w+len(wave)]   # für 3 Sternesystem
            break

    fluxsumme = (ant_A * flux_rv_A + ant_B * flux_rv_B) / Sum
    # fluxsumme = (ant_A * flux_rv_A + ant_B *
    #              flux_rv_B + ant_C * flux_rv_C) / Sum   # für 3 Sternesystem

    # fluxsumme = fluxsumme[w:w+len(wave)]

    # # Saving the theoretical sum spectrum:
    if i == 0:
        name = input('Enter the wavelength range of the spectra for the \
    file name: ')
    header['JD'] = t[i]
    header["CRVAL1"] = wave_A[0]
    header["CRPIX1"] = 1
    header["NAXIS1"] = len(wave_A)
    newfile = 'Summenspektrum' + name + "_OrbitRVcorrected" + filelist[i]
    # header["RV_COR"] = (vr_B, "km/s, corrected")
    fits.writeto(newfile, fluxsumme, header, overwrite=True,
                 output_verify="silentfix")


    ###########################
    # Forming and saving the difference spectra

    difflux = flux - fluxsumme
    header['CRVAL1'] = wave_A[0]
    header['NAXIS1'] = len(wave_A)

    newfile = 'Differenzspectrum' + \
        name + filelist[i].rstrip('.fit')
    fits.writeto(newfile + '.fit', difflux, header, overwrite=True,
                 output_verify="silentfix")

    ascii.write([wave, difflux],
                newfile + '.csv',
                overwrite=True,
                names=['WAVE', 'FLUX'],
                format="csv")
