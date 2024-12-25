#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Skript erzeugt ein künstliches, per gaußbroadening auf eine bestimmte FWHM,
ausgedrückt durch ein R, verbreitertes Linien-Spektrum

Created on Fri Feb 18 15:38:26 2022

@author: lothar
"""

from PyAstronomy import pyasl
import matplotlib.pylab as plt
import numpy as np
from astropy.io import ascii, fits

# Set up an input spectrum
x = np.linspace(5000.0, 5100.0, 10001)
y = np.ones(x.size)


# Introduce some delta-peaked lines
y[1650] = 0.7
y[1870] = 0.3
y[5050] = 0.1
y[6100] = 0.2
y[6150] = 0.7
y[9000] = 0.4

print("Linien bei: ", x[1650], x[1870],
      x[5050], "\n", x[6100], x[6150], x[9000])
X = [x[1650], x[1870], x[5050], x[6100], x[6150], x[9000]]

R = float(input('Geben Sie das gewünschte Auflösungsvermögen R ein: '))
# Apply Gaussian instrumental broadening, setting the resolution to R.
r, fwhm = pyasl.instrBroadGaussFast(
    x, y, R, edgeHandling='firstlast', fullout=True)

# Add some independent, Gaussian noise
gstd = 0.001
r += np.random.normal(0., gstd, len(x))

print("FWHM used for the Gaussian kernel: ", fwhm, " A")

# Plot the output
plt.plot(x, r, "r-", label="Broadened curve (full)")
plt.plot(x, y, "b-", label="Input")
plt.legend(loc=4)
plt.show()


f = open("Kunstspektrum_Linienliste.txt", "w")
for i in X:
    f.write(str(i) + "\n")
f.write("FWHM: " + str(fwhm))
f.close()

header = fits.Header()
header["SIMPLE"] = "T"
header["BITPIX"] = -32
header["NAXIS"] = 1
header["CUNIT1"] = "Angstrom"
header["CTYPE1"] = "Wavelength"
header["CRVAL1"] = x[0]
header["NAXIS1"] = len(x)
header["CRPIX1"] = 1.0
header["CDELT1"] = (x[1] - x[0])
fits.writeto('R'+str(int(R))+'_Kunstspektrum_FWHM' + str(fwhm) + 'A.fits', r, header,
             overwrite=True, output_verify="silentfix")
fits.writeto("Kunstspektrum_Original_ohneVerbreiterung.fits", y, header,
             overwrite=True, output_verify="silentfix")
