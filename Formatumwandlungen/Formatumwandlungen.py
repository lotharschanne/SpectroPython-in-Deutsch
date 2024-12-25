#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Formatumwandlungen.py

Umwandlung JD -> Kalenderdatum
Umwandlung Kalenderdatum -> JD
Umwandlung hexagesimale RA und DEC in float
Umwandlung RA und DEC in float in hexagesimale RA und DEC

20221208

@author: lothar schanne
"""

from PyAstronomy import pyasl
from astropy.time import Time

print("Auswahl der Aufgabe:")

typ = input(
    """Wenn Sie ein Julianisches Datum in Kalenderdaten umwandeln wollen, geben Sie JD ein: \n\
Wenn Sie ein Kalenderdatum in Julianisches Datum umwandeln wollen geben sie\
 KD ein: \n\
Wenn Sie RA und DEC in hexagesimaler Form in float umwandeln wollen, geben Sie \
RA ein: \n\
Wenn Sie RA und DEC von float in hexagesimale Form umwandeln möchten, geben Sie \
DEC ein: \n\
Eingabe:
"""
)

if typ == "JD":
    # Convert JD to calendar date
    jd = float(input("Eingabe des JD: "))
    t = Time(jd, format='jd')
    print("Datum: ", t.strftime('%H:%M:%S %d %b %Y'))
    print()

if typ == "KD":
    # Convert calendar date to JD
    kd = input(
        "Geben Sie das Kalenderdatum ein in der Form 2017-01-19T18:22:45 (=isot): \n")
    t = Time(kd, format='isot', scale='utc')
    print("Corresponding Julian date: ", t.jd)
    print("Corresponding reduced Julian date: ", t.mjd)
    print()

if typ == "RA":
    # Obtain decimal representation
    # The coordinate string. Valid formats are, e.g., “00 05 08.83239 +67 50 24.0135”
    # or “00:05:08.83239 -67:50:24.0135”.
    # Spaces or colons are allowed as separators for the individual components
    # of the coordinates.
    radec = input(
        "Geben Sie RA und DEC in hexagesimaler Form ein. \n\
        Wie 00 05 08.83239 +67 50 24.01: \n"
    )
    print("Eingegebene Koordinaten : ", radec)
    # Wandele sexagesimale Koordinaten in dezimale um
    ra, dec = pyasl.coordsSexaToDeg(radec)
    print("Coordinates  [deg]: %010.6f  %+09.6f" % (ra, dec))

if typ == "DEC":
    ra = float(input("Geben Sie RA als float ein: "))
    dec = float(input("Geben Sie DEC als float ein: "))
    # Convert dezimal into sexagesimal representation
    print("RA, DEC", ra, dec)
    sexa = pyasl.coordsDegToSexa(ra, dec)
    print("Coordinates  [sexa]: ", sexa)
