#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Berechnet die Höhen der Objekte einer Objektliste (ascii-Datei (csv) mit den
Objektnamen, je eines in einer Zeile, siehe BeispielObjektliste.txt) über dem
Horizont für einen geplanten Zeitpunkt an einem Beobachtungsort und plottet sie
 über einen Zeitraum von 24 Stunden. Die Grafiken werden als pdf im
 Arbeitsverzeichnis gespeichert.
Eingabe der Liste der Sterne, des Beobachters und des Datums/Zeitpunkts nötig.
Die Liste der Sterne muß im Arbeitsverzeichnis stehen oder in PYTHONPATH.

20221108

@author: lothar
"""


import matplotlib.pyplot as plt
import astropy.units as u
from astropy.coordinates import EarthLocation
from astroplan import Observer
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astroplan.plots import plot_airmass, plot_altitude
from pytz import timezone
from astroplan import FixedTarget
import csv


# Daten des Objekts
# coordinates = SkyCoord('9h23m55.42s', '+54d55m31.5s', frame='icrs')
# Einlesen der Sternkoordinaten über das Internet
objektliste = input('Geben Sie den Namen der Liste der Objektsterne ein: ')


# gewünschtes Beobachtungsdatum
time = input('Geben Sie den gewünschten Beobachtungszeitpunkt ein.\
 In der Form YYYY-MM-DD hh:mm:ss   : ')
time = Time(time)  # gewünschte Beobachtungszeit

# Standort des Observatoriums
ort = input('Geben Sie den Namen des Observatoriums (Beobachter) ein: ')

# ++++++++++++++   Liste der Standorte der Beobachter  ++++++++++++++
# bitte analog ergänzen, falls nötig
# Daten des Oberservatoriums, hier Beispiel in Israel, wo NRES installiert ist
if ort == 'NRES':
    longitude = '34.763333'
    latitude = '30.595833'
    elevation = 1500 * u.m
    location = EarthLocation.from_geodetic(longitude, latitude, elevation)
    obs = Observer(name='NRES Israel',
                   location=location,
                   timezone=timezone('Asia/Jerusalem'),
                   description="")

# Daten VEGA Observatorium in Anthing bei Innsbruck, Österreich
if ort == 'VEGA':
    longitude = '13.008816'
    latitude = '47.886710'
    elevation = 800 * u.m
    location = EarthLocation.from_geodetic(longitude, latitude, elevation)
    obs = Observer(name='Flechas Innsbruck',
                   location=location,
                   timezone=timezone('Europe/Berlin'),
                   description="")

# Koordinaten vom Berthold Stober (Glan-Münchweiler)
if ort == 'Berthold':
    longitude = '7.4775'
    latitude = '49.47527777778'
    elevation = 200 * u.m
    location = EarthLocation.from_geodetic(longitude, latitude, elevation)
    obs = Observer(name='Berthold Stober',
                        location=location,
                        timezone=timezone('Europe/Berlin'),
                        description="")

# Koordinaten von Siegfried Hold
if ort == 'Siegfried':
    longitude = '5.68461111'
    latitude = '47.00161111'
    elevation = 380 * u.m
    location = EarthLocation.from_geodetic(longitude, latitude, elevation)
    obs = Observer(name='Siegfried Hold',
                   location=location,
                   timezone=timezone('Europe/Berlin'),
                   description="")

reader = csv.reader(open(objektliste))

for object in reader:
    object = object[0]
    coordinates = FixedTarget.from_name(object)
    # plot_airmass(coordinates, obs, time)
    fig = plt.figure()
    plot_altitude(coordinates, obs, time, airmass_yaxis=True)
    plt.title(object + ', ' + obs.name)
    fig.savefig(object+'.pdf')
    plt.pause(3)

print('Zum beenden des Programms in das zuletzt geöffnete Diagramm klicken')
plt.waitforbuttonpress(-1)
plt.close('all')
