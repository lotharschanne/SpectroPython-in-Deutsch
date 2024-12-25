#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Abfrage von Linienwellenlängen bei NIST

Es werden abgefragt:
    + Das Element oder das Ion
    + Der Wellenlängenbereich.
Ausgedruckt werden die gefundenen Linien.
Und diejenigen, für die auch beobachtete Wellenlängen und relative Intensitäten
tabelliert sind, werden auch als Excel-Datei abgespeichert.

Created on Mon Sep 27 2023

@author: lothar schanne
"""

from astroquery.nist import Nist
import astropy.units as u


element = input('Linien von welchem Element/Ion möchten Sie abfragen? Bitte in der\
 Form H I oder Fe II oder C IV (also mit Ionisationsgrad) oder H oder Fe (ohne\
 Ionistaionsgrad) eingeben: ')
print('\nIn welchem Wellenlängenbereich sollen die Linien liegen?')
begin = float(input('Untere Wellenlänge in Angström: '))
ende = float(input('Obere Wellenlänge in Angström: '))

table = Nist.query(begin * u.AA, ende * u.AA,  linename=element,
                   wavelength_type='vac+air')
# Die Wellenlängen sind in der Luft gemessen


print('\n', table, '\n')


wellenlängen = table['Observed', 'Rel.']
wellenlängen_beobachtet = wellenlängen[wellenlängen['Observed'].data > 0]
print('\nGefundene (beobachtete)  Linienwellenlängen [Angström] \
mit relativen Intensitäten. Diese Ergebnisse sind als Excel-Datei abgespeichert:')
print(wellenlängen_beobachtet)

# Abspeichern in einem Execel-Sheet;
pdframe = wellenlängen_beobachtet.to_pandas()
pdframe.to_excel(element+'.xlsx')
