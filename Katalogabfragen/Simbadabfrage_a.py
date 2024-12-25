#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Skript zu Abfragen aus Simbad

noch unvollständig. Experimentell.

Created on Sat Sep  9 16:09:43 2023

@author: lothar
"""

from astroquery.simbad import Simbad


customSimbad = Simbad()  # eigene Simbad-Instanz erzeugen
customSimbad.TIMEOUT = 120
customSimbad.ROW_LIMIT = -1  # alle Daten anzeigen

obj = input('Daten von welchem Objekt möchten Sie auswählen. Objektname eingeben: ')

object_Ids_table = customSimbad.query_objectids(obj)
print('n\Bezeichnerliste für das Objekt')
object_Ids_table.pprint()

result_table = customSimbad.query_object(obj)
result_table.pprint()

# # Ausdruck der in customSimbad enthaltenen Felder
# print('Ausdruck der in customSimbad enthaltenen Felder:\n',
#       customSimbad.get_votable_fields())

# # Welche Felder/Werte man in die Abfrage übernehmen kann, kann man mit folgendem
# # Befehl erfahren.
# # customSimbad.list_votable_fields()
# # Die Beschreibung eines Feldes erhält man mit
# # print('\n Beschreibung des Feldes:')
# print('\nBeschreibung des Parameters:')
# Simbad.get_field_description('main_id')

# # Einbeziehen der gewünschten Datentypen in die Abfrage
# customSimbad.add_votable_fields('fe_h', 'mk')

# # Ausdruck der jetzt in customSimbad enthaltenen Felder
# print('Ausdruck der in customSimbad enthaltenen Felder:\n',
#       customSimbad.get_votable_fields())

# # Simbad-Abfrage mit Bedingungen:
# result = customSimbad.query_criteria('id = bet Cyg')
# print('\nSpaltennamen in result:\n', result.colnames, '\n')

# # Ausdrucken der im nachfolgenden Befehl aufgeführten interessierenden Datenfelder
# print(result['MAIN_ID', 'Fe_H_Teff', 'Fe_H_log_g'])

# Erzeugen einer Pandas-Tabelle aus den interessierenden Werten
# pandastabelle = result['MAIN_ID', 'MK_Spectral_type', 'Fe_H_Teff', 'Fe_H_log_g'].to_pandas()
