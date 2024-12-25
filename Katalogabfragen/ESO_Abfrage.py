#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Skript zu Abfragen aus ESO

noch unvollständig. Experimentell.

Created on Sun Sep 10 13:22:06 2023

@author: lothar
"""

from astroquery.eso import Eso

eso = Eso()

# einloggen
eso.login('lschanne')

eso.ROW_LIMIT = -1   # Return all results
surveys = eso.list_surveys()
print('Alle ESO-Surveys:\n', surveys)

survey = input('\nWelchen survey möchten Sie einsehen? Bitte eingeben: ')
obj = input('Daten von welchem Objekt möchten Sie auswählen. Objektname eingeben: ')

table = eso.query_surveys(surveys=survey, cache=False, target=obj)
print('\nSpalten in der Tabelle:\n', table.colnames)
print('\nErste Zeilen der Tabelle:\n', table[:10])

table_target = table[table['Object'] == obj]
print('\nObjektdaten von'+obj+'\n', table_target)

frage = input(
    '\nDie Objektfiles von ESO herunterladen? Wenn ja dann y eingeben: ')
if frage == 'y':
    data_files = eso.retrieve_data(table_target['ARCFILE'],
                                   destination='~/Downloads')
