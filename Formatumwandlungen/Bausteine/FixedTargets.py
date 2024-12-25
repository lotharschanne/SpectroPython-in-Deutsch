#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 17:27:39 2022

@author: lothar
"""

from astroplan import FixedTarget
from astropy.coordinates import SkyCoord
import astropy.units as u


# sirius = FixedTarget.from_name("Sirius")
# _7and = FixedTarget.from_name("7 and")
# delCep = FixedTarget.from_name("del cep")
# gamCyg = FixedTarget.from_name("gam Cyg")
# betageuze = FixedTarget.from_name("alp Ori")
# the1OriC = FixedTarget.from_name("the1 Ori C")
# oxAur = FixedTarget.from_name("OX Aur")

star = input('Geben Sie den Namen des Objektsterns ein: ')
star_ra = FixedTarget.from_name(star).ra.value
star_dec = FixedTarget.from_name(star).dec.value
