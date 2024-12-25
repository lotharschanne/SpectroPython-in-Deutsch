#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Das Skript erzeugt eine Animation aus png's (z.B. Spektrenplots).

20240130
@author: lothar
"""

from PIL import Image
import glob


files = input('Geben Sie den Namen der png-Bilder an: ')
imgs = glob.glob(files)
imgs.sort()

frames = []

for i in imgs:
    frames.append(Image.open(i))

frames[0].save('animation.gif', format='GIF',
               append_images=frames[1:],
               save_all=True,
               duration=300, loop=0)  # duration anpassen
