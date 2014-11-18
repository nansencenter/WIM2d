import numpy
import os
import sys

dd   = os.path.abspath("..")
sys.path.append(dd+"/Build")

import WIM2d_f2py
wim   = WIM2d_f2py.wim2d_f2py

dirs  = ['out','log']
for j in [0,1]:
   dirj  = dirs[j]
   if not os.path.exists(dirj):
      os.makedirs(dirj)

# run the WIM
wim()
