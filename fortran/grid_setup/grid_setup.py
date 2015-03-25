import fns_grid_setup as gs
import os
import sys
import matplotlib.pyplot as plt

# Set grid configuration (see fns_grid_setup.py)
if 1:
   GRID_OPT = 0 # standard 1d grid configuration
else:
   GRID_OPT = 1 # standard 2d grid configuration

TEST     = 0 # if 1: don't save grid files to ../run/inputs 
             # if 0: save grid files to test/out_py
             # *Both options make a plot of the LAND MASK
LAND_OPT = 2
grid_fields,grid_arrays = gs.grid_setup(GRID_OPT=GRID_OPT,TEST=TEST,LAND_OPT=LAND_OPT)

if 1:
   print('**************************************************')
   print('Testing values of grid_fields\n')
   keys  = grid_fields.keys()

   for key in keys:
      fld   = grid_fields[key]
      if hasattr(fld,'shape'):
         # array:
         print('field: '+key)
         print('field min: '+str(fld.min()))
         print('field max: '+str(fld.max())+'\n')
      else:
         # scalar:
         print(key+' = '+str(fld)+'\n')

   print('**************************************************')
