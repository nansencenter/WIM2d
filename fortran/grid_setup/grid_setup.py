import fns_grid_setup as gs
import os
import sys
import matplotlib.pyplot as plt

CHANGE_GRID    = 1
CHANGE_WAVES   = 1

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

if CHANGE_GRID:
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

if CHANGE_WAVES:
   # change number of wave frequencies and directions
   nfreq = 1
   ndir  = 16
   hfil  = '../header_files/wave_info.h'

   print(' ')
   print('**************************************************')
   print('Changing number of wave frequencies and directions:')
   print('editing '+hfil)
   print('**************************************************')
   print(' ')

   """
   wave_info.h should look like this:
   integer,parameter :: n_wavdir    = 16
   integer,parameter :: n_wave_freq = 1!should be odd (Simpson's rule)
   """

   spc   = 6*' ' # F77 indent
   str1  = 'integer,parameter :: n_wavdir    = '
   str2  = 'integer,parameter :: n_wave_freq = '
   str3  = "!should be odd (Simpson's rule)"

   # write file
   hf = open(hfil,'w')
   hf.write(spc+str1+str(ndir)+'\n')
   hf.write(spc+str2+str(nfreq)+str3)
   hf.close()
