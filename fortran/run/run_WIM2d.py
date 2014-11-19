import numpy as np
import os
import sys
import matplotlib.pyplot as plt

dd   = os.path.abspath("..")
sys.path.append(dd+"/Build")

import WIM2d_f2py as Mwim

# check directories for outputs exist
dirs  = ['out','log','prog']
for j in range(0,len(dirs)):
   dirj  = dirs[j]
   if not os.path.exists(dirj):
      os.makedirs(dirj)

# clear out old progress files
dd    = os.path.abspath("prog")
files = os.listdir("prog")
for f in files:
   os.remove(dd+"/"+f)

if 0:
   # run the "dumb" WIM
   print(" ")
   print("###################################################")
   print("Running WIM without input/output")
   print("###################################################")
   print(" ")

   Mwim.wim2d_run()
   
   print(" ")
   print("###################################################")
   print("Finished call to wim2d_run:")
   print("###################################################")
   print(" ")
else:
   # run it with inputs and outputs
   GRID_OPT       = 1
   nx,ny          = Mwim.get_grid_size()
   X,Y,LANDMASK   = Mwim.retrieve_grid(GRID_OPT,nx,ny)
   xmin           = X.min()
   xmax           = X.max()

   # set ice conc/thickness
   ICE_MASK = np.zeros((nx,ny))
   ICE_MASK[np.logical_and(X>0.7*xmin,LANDMASK<1)]  = 1 # i>=24
   icec  = 0.7*ICE_MASK
   iceh  = 2.0*ICE_MASK
   dfloe = 250*ICE_MASK

   # set wave fields
   WAVE_MASK   = np.zeros((nx,ny))
   WAVE_MASK[X<xmin*0.8]   = 1   # i<=15
   Hs    = 2.0*WAVE_MASK
   Tp    = 12.0*WAVE_MASK
   mwd   = -90.0*WAVE_MASK

   in_arrays   = np.zeros((nx,ny,6))
   in_arrays[:,:,0]  = icec
   in_arrays[:,:,1]  = iceh
   in_arrays[:,:,2]  = dfloe
   in_arrays[:,:,3]  = Hs
   in_arrays[:,:,4]  = Tp
   in_arrays[:,:,5]  = mwd

   # run the WIM
   print(" ")
   print("###################################################")
   print("Running WIM with input/output")
   print("###################################################")
   print(" ")
   
   out_arrays  = Mwim.wim2d_run_io(in_arrays)

   print(" ")
   print("###################################################")
   print("Finished call to wim2d_run_io:")
   print("###################################################")
   print(" ")

   Dmax  = out_arrays[:,:,0]
   Hs    = out_arrays[:,:,1]
   Tp    = out_arrays[:,:,2]
   tau_x = out_arrays[:,:,3]
   tau_y = out_arrays[:,:,4]
