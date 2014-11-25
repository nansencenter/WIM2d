import numpy as np
import os
import sys
import struct
# import matplotlib.rcsetup as rc

##
## NB run from 'run' directory !!
##
dd    = os.path.abspath("..")
dirs  = [dd+"/Build",dd+"/run",dd+"/misc_py"]
for n in range(0,len(dirs)):
   dd2   = dirs[n]
   print('adding path : '+dd2)
   sys.path.append(dd2)

import WIM2d_f2py    as Mwim
import run_WIM2d     as Rwim
import fns_get_data  as Fdat
import fns_plot_data as Fplt

# get grid
gdir        = 'out'
grid_prams  = Fdat.fn_check_grid(gdir)
nx          = grid_prams['nx']
ny          = grid_prams['ny']
dx          = grid_prams['dx']
dy          = grid_prams['dy']
X           = grid_prams['X']
Y           = grid_prams['Y']
LANDMASK    = grid_prams['LANDMASK']
xmin        = X.min()
xmax        = X.max()

#########################################################
if 1:
   # set input ice conc/thickness
   ICE_MASK = np.zeros((nx,ny))
   ICE_MASK[np.logical_and(X>0.7*xmin,LANDMASK<1)]  = 1 # i>=24
   icec  = 0.75*ICE_MASK
   iceh  = 2.0*ICE_MASK
   dfloe = 250*ICE_MASK

   # set input wave fields
   WAVE_MASK   = np.zeros((nx,ny))
   WAVE_MASK[X<xmin*0.8]   = 1   # i<=15
   Hs    = 2.0*WAVE_MASK
   Tp    = 12.0*WAVE_MASK
   mwd   = -90.0*WAVE_MASK
#########################################################

in_fields   = {'icec'   : icec,
               'iceh'   : iceh,
               'dfloe'  : dfloe,
               'Hs'     : Hs,
               'Tp'     : Tp,
               'mwd'    : mwd}
out_fields,outdir = Rwim.do_run(RUN_OPT=2,in_fields=in_fields)

taux_max = np.zeros(ny)
Wmiz     = np.zeros(ny)
for j in range(0,ny):
   Dmax  = out_fields['dfloe'][:,j]
   taux  = out_fields['taux'] [:,j]
   ##
   taux_max[j] = taux.max() # max stress in Pa
   ##
   miz   = np.zeros(nx)
   miz[np.logical_and(Dmax>0.0,Dmax<250.0)]  = 1.0
   Wmiz[j]  = sum(miz*dx/1.0e3)  # MIZ width in km

print(Wmiz)
print(taux_max)
