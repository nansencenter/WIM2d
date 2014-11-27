import numpy as np
import os
import sys
import struct
import matplotlib.pyplot as plt
# import matplotlib.rcsetup as rc

##
## NB run from 'run' directory !!
##
dd    = os.path.abspath("..")
dirs  = [dd+"/bin",dd+"/run",dd+"/misc_py"]
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
v        = np.arange(0,ny)/float(ny-1) # ny points between 0,1
Tp0      = 8.0
Tp1      = 18.0
Tp_vec   = Tp0+(Tp1-Tp0)*v

if 0:
   # do a new calculation:

   ######################################################
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
   mwd   = -90.0*WAVE_MASK

   if 0:
      Tp = 12.0*WAVE_MASK
   else:
      Tp = np.zeros((nx,ny))
      for j in range(0,ny):
         Tp[:,j]  = Tp_vec[j]*WAVE_MASK[:,j]
   ######################################################

   # set in_fields:
   in_fields   = {'icec'   : icec,
                  'iceh'   : iceh,
                  'dfloe'  : dfloe,
                  'Hs'     : Hs,
                  'Tp'     : Tp,
                  'mwd'    : mwd}

   # parameters for advection/attenuation
   SOLVER      = 1
   ADV_DIM     = 1
   int_prams   = np.array([SOLVER,ADV_DIM])

   # do calculation in fortran:
   out_fields,outdir = Rwim.do_run(RUN_OPT=2,
                                   in_fields=in_fields,
                                   int_prams=int_prams)
   ######################################################
else:
   # load results of previous run:
   out_fields,outdir = Rwim.do_run(RUN_OPT=3)
#########################################################

#########################################################
# calculate maxima of tau_x and get Wmiz:
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
#########################################################

#########################################################
def plot_diagnostics(Tp_vec,Wmiz,taux_max):
   # plot MIZ width and taux max

   figdir   = 'out_io/figs_diag'

   # MIZ width:
   fig   = figdir+'/Wmiz.png'
   f     = Fplt.plot_1d(Tp_vec,Wmiz,['$T_p$, s','$W_{MIZ}$, km'])
   plt.savefig(fig,bbox_inches='tight',pad_inches=0.05)
   plt.close()
   f.clf()

   # tau_x
   fig   = figdir+'/taux.png'
   f     = Fplt.plot_1d(Tp_vec,taux_max,['$T_p$, s','stress ($x$ dir), Pa'])
   plt.savefig(fig,bbox_inches='tight',pad_inches=0.05)
   plt.close()
   f.clf()
#########################################################

if 1:
   print("Tp (s) :")
   print(Tp_vec)

   print(" ")
   print("MIZ midth (km) :")
   print(Wmiz)

   print(" ")
   print("max tau_x (Pa) :")
   print(taux_max)

   print(" ")
   print("Plotting diagnostics...")
   print(" ")
   plot_diagnostics(Tp_vec,Wmiz,taux_max)
#########################################################

if 1:
   # plot results from binaries:
   print(" ")
   print("Getting inputs/outputs from binaries...")
   print(" ")

   ## look at initial fields:
   print("Plotting initial conditions...")
   grid_prams              = Fdat.fn_check_grid(outdir) # load grid from binaries
   ice_fields,wave_fields  = Fdat.fn_check_init(outdir) # load initial conditions from binaries
   ##
   figdir   = outdir+'/figs/'
   Fplt.fn_plot_init(grid_prams,ice_fields,wave_fields,figdir) # plot initial conditions
   print("Plots in "+figdir+"/init")
   print(" ")

   ## look at results:
   print("Plotting results...")
   out_fields  = Fdat.fn_check_out_bin(outdir)
   # Fplt.fn_plot_final(grid_prams,out_fields,figdir)
   Fplt.fn_plot_final_V1d(grid_prams,Tp_vec,out_fields,figdir)
   print("Plots in "+figdir+"/final")
elif 0:
   # plot results from outputs:

   ## look at initial fields:
   print("Plotting initial conditions...")
   figdir   = outdir+'/figs/'
   Fplt.fn_plot_init(grid_prams,ice_fields,wave_fields,figdir) # plot initial conditions
   print("Plots in "+figdir+"/init")
   print(" ")

   ## look at results:
   print("Plotting results...")
   # Fplt.fn_plot_final(grid_prams,out_fields,figdir)
   Fplt.fn_plot_final_V1d(grid_prams,Tp_vec,out_fields,figdir)
   print("Plots in "+figdir+"/final")
