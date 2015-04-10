import numpy as np
import os
import sys
import shutil
import struct
# import matplotlib.rcsetup as rc

##
## NB run from 'run' directory !!
##
dd   = os.path.abspath("..")
sys.path.append(dd+"/bin")
sys.path.append(dd+"/py_funs")

import run_WIM2d     as Rwim
import WIM2d_f2py    as Mwim
import fns_get_data  as Fdat
import fns_plot_data as Fplt

RUN_OPT  = 2 # rerun then plot
# RUN_OPT  = 3 # plot saved results

gf          = Fdat.fn_check_grid('inputs')
gfl         = gf['LANDMASK']
ICEMASK     = 1.-gfl
WAVEMASK    = 1.-gfl
grid_prams  = gf

###########################################################################
# set inputs: (icec,iceh,dfloe), (Hs,Tp,mwd)
if 1:

   # ice edge
   xe                   = -220.e3
   ICEMASK              = 1+0*gf['X']
   ICEMASK[gf['X']<xe]  = 0.
   ICEMASK[gfl>0]       = 0.

   # edge of wave mask
   xw                   = -260.e3
   WAVEMASK             = 1+0*gf['X']
   WAVEMASK[gf['X']>xw] = 0.
   WAVEMASK[gfl>0]      = 0.

   in_fields   = {'icec':.7*ICEMASK,'iceh':2.*ICEMASK,'dfloe':300.*ICEMASK,
                  'Hs':3.*WAVEMASK,'Tp':12.*WAVEMASK,'mwd':-90.*WAVEMASK}

int_prams   = None # default integer parameters
real_prams  = None # default real parameters

if 1:
   # change integer parameters:
   SCATMOD      = 1
   ADV_DIM     = 1
   CHECK_FINAL = 1
   CHECK_PROG  = 1
   CHECK_INIT  = 1
   DO_BREAKING = 1
   int_prams   = np.array([SCATMOD,ADV_DIM,
                           CHECK_FINAL,CHECK_PROG,CHECK_INIT,
                           DO_BREAKING])

if 1:
   # change real parameters:
   young          = 5.0e9
   visc_rp        = 0.0
   duration_hours = 12.0
   duration       = duration_hours*60*60
   real_prams     = np.array([young,visc_rp,duration])


out_fields,outdir = Rwim.do_run(RUN_OPT=RUN_OPT,in_fields=in_fields,
                                    int_prams=int_prams,
                                    real_prams=real_prams)

##########################################################################
# Make plots
bindir   = outdir+'/binaries'
figdir   = outdir+'/figs'

##########################################################################
# Look at initial fields:
print("Plotting initial conditions...")
grid_prams              = Fdat.fn_check_grid(bindir) # load grid from binaries
ice_fields,wave_fields  = Fdat.fn_check_init(bindir) # load initial conditions from binaries
#
figdir1  = figdir+'/init/'
Fplt.fn_plot_init(grid_prams,ice_fields,wave_fields,figdir1) # plot initial conditions
print("Plots in "+figdir+"/init")
print(" ")
##########################################################################

################################################################
# Look at end results:
print("Plotting results...")
figdir2  = figdir+'/final/'
Fplt.fn_plot_final(grid_prams,out_fields,figdir2)
print("Plots in "+figdir2+'\n')
print(" ")
################################################################

if 1:
   ################################################################
   # Plot progress files (if they exist)
   figdir3     = figdir+'/prog'
   prog_files  = os.listdir(bindir+'/prog')
   steps       = []
   for pf in prog_files:
      if '.a'==pf[-2:]:
         stepno   = pf[-5:-2]
         steps.append(stepno)

   # make dir for progress plots
   if (not os.path.exists(figdir3)) and len(steps)>0:
      os.mkdir(figdir3)

   # clear old progress plots
   old_dirs = os.listdir(figdir3)
   for od in old_dirs:
      # need this for macs
      if od!='.DS_Store':
         # os.rmdir(figdir3+'/'+od)
         shutil.rmtree(figdir3+'/'+od)

   for stepno in steps:
      print("Plotting results at time step "+stepno+" ...")
      prog_fields = Fdat.fn_check_prog(outdir,int(stepno))
      figdir3_0   = figdir3+'/'+stepno
      Fplt.fn_plot_final(grid_prams,prog_fields,figdir3_0)
      print("Plots in "+figdir3_0+'\n')
   ################################################################
