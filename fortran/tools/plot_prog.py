import numpy as np
import os
import sys
import shutil
import struct
import matplotlib.pyplot as plt

wim2d_path  = os.getenv('WIM2D_PATH')
dd          = os.path.abspath(wim2d_path+"/fortran")
sys.path.append(dd+"/bin")
sys.path.append(dd+"/py_funs")

import fns_get_data  as Fdat
import fns_plot_data as Fplt

# run from root results directory eg out, out_io
outdir   = os.getcwd()
bindir   = outdir+'/binaries'
figdir   = outdir+'/figs'

if not os.path.exists(bindir):
   raise ValueError('No binaries folder in current directory')
else:
   grid_prams  = Fdat.fn_check_grid(bindir)

##########################################################################
# Make plots
if not os.path.exists(figdir):
   os.mkdir(figdir)

PLOT_INIT   = 1
PLOT_FINAL  = 1
PLOT_PROG   = 1

##########################################################################
if PLOT_INIT:
   if not os.path.exists(bindir+'/wim_init.a'):
      print('Initial conditions not outputted: wim_init.[a,b]')
      print('- not plotting')
   else:
      # Look at initial fields:
      print("Plotting initial conditions...")
      figdir1  = figdir+'/init/'
      if not os.path.exists(figdir1):
         os.mkdir(figdir1)

      # new way (more general)
      fields,info = Fdat.fn_read_general_binary(bindir+'/wim_init.a')
      Fplt.fn_plot_gen(grid_prams,fields,figdir1)    # plot initial conditions

      print("Plots in "+figdir1)
      print(" ")
##########################################################################

################################################################
if PLOT_FINAL:
   if not os.path.exists(bindir+'/wim_init.a'):
      print('Final conditions not outputted: wim_out.[a,b]')
      print('- not plotting')
   else:
      # Look at end results:
      print("Plotting final results...")
      figdir2     = figdir+'/final/'
      if not os.path.exists(figdir2):
         os.mkdir(figdir2)

      fields,info = Fdat.fn_read_general_binary(bindir+'/wim_out.a')
      Fplt.fn_plot_gen(grid_prams,fields,figdir2)    # plot final results

      print("Plots in "+figdir2+'\n')
      print(" ")
################################################################

if PLOT_PROG:
   ################################################################
   # Plot progress files (if they exist)
   pdir        = bindir+'/prog'
   prog_files  = os.listdir(pdir)
   figdir3     = figdir+'/prog'
   if not os.path.exists(figdir3):
      os.mkdir(figdir3)

   # Typical limits for parameters
   # NB use the same key names as in the .b file
   zlims  = {'icec':[0,1],      \
             'iceh':[0,5],      \
             'Dmax':[0,300],    \
             'tau_x':[-.5,.5],  \
             'tau_y':[-.05,.05],\
             'Hs':[0,4],        \
             'Tp':[10,20],      \
             'mwd':[-180,180]}

   for pf in prog_files:
      if '.a'==pf[-2:]:
         stepno   = pf.strip('wim_prog').strip('.a')
         #
         afile    = pdir+'/'+pf
         print(afile)
         #
         fields,info = Fdat.fn_read_general_binary(afile)
         figdir3B    = figdir3+'/'+stepno
         Fplt.fn_plot_gen(grid_prams,fields,figdir3B,zlims_in=zlims)

   ################################################################
   print('\n**********************************************************************')
   print('to make movie, type')
   print(wim2d_path+'/fortran/tools/prog2mp4.sh Hs '+outdir+'/figs/prog')
   print('or')
   print(wim2d_path+'/fortran/tools/prog2mp4.sh Dmax '+outdir+'/figs/prog')
   print('**********************************************************************\n')
   ################################################################
