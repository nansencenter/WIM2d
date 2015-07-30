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

grid_prams  = Fdat.fn_check_grid(bindir)

##########################################################################
# Make plots
if not os.path.exists(figdir):
   os.mkdir(figdir)

PLOT_INIT   = 1
PLOT_FINAL  = 1
PLOT_PROG   = 1
OLD_WAY     = 0

##########################################################################
if PLOT_INIT:
   # Look at initial fields:
   print("Plotting initial conditions...")
   figdir1  = figdir+'/init/'
   if not os.path.exists(figdir1):
      os.mkdir(figdir1)

   if OLD_WAY:
      # old way
      ice_fields,wave_fields  = Fdat.fn_check_init(bindir) # load initial conditions from binaries
      Fplt.fn_plot_init(grid_prams,ice_fields,wave_fields,figdir1) # plot initial conditions
   else:
      # new way (more general)
      fields   = Fdat.fn_read_general_binary(bindir+'/wim_init.a')
      Fplt.fn_plot_gen(grid_prams,fields,figdir1)    # plot initial conditions

   print("Plots in "+figdir1)
   print(" ")
##########################################################################

################################################################
if PLOT_FINAL:
   # Look at end results:
   print("Plotting final results...")
   figdir2     = figdir+'/final/'
   if not os.path.exists(figdir2):
      os.mkdir(figdir2)

   if OLD_WAY:
      # old way
      out_fields  = Fdat.fn_check_out_bin(bindir)
      Fplt.fn_plot_final(grid_prams,out_fields,figdir2)
   else:
      fields   = Fdat.fn_read_general_binary(bindir+'/wim_out.a')
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

   steps       = []
   for pf in prog_files:
      if '.a'==pf[-2:]:
         stepno   = pf.strip('wim_prog').strip('.a')
         if not OLD_WAY:
            # only PROG_OPT==1 implemented
            afile    = pdir+'/'+pf
            print(afile)
            fields   = Fdat.fn_read_general_binary(afile)
            figdir3B = figdir3+'/'+stepno
            Fplt.fn_plot_gen(grid_prams,fields,figdir3B)    # plot final results
         else:
            steps.append(stepno)

   if not OLD_WAY:
      ################################################################
      print('\n**********************************************************************')
      print('to make movie, go to figs/prog and type')
      print(wim2d_path+'/fortran/tools/prog2mp4.sh Hs')
      print('or')
      print(wim2d_path+'/fortran/tools/prog2mp4.sh Dmax')
      print('**********************************************************************\n')
      ################################################################
   else:
      Nprogs   = len(steps)
      print(str(Nprogs)+' sets of progress files to plot\n')

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
         prog_fields = Fdat.fn_check_prog(outdir,stepno)
         figdir3_0   = figdir3+'/'+stepno
         Fplt.fn_plot_final(grid_prams,prog_fields,figdir3_0)
         print("Plots in "+figdir3_0+'\n')

      print('\n**********************************************************************')
      print('to make movie, go to figs/prog and type')
      print(wim2d_path+'/fortran/tools/prog2mp4.sh Hs')
      print('or')
      print(wim2d_path+'/fortran/tools/prog2mp4.sh Dmax')
      print('**********************************************************************\n')
      ################################################################
