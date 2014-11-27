import numpy as np
import os
import sys
import struct
# import matplotlib.rcsetup as rc

##
## NB run from 'run' directory !!
##
dd   = os.path.abspath("..")
sys.path.append(dd+"/bin")
sys.path.append(dd+"/misc_py")

import run_WIM2d     as Rwim
import WIM2d_f2py    as Mwim
import fns_get_data  as Fdat
import fns_plot_data as Fplt

if 1:
   RUN_OPT     = 2
   out_fields  = Rwim.do_run(RUN_OPT)
elif 1:
   # check passing in of 'in_fields'
   # - read in inputs from saved files:
   in_dir                  = 'out'
   grid_prams              = Fdat.fn_check_grid(in_dir)
   ice_fields,wave_fields  = Fdat.fn_check_init(in_dir)

   # merge ice and wave fields:
   ice_fields.update(wave_fields)
   in_fields   = ice_fields

   out_fields  = Rwim.do_run(RUN_OPT=2,in_fields=in_fields)
elif 0:
   out_fields  = Rwim.do_run(0)
   out_fields2 = Rwim.do_run(2)
elif 0:
   out_fields  = Rwim.do_run(1)
   out_fields2 = Rwim.do_run(3)

if 0:
   # plot results
   ## look at initial fields:
   print("Plotting initial conditions...")
   grid_prams              = Fdat.fn_check_grid(outdir) # load grid from binaries
   ice_fields,wave_fields  = Fdat.fn_check_init(outdir) # load initial conditions from binaries
   ##
   Fplt.fn_plot_init(grid_prams,ice_fields,wave_fields,figdir) # plot initial conditions
   print("Plots in "+figdir+"/init")
   print(" ")

   ## look at results:
   print("Plotting results...")
   Fplt.fn_plot_final(grid_prams,out_fields,figdir)
   print("Plots in "+figdir+"/final")
elif 0:
   # compare binaries from different runs (wim2d_run & wim2d_io)
   # NB need same grid & initial conditions
   outdir1  = 'out'
   outdir2  = 'out_io'

   gp1   = Fdat.fn_check_grid(outdir1)
   gp2   = Fdat.fn_check_grid(outdir2)
   if 0:
      # check grids are the same
      print('Checking grids are the same...')
      ##
      keys  = ['X','Y','scuy','scvx','scp2','scp2i','LANDMASK']
      Key   = '         '
      for key in keys:
         diff     = np.abs(gp1[key]-gp2[key])
         diff_max = diff.max()
         diff_min = diff.min()
         #
         ss = ' max/min |difference|: %f %f' % (diff_max,diff_min)
         print(' '+key+Key[len(key):]+ss)
   elif 0:
      # check initial fields are the same
      print('Checking initial fields are the same...')
      ##
      if1,wf1  = Fdat.fn_check_init(outdir1)
      if2,wf2  = Fdat.fn_check_init(outdir2)

      Key   = '         '
      for key in if1.keys():
         diff     = np.abs(if1[key]-if2[key])
         diff_max = diff.max()
         diff_min = diff.min()
         #
         ss = ' max/min |difference|: %f %f' % (diff_max,diff_min)
         print(' '+key+Key[len(key):]+ss)

      for key in wf1.keys():
         diff     = np.abs(wf1[key]-wf2[key])
         diff_max = diff.max()
         diff_min = diff.min()
         #
         ss = ' max/min |difference|: %f %f' % (diff_max,diff_min)
         print(' '+key+Key[len(key):]+ss)

   elif 1:
      # check out fields are the same
      print('Checking final fields are the same...')
      ##
      of1   = Fdat.fn_check_out_bin(outdir1)
      of2   = Fdat.fn_check_out_bin(outdir2)

      Key   = '         '
      for key in of1.keys():
         diff     = np.abs(of1[key]-of2[key])
         diff_max = diff.max()
         diff_min = diff.min()
         #
         ss = ' max/min |difference|: %f %f' % (diff_max,diff_min)
         print(' '+key+Key[len(key):]+ss)

# elif 1:
#    # compare exponential decay of Hs for SOLVER = 1,0
# 
