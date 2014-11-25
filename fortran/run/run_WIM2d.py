import numpy as np
import os
import sys
import struct
# import matplotlib.rcsetup as rc

dd   = os.path.abspath("..")
sys.path.append(dd+"/Build")
sys.path.append(dd+"/misc_py")

import WIM2d_f2py as Mwim
import fns_get_data as Fdat
import fns_plot_data as Fplt


################################################################
def do_run(RUN_OPT=0,in_fields=None):
   run_dict = {0: 'old version (no in/out)',
               1: 'look at saved results of no in/out run',
               2: 'in/out',
               3: 'look at saved results of in/out run'}
   print("###################################################")
   print("Run option:")
   print(" "+str(RUN_OPT)+' : '+run_dict[RUN_OPT])
   print("###################################################")

   # check directories for outputs exist
   if (RUN_OPT < 2):
      outdir   = 'out'
   else:
      outdir   = 'out_io'

   figdir   = outdir+'/figs'
   dirs  = [outdir,outdir+'/log',
            outdir+'/binaries',outdir+'/binaries/prog',
            figdir,figdir+'/init',figdir+'/final']
   for j in range(0,len(dirs)):
      dirj  = dirs[j]
      if not os.path.exists(dirj):
         os.makedirs(dirj)

   if RUN_OPT%2 is 0:
      # clear out old progress files
      dd    = os.path.abspath(outdir+"/binaries/prog")
      files = os.listdir(dd)
      for f in files:
         os.remove(dd+"/"+f)
   else:
      # load out_fields
      print(" ")
      print("###################################################")
      print("Reading results from: "+outdir)
      print("###################################################")
      print(" ")
      out_fields  = Fdat.fn_check_out_bin(outdir)

   #############################################################
   if RUN_OPT is 0:
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

      # load results from binary files
      out_fields  = Fdat.fn_check_out_bin(outdir)

   #############################################################
   elif RUN_OPT is 2:
      # run wim2d with inputs and outputs

      if not (in_fields is None):
         # 'in_fields' is given as input
         # - put data into 'in_arrays':
         keys  = ['icec','iceh','dfloe','Hs','Tp','mwd']

         nx,ny       = in_fields[keys[0]].shape
         in_arrays   = np.zeros((nx,ny,6))

         n  = 0
         for key in keys:
            in_arrays[:,:,n]  = in_fields[key]
            n                 = n+1

         del in_fields

      elif 1:
         # 'in_fields' not given as input
         # - read in inputs from saved files:
         in_dir                  = 'out'
         grid_prams              = Fdat.fn_check_grid(in_dir)
         ice_fields,wave_fields  = Fdat.fn_check_init(in_dir)

         # merge ice and wave fields:
         ice_fields.update(wave_fields)
         in_fields   = ice_fields
         del wave_fields

         # put data into 'in_arrays':
         keys  = ['icec','iceh','dfloe','Hs','Tp','mwd']

         nx,ny       = in_fields[keys[0]].shape
         in_arrays   = np.zeros((nx,ny,6))

         n  = 0
         for key in keys:
            in_arrays[:,:,n]  = in_fields[key]
            n                 = n+1

         del in_fields,ice_fields

      elif 1:
         # 'in_fields' not given as input
         # - specify 'in_arrays' manually

         GRID_OPT       = 1
         nx,ny          = Mwim.get_grid_size()
         X,Y,LANDMASK   = Mwim.retrieve_grid(GRID_OPT,nx,ny)
         xmin           = X.min()
         xmax           = X.max()

         # set ice conc/thickness
         ICE_MASK = np.zeros((nx,ny))
         ICE_MASK[np.logical_and(X>0.7*xmin,LANDMASK<1)]  = 1 # i>=24
         icec  = 0.75*ICE_MASK
         iceh  = 2.0*ICE_MASK
         dfloe = 250*ICE_MASK

         # set wave fields
         WAVE_MASK   = np.zeros((nx,ny))
         WAVE_MASK[X<xmin*0.8]   = 1   # i<=15
         Hs    = 2.0*WAVE_MASK
         Tp    = 12.0*WAVE_MASK
         mwd   = -90.0*WAVE_MASK

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

      # Dmax  = out_arrays[:,:,0]
      # Hs    = out_arrays[:,:,1]
      # Tp    = out_arrays[:,:,2]
      # tau_x = out_arrays[:,:,3]
      # tau_y = out_arrays[:,:,4]

      # convert out_arrays to Out_Fields object
      out_fields  = Fdat.fn_check_out_arr(out_arrays)
      del out_arrays
   #############################################################

   return out_fields,outdir
################################################################

if 0:
   RUN_OPT     = 2
   out_fields  = do_run(RUN_OPT)
elif 1:
   # check passing in of 'in_fields'
   # - read in inputs from saved files:
   in_dir                  = 'out'
   grid_prams              = Fdat.fn_check_grid(in_dir)
   ice_fields,wave_fields  = Fdat.fn_check_init(in_dir)

   # merge ice and wave fields:
   ice_fields.update(wave_fields)
   in_fields   = ice_fields

   out_fields  = do_run(RUN_OPT=2,in_fields=in_fields)
elif 0:
   out_fields  = do_run(0)
   out_fields2 = do_run(2)
elif 0:
   out_fields  = do_run(1)
   out_fields2 = do_run(3)

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
   Fplt.fn_plot_final(grid_prams,out_fields,figdir) # plot TODO - change fn from out_arrays to out_fields
   print("Plots in "+figdir+"/final")
elif 1:
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
