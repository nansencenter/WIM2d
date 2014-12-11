import numpy as np
import os
import sys

dd   = os.path.abspath("..")
sys.path.append(dd+"/bin")
sys.path.append(dd+"/misc_py")

import WIM2d_f2py as Mwim
import fns_get_data as Fdat
import fns_plot_data as Fplt

################################################################
def do_run(RUN_OPT=0,in_fields=None,int_prams=None):
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
      out_fields  = Fdat.fn_check_out_bin(outdir+'/binaries')

   #############################################################
   if RUN_OPT is 0:
      # run the "dumb" WIM
      print(" ")
      print("###################################################")
      print("Running WIM without input/output")
      print("###################################################")
      print(" ")

      Mwim.py_wim2d_run()
      
      print(" ")
      print("###################################################")
      print("Finished call to wim2d_run:")
      print("###################################################")
      print(" ")

      # load results from binary files
      out_fields  = Fdat.fn_check_out_bin(outdir+'/binaries')

   #############################################################
   elif RUN_OPT is 2:
      # run wim2d with inputs and outputs

      ##########################################################
      if int_prams is None:
         SOLVER      = 1
         ADV_DIM     = 2
         int_prams   = np.array([SOLVER,ADV_DIM])
      ##########################################################

      ##########################################################
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

      ##########################################################
      elif 1:
         # 'in_fields' not given as input
         # - read in inputs from saved files:
         in_dir                  = 'out/binaries'
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

      ##########################################################
      elif 1:
         # 'in_fields' not given as input
         # - specify 'in_arrays' manually

         in_dir      = 'out/binaries'
         grid_prams  = Fdat.fn_check_grid(in_dir)

         n        = grid_prams['nx']
         ny       = grid_prams['ny']
         X        = grid_prams['X']
         Y        = grid_prams['Y']
         LANDMASK = grid_prams['LANDMASK']
         xmin     = X.min()
         xmax     = X.max()

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
      ##########################################################

      # run the WIM
      print(" ")
      print("###################################################")
      print("Running WIM with input/output")
      print("###################################################")
      print(" ")
      
      out_arrays  = Mwim.py_wim2d_run_io(in_arrays,int_prams)

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
