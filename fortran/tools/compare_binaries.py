import numpy as np
import os,sys
import shutil
import struct

##
## NB run from 'run' directory !!
##
w2d   = os.getenv('WIM2D_PATH')
dd    = w2d+'/fortran'
sys.path.append(dd+"/bin")
sys.path.append(dd+"/py_funs")

import fns_get_data  as Fdat
import fns_plot_data as Fplt

##########################################################################
mdir  = w2d+'/matlab/main/'
fdir  = w2d+'/fortran/run/'
cdir  = w2d+'/CXX/'
if 1:
   # odir  = [cdir+'outputs',fdir+'out_io']
   odir  = [cdir+'outputs',fdir+'out']
elif 0:
   odir  = [fdir+'out',fdir+'out_io']
   # odir  = [fdir+'out',fdir+'out_2']
   # odir  = [fdir+'out_io',fdir+'out_2']
elif 1:
   # can compare pure matlab results with results saved from mex functions
   # odir  = [mdir+'out_2',mdir+'out_io']
   # odir  = [mdir+'m_out',mdir+'out_io']
   odir  = [mdir+'m_out',mdir+'out_2']
elif 1:
   # matlab vs py interfaces
   odir  = [fdir+'out_io',mdir+'out_io']

print('Comparing directories:')
print(odir)
print('\n')
##########################################################################

##########################################################################
OPT   = 3   # 1: initial conditions; 2: final results; 3: prog results
if OPT==1:
   print("\n************************************************")
   print("Checking initial conditions...")
   print("************************************************\n")
elif OPT==2:
   print("\n************************************************")
   print("Checking final results...")
   print("************************************************\n")
elif OPT==3:
   print("\n************************************************")
   print("Checking progress files...")
   print("************************************************\n")
   n_prog   = 5
##########################################################################

##########################################################################
arrays   = 2*[0]
for j in range(2):
   outdir   = odir[j]
   bindir   = outdir+'/binaries'

   if OPT==1:

      ##########################################################################
      # Look at initial fields:
      ice_fields,wave_fields  = Fdat.fn_check_init(bindir) # load initial conditions from binaries
      arrays[j]               = ice_fields
      arrays[j].update(wave_fields)
      ##########################################################################

   elif OPT==2:

      ################################################################
      # Look at end results:
      out_fields  = Fdat.fn_check_out_bin(bindir)
      arrays[j]   = out_fields
      ################################################################

   elif OPT==3:

      ################################################################
      # Plot progress files (if they exist)
      # - as colormaps
      prog_files  = os.listdir(bindir+'/prog')
      steps       = []
      for pf in prog_files:
         if '.a'==pf[-2:]:
            stepno   = pf.strip('wim_prog').strip('.a')
            steps.append(stepno)

      print('Available time steps:')
      print(steps)

      stepno   = steps[n_prog]
      print("Checking results at time step "+stepno+" ...")
      print(outdir)
      prog_fields = Fdat.fn_check_prog(outdir,stepno) # dictionary eg {'Hs':Hs_array,...}
      arrays[j]   = prog_fields
      #############################################################

keys  = arrays[0].keys()
for key in keys:
   print('Comparing field: '+key)
   diff  = abs(arrays[0][key]-arrays[1][key])
   print('Maximum difference = '+str(np.max(diff))+'\n')
