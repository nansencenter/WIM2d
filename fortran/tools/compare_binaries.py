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

import run_WIM2d     as Rwim
import fns_get_data  as Fdat
import fns_plot_data as Fplt

##########################################################################
mdir  = '../../matlab/main/'
if 0:
   odir  = ['out','out_io']
elif 0:
   odir  = ['out','out_2']
elif 1:
   # can compare matlab results saved with mex functions also
   odir  = ['out',mdir+'out_2']
else:
   # can compare matlab results saved with mex functions also
   odir  = ['out',mdir+'out_io']

print('Comparing directories:')
print(odir)
print('\n')
##########################################################################

##########################################################################
OPT   = 2   # 1: initial results; 2: final results; 3: prog results
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
   n_prog   = 1
##########################################################################

##########################################################################
arrays   = 2*[0]
for j in range(2):
   bindir   = odir[j]+'/binaries'

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
      prog_fields = Fdat.fn_check_prog(outdir,stepno)
      arrays[j]   = prog_fields
      #############################################################

keys  = arrays[0].keys()
for key in keys:
   print('Comparing field: '+key)
   diff  = abs(arrays[0][key]-arrays[1][key])
   print('Maximum difference = '+str(np.max(diff))+'\n')
