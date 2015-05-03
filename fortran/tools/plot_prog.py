import numpy as np
import os
import sys
import shutil
import struct
# import matplotlib.rcsetup as rc

wim2d_path=os.getenv('WIM2D_PATH')
dd   = os.path.abspath(wim2d_path+"/fortran")
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

##########################################################################
if 1:
   # Look at initial fields:
   print("Plotting initial conditions...")
   ice_fields,wave_fields  = Fdat.fn_check_init(bindir) # load initial conditions from binaries
   #
   figdir1  = figdir+'/init/'
   Fplt.fn_plot_init(grid_prams,ice_fields,wave_fields,figdir1) # plot initial conditions
   print("Plots in "+figdir+"/init")
   print(" ")
##########################################################################

################################################################
if 1:
   # Look at end results:
   print("Plotting final results...")
   figdir2     = figdir+'/final/'
   out_fields  = Fdat.fn_check_out_bin(bindir)
   Fplt.fn_plot_final(grid_prams,out_fields,figdir2)
   print("Plots in "+figdir2+'\n')
   print(" ")
################################################################

if 1:
   ################################################################
   # Plot progress files (if they exist)
   prog_files  = os.listdir(bindir+'/prog')
   steps       = []
   for pf in prog_files:
      if '.a'==pf[-2:]:
         stepno   = pf[-5:-2]
         steps.append(stepno)

   Nprogs   = len(steps)
   print(str(Nprogs)+' sets of progress files to plot\n')

   if 1:
      # plot as colormaps
      figdir3     = figdir+'/prog'

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

   elif 1:
      ################################################################
      # Plot as profiles
      steps2plot  = range(0,600,40)
      # steps2plot  = range(0,140,40)
      cols     = ['k','b','r','g','m','c']
      lstil    = ['-','--','-.',':']
      Nc       = len(cols)
      loop_c   = -1
      loop_s   = 0

      for nstep in steps2plot:
         out_fields  = Fdat.fn_check_prog(outdir,nstep) # load ice/wave conditions from binaries
         Hs_n        = out_fields['Hs']
         #
         if loop_c==Nc-1:
            loop_c   = 0
            loop_s   = loop_s+1
         else:
            loop_c   = loop_c+1

         fig      = Fplt.plot_1d(xx,Hs_n,labs=labs,f=fig,color=cols[loop_c],linestyle=lstil[loop_s])
      #
      figname  = figdir+'/convergence2steady.png'
      print('saving to '+figname+'...')
      plt.savefig(figname,bbox_inches='tight',pad_inches=0.05)
      plt.close()
      fig.clf()
