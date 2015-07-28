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

PLOT_PROG   = 1
PROG_OPT    = 1

##########################################################################

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
      steps.append(stepno)

if PROG_OPT==1:
   ################################################################
   # Plot as profiles on different graphs (eg test advection)
   cols     = ['k']
   lstil    = ['-']

   figdir3     = figdir+'/prog'
   first_call  = True
   xx          = grid_prams['X'][:,0]
   labs        = ['$x$, km','$H_s$, m']

   # make dir for progress plots
   if (not os.path.exists(figdir3)) and len(steps)>0:
      os.mkdir(figdir3)

   for stepno in steps:
      print("Plotting results at time step "+stepno+" ...")
      prog_fields = Fdat.fn_check_prog(outdir,stepno)
      figdir3_0   = figdir3+'/'+stepno

      Hs_n  = prog_fields['Hs'][:,0]
      if first_call:
         Hs_0        = prog_fields['Hs'][:,0]
         fig         = None
         first_call  = False

      fig   = Fplt.plot_1d(xx/1.e3,Hs_n,labs=labs,f=fig,color='k',linestyle='-')
      if 0:
         # TODO: calc & plot "Hs_exact" from speed in x direction

         ################################################
         # set manually
         T        = 12. # wave period (s)
         wavdir0  = -90.# waves-from dirn (degrees)
         CFL      = 0.7
         ################################################

         ################################################
         dir0     = -np.pi/180.*(90.+wavdir0)   # waves-to dirn (radians)
         om       = 2*np.pi/T
         g        = 9.81
         k        = om*om/g
         cp       = om/k
         cg       = cp/2.
         ug       = cg*np.cos(dir0)
         #
         dx    = xx[1]-xx[0]
         dt    = CFL*dx/cg
         xt    = xx+ug*int(stepno)*dt
         fig   = Fplt.plot_1d(xt/1.e3,Hs_0,labs=labs,f=fig,color='r',linestyle='--')
         ################################################

      plt.xlim([xx.min()/1.e3,xx.max()/1.e3])
      plt.ylim([0.,2.*Hs_0.max()])

      if not os.path.exists(figdir3_0):
         os.mkdir(figdir3_0)
      figname  = figdir3_0+'/Hs.png'
      print('saving to '+figname+'...\n')
      plt.savefig(figname,bbox_inches='tight',pad_inches=0.05)
      plt.close()
      fig.clf()

   print('**********************************************************************')
   print('to make movie, go to figs/prog and type')
   print(wim2d_path+'/fortran/tools/prog2mp4.sh Hs')
   print('**********************************************************************\n')
   ################################################################

elif PROG_OPT==2:
   ################################################################
   # Plot as profiles on same graph (test convergence to steady)
   print(steps)
   steps2plot  = range(len(steps))
   #
   cols     = ['k','b','r','g','m','c']
   lstil    = ['-','--','-.',':']
   Nc       = len(cols)
   loop_c   = -1
   loop_s   = 0
   xx       = grid_prams['X'][:,0]
   fig      = None
   labs     = ['$x$, km','$H_s$, m']

   for nstep in steps2plot:
      out_fields  = Fdat.fn_check_prog(outdir,steps[nstep]) # load ice/wave conditions from binaries
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
