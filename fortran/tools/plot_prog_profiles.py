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

PROG_OPT = 1
##########################################################################

################################################################
# Plot progress files (if they exist)
pdir        = bindir+'/prog'
prog_files  = os.listdir(pdir)
figdir3     = figdir+'/prog_profiles'
if not os.path.exists(figdir3):
   os.mkdir(figdir3)

if PROG_OPT==1:
   # plot each variable, save different fig for each time step
   for pf in prog_files:
      if '.a'==pf[-2:]:
         stepno   = pf.strip('wim_prog').strip('.a')
         # steps.append(stepno)
         #
         afile    = pdir+'/'+pf
         print(afile)
         fields,info = Fdat.fn_read_general_binary(afile)
         figdir3B = figdir3+'/'+stepno
         #
         field_profiles = {}
         for key in fields.keys():
            F  = fields[key][:,0]
            field_profiles.update({key:F})
         #
         grid_profiles = {}
         for key in grid_prams.keys():
            try:
               F  = grid_prams[key][:,0]
            except:
               F  = grid_prams[key]
            grid_profiles.update({key:F})

         grid_profiles['ny']  = 1
         Fplt.fn_plot_gen(grid_profiles,field_profiles,figdir3B)

   print('\n**********************************************************************')
   print('to make movie, type')
   print(wim2d_path+'/fortran/tools/prog2mp4.sh Hs figs/prog_profiles')
   print('or')
   print(wim2d_path+'/fortran/tools/prog2mp4.sh Dmax figs/prog_profiles')
   print('**********************************************************************\n')

elif PROG_OPT==2:
   ################################################################
   # Plot 1 variable as profiles on same graph (eg to test convergence to steady)
   #
   cols     = ['k','b','r','g','m','c']
   lstil    = ['-','--','-.',':']
   Nc       = len(cols)
   Ns       = len(lstil)
   loop_c   = -1
   loop_s   = -1
   xx       = grid_prams['X'][:,0]/1.e3
   if 1:
      vname = 'Hs'
      labs  = ['$x$, km','$H_s$, m']
   else:
      vname = 'tau_x'
      labs  = ['$x$, km',r'$\tau_x$, Pa']
   #
   fig   = plt.figure()
   ax    = fig.add_subplot(1,1,1)
   lines = []
   legs  = []


   for pf in prog_files:
      if '.a'==pf[-2:]:
         stepno   = pf.strip('wim_prog').strip('.a')
         afile    = pdir+'/'+pf
         print(afile)
         #
         out_fields,info   = Fdat.fn_read_general_binary(afile) # load ice/wave conditions from binaries
         Hs_n              = out_fields[vname][:,0]
         t_prog            = info['t_out']/3600.   # time in h
         leg               = '%5.1fh' % (t_prog)
         #
         loop_c   = np.mod(loop_c+1,Nc)
         if loop_c==0:
            loop_s   = np.mod(loop_s+1,Ns)

         fig,ax,line = Fplt.plot_1d(xx,Hs_n,labs=labs,pobj=[fig,ax],color=cols[loop_c],linestyle=lstil[loop_s])
         lines.append(line)
         legs.append(leg)
   #
   ax.legend(lines,legs)
   figname  = figdir+'/time_dep_'+vname+'.png'
   print('Saving to '+figname+'...')
   fig.savefig(figname,bbox_inches='tight',pad_inches=0.05)
   plt.close(fig)
