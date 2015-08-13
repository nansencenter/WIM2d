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

RUN_OPT  = 2 # rerun then plot
# RUN_OPT  = 3 # plot saved results

gf          = Fdat.fn_check_grid('inputs')
gfl         = gf['LANDMASK']
ICEMASK     = 1.-gfl
WAVEMASK    = 1.-gfl
grid_prams  = gf

###########################################################################
# set inputs: (icec,iceh,dfloe), (Hs,Tp,mwd)
if 1:

   if 1:
      # ice edge
      xe       = .5*(gf['X'].min()+gf['X'].max())\
                  -.7*.5*(-gf['X'].min()+gf['X'].max())
      ICEMASK  = 1+0*gf['X']
      #
      ICEMASK[gf['X']<xe]  = 0.
      ICEMASK[gfl>0]       = 0.
      #
      c_in  = .7
      h_in  = 2.
      D_in  = 300.
   else:
      # strip
      strip_width = 100.e3
      xe          = .5*(gf['X'].min()+gf['X'].max())\
                     -.7*.5*(-gf['X'].min()+gf['X'].max())
      ICEMASK     = 1+0*gf['X']
      #
      ICEMASK[abs(gf['X'])<xe]               = 0.
      ICEMASK[abs(gf['X'])>xe+strip_width]   = 0.
      ICEMASK[gfl>0]                         = 0. # 0 on land
      #
      c_in  = .7
      h_in  = 2.
      D_in  = 100.

   ice_fields0 = {'icec':c_in*ICEMASK,'iceh':h_in*ICEMASK,'dfloe':D_in*ICEMASK}

   # edge of wave mask
   xw                   = .5*(gf['X'].min()+gf['X'].max())\
                           -.8*.5*(-gf['X'].min()+gf['X'].max())
   WAVEMASK             = 1+0*gf['X']
   WAVEMASK[gf['X']>xw] = 0.
   WAVEMASK[gfl>0]      = 0.

   Hs_in          = 2.
   Tp_in          = 12.
   mwd_in         = -90.
   wave_fields0   = {'Hs':Hs_in*WAVEMASK,'Tp':Tp_in*WAVEMASK,'mwd':mwd_in*WAVEMASK}

int_prams   = None # default integer parameters
real_prams  = None # default real parameters

if 1:
   # change integer parameters:
   SCATMOD     = 1
   ADV_DIM     = 1
   ADV_OPT     = 2
   CHECK_FINAL = 1
   CHECK_PROG  = 1
   CHECK_INIT  = 1
   DO_BREAKING = 1
   STEADY      = 1
   int_prams   = np.array([SCATMOD,ADV_DIM,ADV_OPT,
                           CHECK_FINAL,CHECK_PROG,CHECK_INIT,
                           DO_BREAKING,STEADY])

if 1:
   # change real parameters:
   young          = 5.49e9
   visc_rp        = 0.0
   duration_hours = 6.0
   duration       = duration_hours*60*60
   CFL            = 0.7
   real_prams     = np.array([young,visc_rp,duration,CFL])

# call gateway between python and pre-compiled f2py module
out_fields,outdir = Rwim.do_run_vSdir(ice_fields=ice_fields0,wave_fields=wave_fields0,\
                                       int_prams=int_prams,real_prams=real_prams)

##########################################################################
# Make plots
bindir   = outdir+'/binaries'
figdir   = outdir+'/figs'

##########################################################################
# Look at initial fields:
print("Plotting initial conditions...")
grid_prams              = Fdat.fn_check_grid(bindir) # load grid from binaries
ice_fields,wave_fields  = Fdat.fn_check_init(bindir) # load initial conditions from binaries
#
figdir1  = figdir+'/init/'
Fplt.fn_plot_init(grid_prams,ice_fields,wave_fields,figdir1) # plot initial conditions
print("Plots in "+figdir+"/init")
print(" ")
##########################################################################

################################################################
# Look at end results:
print("Plotting results...")
figdir2  = figdir+'/final/'
Fplt.fn_plot_final(grid_prams,out_fields,figdir2)
print("Plots in "+figdir2+'\n')
print(" ")
################################################################

if 1:
   ################################################################
   # Plot progress files (if they exist)
   # - as colormaps
   figdir3     = figdir+'/prog'
   prog_files  = os.listdir(bindir+'/prog')
   steps       = []
   for pf in prog_files:
      if '.a'==pf[-2:]:
         stepno   = pf.strip('wim_prog').strip('.a')
         steps.append(stepno)

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
   ################################################################

   print(' ')
   print("To make a movie of progress images:")
   print("cd "+figdir+'/prog')
   print(dd+'/tools/prog2mp4.sh Hs (or Dmax,taux,tauy)')

elif 0:
   ################################################################
   # Plot progress files (if they exist)
   # - as profiles
   cols     = ['k','b','r','g','m','c']
   lstil    = ['-','--','-.',':']
   Nc       = len(cols)
   loop_c   = -1
   loop_s   = 0

   figdir3     = figdir+'/prog'
   prog_files  = os.listdir(bindir+'/prog')
   steps       = []
   for pf in prog_files:
      if '.a'==pf[-2:]:
         stepno   = pf.strip('wim_prog').strip('.a')
         steps.append(stepno)

   for step in steps:
      out_fields  = Fdat.fn_check_prog(outdir,step) # load ice/wave conditions from binaries
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
