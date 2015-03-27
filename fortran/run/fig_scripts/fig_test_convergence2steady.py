import numpy as np
import os
import sys
import shutil
import struct
from matplotlib import pyplot as plt
# import matplotlib.rcsetup as rc

##
## NB run from 'run' directory !!
##
dd   = os.path.abspath("..")
sys.path.append(dd+"/bin")
sys.path.append(dd+"/py_funs")

# interface to fortran code
# - compile with single frequency
import run_WIM2d     as Rwim 

# other python modules
import WIM2d_f2py    as Mwim
import fns_get_data  as Fdat
import fns_plot_data as Fplt

# steady state results to compare to:
import fns_boltzmann_steady   as Fbs

# RUN_OPT  = 2 # rerun then plot
RUN_OPT  = 3 # plot saved results

gf          = Fdat.fn_check_grid('inputs')
gfl         = gf['LANDMASK']
ICEMASK     = 1.-gfl
WAVEMASK    = 1.-gfl
grid_prams  = gf

###########################################################################
# set inputs: (icec,iceh,dfloe), (Hs,Tp,mwd)
if 1:

   # ice band of finite width
   xe                   = -220.e3
   ICEMASK              = 1+0*gf['X']
   ICEMASK[gf['X']<xe]  = 0.
   ICEMASK[gfl>0]       = 0.

   # right hand ice edge
   ice_width                        = 100.e3
   ICEMASK[gf['X']>(xe+ice_width)]  = 0.

   # edge of wave mask
   xw                   = -260.e3
   WAVEMASK             = 1+0*gf['X']
   WAVEMASK[gf['X']>xw] = 0.
   WAVEMASK[gfl>0]      = 0.

   in_fields   = {'icec':.7*ICEMASK,'iceh':2.*ICEMASK,'dfloe':100.*ICEMASK,
                  'Hs':3.*WAVEMASK,'Tp':12.*WAVEMASK,'mwd':-90.*WAVEMASK}

int_prams   = None # default integer parameters
real_prams  = None # default real parameters

if 1:
   # change integer parameters:
   SOLVER      = 1
   ADV_DIM     = 1
   CHECK_FINAL = 1
   CHECK_PROG  = 1
   CHECK_INIT  = 1
   DO_BREAKING = 0 # no breaking - testing convergence of Hs to steady-state
   int_prams   = np.array([SOLVER,ADV_DIM,
                           CHECK_FINAL,CHECK_PROG,CHECK_INIT,
                           DO_BREAKING])

if 1:
   # change real parameters:
   young          = 5.0e9
   visc_rp        = 0.0
   duration_hours = 12.0
   duration       = duration_hours*60*60
   real_prams     = np.array([young,visc_rp,duration])


out_fields,outdir = Rwim.do_run(RUN_OPT=RUN_OPT,in_fields=in_fields,
                                    int_prams=int_prams,
                                    real_prams=real_prams)

##########################################################################
if 0:
   # get steady state solution # TODO get this working
   out   = Fbs.solve_boltzmann_ft(width=width,
            alp=alp,N=N,alp_dis=alp_dis,cg=cg,f_inc=Fbs.dirspec_inc_spreading,Hs=Hs)
##########################################################################

##########################################################################
# Make plots
bindir   = outdir+'/binaries'
figdir   = 'fig_scripts/figs/TC2S'
if not os.path.exists(figdir):
   os.mkdir(figdir)



##########################################################################
grid_prams  = Fdat.fn_check_grid(bindir) # load grid from binaries
xx          = 1.e-3*grid_prams['X'][:,0]
labs        = ['x, km','$H_s$, m']

steps2plot  = range(0,300,40)
cols        = ['k','b','r','g',
               'k','b','r','g',
               'k','b','r','g']
fig         = None

loop_i   = -1
for nstep in steps2plot:
   out_fields  = Fdat.fn_check_prog(outdir,nstep) # load ice/wave conditions from binaries
   Hs_n        = out_fields['Hs']
   #
   loop_i   = loop_i+1
   fig      = Fplt.plot_1d(xx,Hs_n,labs=labs,f=fig,color=cols[loop_i])
#
figname  = figdir+'/convergence2steady.png'
print('saving to '+figname+'...')
plt.savefig(figname,bbox_inches='tight',pad_inches=0.05)
plt.close()
fig.clf()

sys.exit('fig_test_convergence2steady.py, L94')
##########################################################################

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
   figdir3     = figdir+'/prog'
   prog_files  = os.listdir(bindir+'/prog')
   steps       = []
   for pf in prog_files:
      if '.a'==pf[-2:]:
         stepno   = pf[-5:-2]
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
      prog_fields = Fdat.fn_check_prog(outdir,int(stepno))
      figdir3_0   = figdir3+'/'+stepno
      Fplt.fn_plot_final(grid_prams,prog_fields,figdir3_0)
      print("Plots in "+figdir3_0+'\n')
   ################################################################
