import numpy as np
import os,sys
import shutil
import struct
from matplotlib import pyplot as plt

w2d   = os.getenv('WIM2D_PATH')
sys.path.append("bin")

import run_WIM2d     as Rwim
import fns_get_data  as Fdat
import fns_plot_data as Fplt

TEST_MESH   = 1
DO_PLOTTING = 1

gf          = Fdat.fn_check_grid('grid')
gfl         = gf['LANDMASK']
ICEMASK     = 1.-gfl
WAVEMASK    = 1.-gfl
grid_prams  = gf

###########################################################################
# set inputs: (icec,iceh,dfloe), (Hs,Tp,mwd)
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
   h_in  = 1.
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
   h_in  = 1.
   D_in  = 100.

ice_fields0 = {'icec'   :c_in*ICEMASK,\
               'iceh'   :h_in*ICEMASK,\
               'dfloe'  :D_in*ICEMASK}
if 0:
   # test plot of inputs:
   fig   = plt.figure()
   ax1   = fig.add_subplot(3,1,1)
   ax2   = fig.add_subplot(3,1,2)
   ax3   = fig.add_subplot(3,1,3)
   ax1.plot(gf['X']/1.e3,ice_fields0['icec' ])
   ax2.plot(gf['X']/1.e3,ice_fields0['iceh' ])
   ax3.plot(gf['X']/1.e3,ice_fields0['dfloe'])
   plt.show(fig)
   sys.exit()

# edge of wave mask
xw                   = .5*(gf['X'].min()+gf['X'].max())\
                        -.8*.5*(-gf['X'].min()+gf['X'].max())
WAVEMASK             = 1+0*gf['X']
WAVEMASK[gf['X']>xw] = 0.
WAVEMASK[gfl>0]      = 0.

Hs_in          = 3.
Tp_in          = 12.
mwd_in         = -90.
wave_fields0   = {'Hs'  :Hs_in*WAVEMASK,\
                  'Tp'  :Tp_in*WAVEMASK,\
                  'mwd' :mwd_in*WAVEMASK}

params_in   = {}
params_in.update({'Hs_init'   :Hs_in})
params_in.update({'T_init'    :Tp_in})
params_in.update({'mwd_init'  :mwd_in})
params_in.update({'conc_init' :c_in})
params_in.update({'h_init'    :h_in})
params_in.update({'Dmax_init' :D_in})


if 1:
   # change integer parameters:
   params_in.update({'SCATMOD'         : 0})
   params_in.update({'ADV_DIM'         : 2})
   params_in.update({'ADV_OPT'         : 2})
   params_in.update({'DO_CHECK_FINAL'  : 1})
   params_in.update({'DO_CHECK_PROG'   : 1})
   params_in.update({'DO_CHECK_INIT'   : 1})
   params_in.update({'STEADY'          : 1})
   params_in.update({'DO_BREAKING'     : 1})
   params_in.update({'DO_ATTEN'        : 1})

if 1:
   # change real parameters:
   duration_hours = 6.0
   params_in.update({'young'    : 5.49e9})
   params_in.update({'visc_rp'  : 13.0})
   params_in.update({'duration' : duration_hours*60*60})
   params_in.update({'CFL'      : 0.7})

if 1:
   # change real parameters:
   params_in.update({'dumpfreq': 10})
   params_in.update({'itest'   : 25})
   params_in.update({'jtest'   : 5})

if TEST_MESH:
   # inputs on test mesh
   ny       = gf['ny']
   nx       = gf['nx']
   mesh_e   = {}
   if ny==1:
      xx = .5*(gf['X'][:-1]+gf['X'][1:])
      cc = .5*(ice_fields0['icec'][:-1]+ice_fields0['icec'][1:])
      hh = .5*(ice_fields0['iceh'][:-1]+ice_fields0['iceh'][1:])
      dd = .5*(ice_fields0['dfloe'][:-1]+ice_fields0['dfloe'][1:])
   else:
      ninterp  = int(np.ceil(ny/2.))
      xx = .5*(gf['X'][:-1,ninterp]+gf['X'][1:,ninterp])
      cc = .5*(ice_fields0['icec'] [:-1,ninterp]+ice_fields0['icec'] [1:,ninterp])
      hh = .5*(ice_fields0['iceh'] [:-1,ninterp]+ice_fields0['iceh'] [1:,ninterp])
      dd = .5*(ice_fields0['dfloe'][:-1,ninterp]+ice_fields0['dfloe'][1:,ninterp])

   yy             = gf['Y'][:-1,ninterp]
   bb             = 0*xx
   Nfloes         = 0*xx
   Nfloes[dd>0]   = cc[dd>0]/dd[dd>0]**2
   mesh_e.update({'x'      :xx})
   mesh_e.update({'y'      :yy})
   mesh_e.update({'conc'   :cc})
   mesh_e.update({'thick'  :hh})
   mesh_e.update({'broken' :bb})
   mesh_e.update({'Nfloes' :Nfloes})
   # print(mesh_e)
   # sys.exit()

   if 0:
      # test plot of grid inputs:
      fig   = plt.figure()
      ax1   = fig.add_subplot(3,1,1)
      ax2   = fig.add_subplot(3,1,2)
      ax3   = fig.add_subplot(3,1,3)
      ax1.plot(gf['X']/1.e3,ice_fields0['icec' ])
      ax2.plot(gf['X']/1.e3,ice_fields0['iceh' ])
      ax3.plot(gf['X']/1.e3,ice_fields0['dfloe'])
      plt.show(fig)
      sys.exit()
   elif 0:
      # test plot of mesh inputs:
      fig   = plt.figure()
      ax1   = fig.add_subplot(3,1,1)
      ax2   = fig.add_subplot(3,1,2)
      ax3   = fig.add_subplot(3,1,3)
      ax1.plot(xx/1.e3,cc)
      ax2.plot(xx/1.e3,hh)
      ax3.plot(xx/1.e3,dd)
      plt.show(fig)
      sys.exit()
else:
   mesh_e   = None

# call gateway between python and pre-compiled f2py module
out = Rwim.do_run_vSdir(ice_fields=ice_fields0,wave_fields=wave_fields0,\
                                       params_in=params_in,mesh_e=mesh_e)
out_fields,results = out[:2]


##########################################################################
# Make plots
if mesh_e is not None:
   mesh_out = out[2]
   fig      = plt.figure()
   ax       = fig.add_subplot(1,1,1)

   # input on mesh
   lines = ax.plot(xx/1.e3,dd,'--b') # tuple of length 1

   # output on original grid
   lines+= ax.plot(gf['X'][:,0]/1.e3,out_fields['dfloe'][:,ninterp],'k') # tuple of length 2

   # output on mesh
   Nfloes2        = mesh_out['Nfloes']
   dd2            = 0*xx
   dd2[Nfloes2>0] = np.sqrt(cc[Nfloes2>0]/Nfloes2[Nfloes2>0])
   lines+= ax.plot(xx/1.e3,dd2,'r') # tuple of length 3

   ax.legend(lines,['mesh in','grid out','mesh out'])

   if 1:
      # plot figure and exit
      plt.show(fig)
      sys.exit()
   else:
      testdir  = results.figdir+'/test_mesh/'
      if not os.path.exists(results.figdir):
         os.mkdir(results.figdir)
      if not os.path.exists(testdir):
         os.mkdir(testdir)
      testfig1 = testdir+'test_mesh_interp.png'
      print('Saving '+testfig1+'...')
      fig.savefig(testfig1)
##########################################################################

if DO_PLOTTING==0:
   # stop here
   
   w2d   = os.getenv('WIM2D_PATH')
   print('#####################################################################')
   print("To plot init, final, and progress files, type:")
   print(w2d+"/fortran/tools/plot_prog.sh 0 "+results.rootdir)
   print("(0 = no movie; change to 1 to make a movie)\n")

   print("Or, to just plot 1d slices, type:")
   print(w2d+"/fortran/tools/plot_prog_profiles.sh 0 "+results.rootdir)
   print('#####################################################################')

else:
   ##########################################################################
   # Look at initial fields:
   results.plot_initial()
   ##########################################################################

   ################################################################
   # Look at end results:
   results.plot_final()
   ################################################################

   #########################
   # Plot progress files (if they exist)
   # - as colormaps
   results.plot_prog()
   print(' ')
   print("To make a movie of progress images:")
   print(w2d+'/fortran/tools/prog2mp4.sh Hs '+results.figdir+'/prog')
   print('(or try with eg Dmax,taux,tauy)')
   ################################################################

if 0:
   # TODO implement this as part of wim_results class
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
