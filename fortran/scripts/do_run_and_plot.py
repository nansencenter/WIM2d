import numpy as np
import os,sys
import shutil
import struct
from matplotlib import pyplot as plt

sys.path.append("bin")

import run_WIM2d     as Rwim
import fns_get_data  as Fdat
import fns_plot_data as Fplt

# RUN_OPT     = 0 # non-IO version
RUN_OPT     = 1 # rerun then plot
DO_PLOTTING = 0 # cancel plotting if desired by setting to 0

gf          = Fdat.fn_check_grid('grid')
gfl         = gf['LANDMASK']
ICEMASK     = 1.-gfl
WAVEMASK    = 1.-gfl
grid_prams  = gf

###########################################################################
if RUN_OPT==0:
   results  = Rwim.do_run(RUN_OPT=RUN_OPT)
else:
   # set inputs: (icec,iceh,dfloe), (Hs,Tp,mwd)

   xmin  = gf['X'].min()
   xmax  = gf['X'].max()
   xav   = .5*(xmin+xmax)
   xm    = .5*(gf['dx']+xmax-xmin)
   if 1:
      # "semi-infinite" ice sheet

      # ice edge
      xe       = xav -.7*.5*(xmax-xmin)
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
      x0          = xmin-.5*gf['dx']
      xe          = xav -.7*xm
      ICEMASK     = 0*gf['X']
      print(1e-3*np.array([xmin,xmax,x0,xe,xav]))
      #
      critter           = np.logical_and(abs(gf['X'])>xe,\
                                         abs(gf['X'])<xe+strip_width)
      ICEMASK[critter]  = 1.
      ICEMASK[gfl>0]    = 0. # 0 on land
      #
      c_in  = .7
      h_in  = 1.
      D_in  = 100.

   # edge of wave mask
   xw                   = xav -.8*xm
   WAVEMASK             = 1+0*gf['X']
   WAVEMASK[gf['X']>xw] = 0.
   WAVEMASK[gfl>0]      = 0.

   Hs_in       = 3.
   Tp_in       = 12.
   mwd_in      = -90.
   in_fields   = {'icec':c_in*ICEMASK,'iceh':h_in*ICEMASK,'dfloe':D_in*ICEMASK,
                  'Hs':Hs_in*WAVEMASK,'Tp':Tp_in*WAVEMASK,'mwd':mwd_in*WAVEMASK}

   # print(WAVEMASK.max())
   # print(WAVEMASK.sum())
   # sys.exit()
   print(in_fields['Tp'].max())
   print(in_fields['Tp'].min())

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
      params_in.update({'young'   : 5.49e9})
      params_in.update({'visc_rp' : 13.0})
      params_in.update({'duration': duration_hours*60*60})
      params_in.update({'CFL'     : .7})

   if 1:
      # change real parameters:
      params_in.update({'dumpfreq': 10})
      params_in.update({'itest'   : 25})
      params_in.update({'jtest'   : 10})

   # call gateway between python and pre-compiled f2py module
   out_fields,results   = Rwim.do_run(RUN_OPT=RUN_OPT,in_fields=in_fields,
                                       params_in=params_in)

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
   results.plot_initial()
   results.plot_final()

   # clear old progress plots
   pfdir = results.figdir+'/prog'
   if os.path.exists(pfdir):
      old_dirs = os.listdir(pfdir)
      for od in old_dirs:
         # need this for macs
         if od!='.DS_Store':
            # os.rmdir(pfdir+'/'+od)
            shutil.rmtree(pfdir+'/'+od)

   # Make plots
   results.plot_prog()
   
   if 0:
      ################################################################
      # Plot progress files (if they exist)
      # - as profiles on one graph
      cols     = ['k','b','r','g','m','c']
      lstil    = ['-','--','-.',':']
      Nc       = len(cols)
      Ns       = len(lstil)
      loop_c   = -1
      loop_s   = -1
   
      pdir        = bindir+'/prog/'
      prog_files  = os.listdir(pdir)
   
      ################################################################
      # checks
      if len(prog_files)>0:
         # make dir for progress plots
         if (not os.path.exists(figdir)):
            os.mkdir(figdir)
      else:
         raise ValueError('No progress files to plot')
      ################################################################
   
      fig   = plt.figure()
      ax    = fig.add_subplot(1,1,1)
      xx    = gf['X'][:,0]/1.e3
      labs  = ['$x$, km','$H_s$, m']
   
      ################################################################
      for pf in prog_files:
         if '.a'==pf[-2:]:
            print("Reading "+pf+" ...")
            out_fields,info   = Fdat.fn_read_general_binary(pdir+pf)
   
            Hs_n     = out_fields['Hs'][:,0]
            loop_c   = np.mod(loop_c+1,Nc)
            if loop_c==0:
               loop_s   = np.mod(loop_s+1,Ns)
   
            fig,ax,line = Fplt.plot_1d(xx,Hs_n,labs=labs,pobj=[fig,ax],color=cols[loop_c],linestyle=lstil[loop_s])
      #
      figname  = figdir+'/convergence2steady.png'
      print('saving to '+figname+'...')
      fig.savefig(figname,bbox_inches='tight',pad_inches=0.05)
      plt.close(fig)
