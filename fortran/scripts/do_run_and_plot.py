import numpy as np
import os,sys
import shutil
import struct
from matplotlib import pyplot as plt

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

RUN_OPT     = 2 # rerun then plot
# RUN_OPT     = 3 # plot saved results
DO_PLOTTING = 1 # cancel plotting if desired by setting to 0

gf          = Fdat.fn_check_grid('inputs')
gfl         = gf['LANDMASK']
ICEMASK     = 1.-gfl
WAVEMASK    = 1.-gfl
grid_prams  = gf

###########################################################################
# set inputs: (icec,iceh,dfloe), (Hs,Tp,mwd)
if 1:

   if 0:
      # "semi-infinite" ice sheet

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
      xmin        = gf['X'].min()
      xmax        = gf['X'].max()
      xav         = .5*(xmin+xmax)
      x0          = xmin-.5*gf['dx']
      xm          = .5*(gf['dx']+xmax-xmin)
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
      h_in  = 2.
      D_in  = 100.

   # edge of wave mask
   xw                   = xav -.8*xm
   WAVEMASK             = 1+0*gf['X']
   WAVEMASK[gf['X']>xw] = 0.
   WAVEMASK[gfl>0]      = 0.

   Hs_in       = 2.
   Tp_in       = 12.
   mwd_in      = -90.
   in_fields   = {'icec':c_in*ICEMASK,'iceh':h_in*ICEMASK,'dfloe':D_in*ICEMASK,
                  'Hs':Hs_in*WAVEMASK,'Tp':Tp_in*WAVEMASK,'mwd':mwd_in*WAVEMASK}

int_prams   = None # default integer parameters
real_prams  = None # default real parameters

if 1:
   # change integer parameters:
   SCATMOD     = 1
   ADV_DIM     = 2
   ADV_OPT     = 2
   CHECK_FINAL = 1
   CHECK_PROG  = 1
   CHECK_INIT  = 1
   STEADY      = 1
   DO_BREAKING = 0
   DO_ATTEN    = 1
   int_prams   = np.array([SCATMOD,ADV_DIM,ADV_OPT,
                           CHECK_FINAL,CHECK_PROG,CHECK_INIT,
                           STEADY,DO_BREAKING,DO_ATTEN])

if 1:
   # change real parameters:
   young          = 5.49e9
   visc_rp        = 0.0
   duration_hours = 24.0
   duration       = duration_hours*60*60
   CFL            = .7
   real_prams     = np.array([young,visc_rp,duration,CFL])


# call gateway between python and pre-compiled f2py module
out_fields,outdir = Rwim.do_run(RUN_OPT=RUN_OPT,in_fields=in_fields,
                                    int_prams=int_prams,
                                    real_prams=real_prams)

##########################################################################
if DO_PLOTTING==0:
   # stop here
   
   print('#####################################################################')
   print("To plot init, final, and progress files, type:")
   print(w2d+"/fortran/tools/plot_prog.sh 0 out_io")
   print("(0 = no movie; change to 1 to make a movie)\n")

   print("Or, to just plot 1d slices, type:")
   print(w2d+"/fortran/tools/plot_prog_profiles.sh 0 out_io")
   print('#####################################################################')

else:
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
      figdir3     = figdir+'/prog/'
      pdir        = bindir+'/prog/'
      prog_files  = os.listdir(pdir)
   
      ################################################################
      # checks
      if len(prog_files)>0:
         # make dir for progress plots
         if (not os.path.exists(figdir3)):
            os.mkdir(figdir3)
         else:
            # clear old progress plots
            old_dirs = os.listdir(figdir3)
            for od in old_dirs:
               # need this for macs
               if od!='.DS_Store':
                  # os.rmdir(figdir3+'/'+od)
                  shutil.rmtree(figdir3+'/'+od)
      else:
         raise ValueError('No progress files to plot')
      ################################################################
   
      # Typical limits for parameters
      # NB use the same key names as in the .b file
      zlims  = {'icec':[0,1],      \
                'iceh':[0,5],      \
                'Dmax':[0,300],    \
                'tau_x':[-.5,.5],  \
                'tau_y':[-.05,.05],\
                'Hs':[0,4],        \
                'Tp':[10,20],      \
                'mwd':[-180,180]}
   
      for pf in prog_files:
         if '.a'==pf[-2:]:
            print("Reading "+pf+" ...")
            prog_fields,info  = Fdat.fn_read_general_binary(pdir+pf)
            #
            stepno      = pf.strip('wim_prog').strip('.a')
            figdir3_0   = figdir3+'/'+stepno
            Fplt.fn_plot_gen(grid_prams,prog_fields,figdir3_0,zlims_in=zlims)
            print("Plots in "+figdir3_0+'\n')
      ################################################################
   
      print(' ')
      print("To make a movie of progress images:")
      print(dd+'/tools/prog2mp4.sh Hs '+figdir3)
   
   elif 1:
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
