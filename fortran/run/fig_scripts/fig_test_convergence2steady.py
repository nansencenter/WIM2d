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
import fns_boltzmann_steady as Fbs

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

   # ice band of finite width
   xe                   = 50.e3
   ICEMASK              = 1+0*gf['X']
   ICEMASK[gf['X']<xe]  = 0.
   ICEMASK[gfl>0]       = 0.

   # right hand ice edge
   ice_width                        = 100.e3
   ICEMASK[gf['X']>(xe+ice_width)]  = 0.

   # edge of wave mask
   xw                   = 30.e3
   WAVEMASK             = 1+0*gf['X']
   WAVEMASK[gf['X']>xw] = 0.
   WAVEMASK[gfl>0]      = 0.

   ########################################################################
   c_in     = .7     # conc
   h_in     = 2.     # thickness
   D_in     = 100.   # initial floe size
   Hs_in    = 3.     # initial Hs
   Tp_in    = 12.    # initial Tp
   mwd_in   = -90    # initial mean wave direction

   in_fields   = {'icec':c_in*ICEMASK,'iceh':h_in*ICEMASK,'dfloe':D_in*ICEMASK,
                  'Hs':Hs_in*WAVEMASK,'Tp':Tp_in*WAVEMASK,'mwd':mwd_in*WAVEMASK}
   ########################################################################

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
   DO_BREAKING = 0 # no breaking - testing convergence of Hs to steady-state
   STEADY      = 1
   int_prams   = np.array([SCATMOD,ADV_DIM,ADV_OPT,
                           CHECK_FINAL,CHECK_PROG,CHECK_INIT,
                           DO_BREAKING,STEADY])

if 1:
   # change real parameters:
   young          = 5.49e9
   visc_rp        = 0.0
   duration_hours = 72.0
   duration       = duration_hours*60*60
   real_prams     = np.array([young,visc_rp,duration])


out_fields,outdir = Rwim.do_run(RUN_OPT=RUN_OPT,in_fields=in_fields,
                                    int_prams=int_prams,
                                    real_prams=real_prams)

##########################################################################
# Make plots
bindir   = outdir+'/binaries'
figdir0  = 'fig_scripts/figs/'
figdir   = figdir0+'TC2S'

if not os.path.exists(figdir0):
   os.mkdir(figdir0)
if not os.path.exists(figdir):
   os.mkdir(figdir)

##########################################################################
grid_prams  = Fdat.fn_check_grid(bindir) # load grid from binaries
xx          = 1.e-3*grid_prams['X'][:,0]
labs        = ['x, km','$H_s$, m']

##########################################################################
if 1:
   # get steady state solution # TODO get this working
   om          = 2*np.pi/Tp_in
   # gravity     = 9.81
   atten_in    = np.array([h_in,om,young,visc_rp])
   atten_out   = Mwim.atten_youngs(atten_in)
   alp         = c_in/D_in*atten_out[4]   # scattering "attenuation" [m^{-1}]
   alp_dis     = 2*c_in*atten_out[0]      # damping [m^{-1}]
   kwtr        = atten_out[2]
   cp          = om/kwtr # phase vel (open water) [m/s]
   cg          = cp/2.   # group vel (open water, inf depth relation) [m/s]
   #
   N     = pow(2,7);
   nx0         = 1.e3
   x_ice       = np.linspace(0.,ice_width,nx0)
   if 1:
      # solve in Fourier space
      meth  = '_FT'
      out   = Fbs.solve_boltzmann_ft(width=ice_width,
               alp=alp,N=N,alp_dis=alp_dis,cg=cg,f_inc=Fbs.dirspec_inc_spreading,Hs=Hs_in)
      #
      E_n         = Fbs.calc_energy(out,x_ice,L=ice_width,n_test=[0])
      Hs_steady   = 4*np.sqrt(E_n[:,0].real)
   else:
      # solve in position space
      meth  = ''
      out   = Fbs.solve_boltzmann(width=ice_width,
               alp=alp,N=N,alp_dis=alp_dis,cg=cg,f_inc=Fbs.dirspec_inc_spreading,Hs=Hs_in)
      #
      E_th        = Fbs.calc_energy(out,x_ice,L=ice_width,n_test=range(N))
      dtheta      = 2*np.pi/float(N)
      E0          = dtheta*E_th.sum(1) # integrate over directions
      Hs_steady   = 4*np.sqrt(E0.real)

   # alp2  = out['inputs'][0]
   # print(alp,alp2)
   # sys.exit('fig_test_convergence2steady')


   # print(x_ice)
   # print(Hs_steady)
   fig      = Fplt.plot_1d(1.e-3*(xe+x_ice),Hs_steady,labs=labs,color='y',linewidth=3)
else:
   fig   = None
##########################################################################

pfiles   = os.listdir(bindir+'/prog')
pf       = pfiles[-1]
Nprog    = int(pf[8:-2])

if 0:
   # compare every stp steps
   stp         = 40
   steps2plot  = range(0,Nprog,stp)
else:
   # just compare final prog file
   steps2plot  = [Nprog]

cols     = ['k','b','r','g','m','c']
lstil    = ['-','--','-.',':']
Nc       = len(cols)
loop_c   = -1
loop_s   = 0

for nstep in steps2plot:
   out_fields  = Fdat.fn_check_prog(outdir,nstep) # load ice/wave conditions from binaries
   Hs_n        = out_fields['Hs'][:,0]
   #
   if loop_c==Nc-1:
      loop_c   = 0
      loop_s   = loop_s+1
   else:
      loop_c   = loop_c+1

   fig      = Fplt.plot_1d(xx,Hs_n,labs=labs,f=fig,color=cols[loop_c],linestyle=lstil[loop_s])
#
figname  = figdir+'/convergence2steady'+meth+'.png'
print('saving to '+figname+'...')
plt.savefig(figname,bbox_inches='tight',pad_inches=0.05)
plt.close()
fig.clf()

if 1:
   #####################################################
   # print  time-dep results to dat-file
   dfil  = figdir+'/test_steady1.dat'
   blk   = 4*' '
   fid   = open(dfil,'w')

   # header:
   lin   = []
   lin.append('# Time-dependant results (end of simulation)\n')
   lin.append('# Time (h): '+str(duration_hours)+'\n')
   lin.append('# x (m), Hs (m)'+'\n')
   lin.append('#####################################################\n')
   lin.append('\n')

   for loop_h in range(len(lin)):
      fid.write(lin[loop_h])
   
   # results:
   for loop_x in range(len(xx)):
      lin   =           ('%f' % (1.e3*xx[loop_x])  )
      lin   = lin+blk + ('%f' % (Hs_n[loop_x])     )
      fid.write(lin+'\n')
   fid.close()
   print('saving to '+dfil+'...')
   #####################################################

   #####################################################
   # print  steady-state results to dat-file
   dfil  = figdir+'/test_steady2'+meth+'.dat'
   fid   = open(dfil,'w')

   # header:
   lin   = []
   lin.append('# Steady-state results\n')
   lin.append('# x (m), Hs (m)'+'\n')
   lin.append('#####################################################\n')
   lin.append('\n')

   for loop_h in range(len(lin)):
      fid.write(lin[loop_h])

   # results:
   for loop_x in range(len(x_ice)):
      lin   =           ('%f' % (xe+x_ice[loop_x])  )
      lin   = lin+blk + ('%f' % (Hs_steady[loop_x])     )
      fid.write(lin+'\n')
   fid.close()
   print('saving to '+dfil+'...')
   #####################################################
