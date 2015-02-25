import numpy as np
import os
import sys
import struct
import matplotlib.pyplot as plt
# import matplotlib.rcsetup as rc

# ##
# ## NB run from 'run' directory !!
# ##
# dd    = os.path.abspath("..")
# dirs  = [dd+"/bin",dd+"/run",dd+"/py_funs"]
# for n in range(0,len(dirs)):
#    dd2   = dirs[n]
#    if not(dd2 in sys.path):
#       print('adding path : '+dd2)
#       sys.path.append(dd2)

import WIM2d_f2py             as Mwim
import run_WIM2d              as Rwim
import fns_get_data           as Fdat
import fns_plot_data          as Fplt
import fns_boltzmann_steady   as Fbs

# linear algebra
from scipy import linalg as LA
#eig   = LA.eig          # 
#solve = np.linalg.solve # x=solve(A,b) -> solves A*x=b

# get grid
gdir        = 'inputs'
grid_prams  = Fdat.fn_check_grid(gdir)
nx          = grid_prams['nx']
ny          = grid_prams['ny']
dx          = grid_prams['dx']
dy          = grid_prams['dy']
X           = grid_prams['X']
Y           = grid_prams['Y']
LANDMASK    = grid_prams['LANDMASK']
xmin        = X.min()
xmax        = X.max()
################################################

##############################################
N        = 2**9
h        = 1.   # thickness [m]
period   = 10.  # period [s]
youngs   = 5.e9 # period [s]
visc_rp  = 0.   # R-P damping parameter [m^s/s]
#
Hs       = 1.   # sig wave height [m]
conc     = .7   # concentration [0-1]
dmean    = 30   # mean floe size [m]
#
om          = np.pi/period
gravity     = 9.81
atten_in    = np.array([h,om,youngs,visc_rp])
atten_out   = Mwim.atten_youngs(atten_in)
alp         = conc/dmean*atten_out[4]  # scattering "attenuation" [m^{-1}]
alp_dis     = 2*conc*atten_out[0]      # damping [m^{-1}]
kwtr        = atten_out[2]
cp          = om/kwtr # phase vel (open water) [m/s]
cg          = cp/2.   # group vel (open water, inf depth relation) [m/s]

if 0:
   #semi-infinite:
   out   = Fbs.solve_boltzmann_ft_semiinf(
            alp=alp,N=N,alp_dis=alp_dis,cg=cg,f_inc=None,Hs=Hs)

   cn       = out['eig_coeffs']
   M_c2ft   = out['edge_lhs']['M_c2ft']
   M_c2th   = out['edge_lhs']['M_c2th']
   semiinf  = True
   width    = None

   # test edge conditions
   Fbs.test_edge_cons(out,semiinf=semiinf,lhs=True)

   # plot energy
   out_plot = Fbs.plot_energy(out,width=None,n_test=0,Hs=Hs)
   xx       = out_plot['x']
   E0       = out_plot['E_n']
   E_coh    = out_plot['E_coh']#TODO

   Hs       = 4*np.sqrt(E0.real)
   Hs_f     = 4*np.sqrt(E_coh)

   ddir  = 'fig_scripts/datfiles'
   if not os.path.exists(ddir):
      os.mkdir(ddir)

   if width is None:
      out_file = ddir+'/steady_semiinf.dat'
   else:
      out_file = ddir+'/steady_L'+str(width)+'.dat'

   print('\n')
   print('Saving data points to '+out_file)
   of1   = open(out_file,'w')
   of1.write('x, m       Hs, m       Hs (coherent), m\n')
   for n in range(len(xx)):
      of1.write('%f   %f    %f\n'%(xx[n],Hs[n],Hs_f[n]))

   of1.close()
elif 1:
   #finite width [width in metres]:
   #for width in [50.,500.,5000.,50.e3,500.e3,5000.e3]:
   for width in [500.e3]:

      print('\n')
      print('width = '+str(width)+'m')
      print('\n')
      out   = Fbs.solve_boltzmann_ft(width=width,
               alp=alp,N=N,alp_dis=alp_dis,cg=cg,f_inc=None,Hs=Hs)

      if 0:
         # test edge conditions
         semiinf  = False
         Fbs.test_edge_cons(out,semiinf=semiinf,lhs=True)
         Fbs.test_edge_cons(out,semiinf=semiinf,lhs=False)

      # plot energy
      out_plot = Fbs.plot_energy(out,width=width,n_test=0)
      xx       = out_plot['x']
      E0       = out_plot['E_n']
      E_coh    = out_plot['E_coh']

      # convert energy to sig wave height
      Hs       = 4*np.sqrt(E0.real)
      Hs_f     = 4*np.sqrt(E_coh)

      ddir  = 'fig_scripts/datfiles'
      if not os.path.exists(ddir):
         os.mkdir(ddir)

      if width is None:
         out_file = ddir+'/steady_semiinf.dat'
      else:
         out_file = ddir+'/steady_L'+str(width)+'.dat'

      print('\n')
      print('Saving data points to '+out_file)
      of1   = open(out_file,'w')
      of1.write('x, m       Hs, m       Hs (coherent), m\n')
      for n in range(len(xx)):
         of1.write('%f   %f    %f\n'%(xx[n],Hs[n],Hs_f[n]))

      of1.close()
