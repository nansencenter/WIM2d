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
N        = 2**5
h        = 2.   # thickness [m]
period   = 12.  # period [s]
youngs   = 5.e9 # period [s]
visc_rp  = 0.   # R-P damping parameter [m^s/s]
#
Hs       = 3.   # sig wave height [m]
conc     = .7   # concentration [0-1]
dmean    = 100  # mean floe size [m]
#
om          = 2*np.pi/period
gravity     = 9.81
atten_in    = np.array([h,om,youngs,visc_rp])
atten_out   = Mwim.atten_youngs(atten_in)
# print(atten_in)
# print(atten_out)

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

elif 1:
   #finite width [width in metres]:
   #for width in [50.,500.,5000.,50.e3,500.e3,5000.e3]:
   for width in [500.e3]:

      print('\n')
      print('width = '+str(width)+'m')
      print('\n')
      out   = Fbs.solve_boltzmann_ft(width=width,
               alp=alp,N=N,alp_dis=alp_dis,cg=cg,f_inc=Fbs.dirspec_inc_spreading,Hs=Hs)
      out2  = Fbs.solve_boltzmann_ft(width=width,
               alp=alp,N=N,alp_dis=alp_dis,cg=cg,f_inc=Fbs.dirspec_inc_spreading,Hs=Hs)

      if 1:
         # test edge conditions
         semiinf  = False
         Fbs.test_edge_cons(out,semiinf=semiinf,lhs=True)
         Fbs.test_edge_cons(out,semiinf=semiinf,lhs=False)

      # plot energy
      out_plot = Fbs.plot_energy(out,width=width,n_test=0)
