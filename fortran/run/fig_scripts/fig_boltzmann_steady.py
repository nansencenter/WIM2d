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
print('alp_scat,alp_dis (m$^{-1}$): '+str(alp)+','+str(alp_dis))

if 1:
   #semi-infinite:
   out   = Fbs.solve_boltzmann(alp=alp,N=N,alp_dis=alp_dis,Hs=Hs,cg=cg,
                               f_inc=Fbs.dirspec_inc_spreading)

   cn       = out['eig_coeffs']
   M_c2th   = out['edge_lhs']['M_c2th']
   semiinf  = True
   width    = None

   # test edge conditions
   Fbs.test_edge_cons(out,semiinf=semiinf,lhs=True)

   # # plot energy
   # out_plot = Fbs.plot_energy(out,width=None,n_test=0,Hs=Hs)

elif 1:
   #finite width [width in metres]:
   #for width in [50.,500.,5000.,50.e3,500.e3,5000.e3]:
   for width in [150.e3]:

      print('\n')
      print('width = '+str(width)+'m')
      print('\n')

      out   = Fbs.solve_boltzmann(width=width,
                                  alp=alp,N=N,alp_dis=alp_dis,cg=cg,Hs=Hs,
                                  f_inc=Fbs.dirspec_inc_spreading)

      th_vec   = out['angles']
      dth      = 2*np.pi/N
      if 0:
         # check edge conditions
         fig   = plt.figure()
         ax1   = fig.add_subplot(1,2,1)
         ax2   = fig.add_subplot(1,2,2)
         #
         # ax1.plot(out['edge_lhs']['E_edge'])
         ax1.plot(th_vec/np.pi,out['edge_lhs']['E_edge'])
         ax1.plot(th_vec/np.pi,out['edge_lhs']['dirspec_inc'],'--r')
         Hs0   = 4*np.sqrt(dth*np.sum(out['edge_lhs']['E_edge']))
         ax1.set_title('LH edge, $H_s=$'+str(Hs0.real)+'m')
         ax1.set_xlabel(r'$\theta/\pi$')
         ax1.set_ylabel(r'$E$, m$^2$')
         #
         # ax2.plot(out['edge_rhs']['E_edge'])
         ax2.plot(th_vec/np.pi,out['edge_rhs']['E_edge'])
         ax2.plot(th_vec/np.pi,out['edge_rhs']['dirspec_inc'],'--r')
         Hs1   = 4*np.sqrt(dth*np.sum(out['edge_rhs']['E_edge']))
         ax2.set_title('RH edge, $H_s=$'+str(Hs1.real)+'m')
         ax2.set_xlabel(r'$\theta/\pi$')
         ax2.set_ylabel(r'$E$, m$^2$')
         #
         plt.show(fig)
         ax1.cla()
         ax2.cla()
         plt.close(fig)

      if 1:
         xx       = np.linspace(0,width,num=800)
         E_all    = Fbs.calc_expansion(out,xx,L=width)
         En_all   = E_all.dot(out['M_E2En'].transpose())
         #
         fig   = plt.figure()
         ax1   = fig.add_subplot(1,2,1)
         ax1.plot(xx/1.e3,4*np.sqrt(En_all[:,0]))
         ax1.set_xlabel('$x$, km',fontsize=16)
         ax1.set_ylabel('$H_s$, m',fontsize=16)

         if 1:
            # compare to numerical scheme from matlab
            tdir  = '/Users/timill/GITHUB-REPOSITORIES/WIM2d/matlab/boltzmann/out/'
            tfil  = tdir+'boltzmann_steady_numeric.dat'
            #
            t_out = Fdat.read_datfile(tfil)[0]
            xm    = t_out['x'] .data
            Hm    = t_out['Hs'].data
            ax1.plot(xm/1.e3,Hm,'--r')

         # check ODEs
         Lmat  = out['solution']['Lmat']
         Rmat  = out['solution']['Rmat']
         # Lmat2 = out ['solution']['Lmat']
         # Rmat2 = out ['solution']['Rmat']
         #
         dx       = xx[1]-xx[0]
         xmid     = xx[1:-1]
         ax2      = fig.add_subplot(1,2,2)

         if 0:
            # plot integral with cos(n\theta)
            dE_mid   = (E_all[2:,:]-E_all[:-2,:])/(2.*dx)
            LHS      = dE_mid.dot(Lmat.transpose()) # Lmat*[num deriv]
            RHS      =  E_all.dot(Rmat.transpose()) # Rmat*E

            nc = 0 # nc>=2 agrees well
            wt = dth*np.cos(nc*th_vec)
            ax2.plot(xmid/1.e3,LHS.dot(wt))
            ax2.plot(xx/1.e3  ,RHS.dot(wt),'--r')
            ax2.set_title(r'$\int S\cos(n\theta)'+'d'+r'\theta$, n = '+str(nc))

         elif 1:
            # plot one element of E
            dE_mid   = (E_all[2:,:]-E_all[:-2,:])/(2.*dx)
            dE_all   = Fbs.calc_expansion_deriv(out,xx,L=width)
            LHS      = dE_mid.dot(Lmat.transpose()) # Lmat*[num   deriv]
            LHS2     = dE_all.dot(Lmat.transpose()) # Lmat*[exact deriv]
            RHS      =  E_all.dot(Rmat.transpose()) # Rmat*E
            #
            Ntst  = 0
            ax2.plot(xx  /1.e3,RHS[:,Ntst],'b')
            ax2.plot(xmid/1.e3,LHS[:,Ntst],'g')
            ax2.plot(xx  /1.e3,LHS2[:,Ntst],'--r')
            ax2.set_title(r'testing ODE at $\theta=$'+str(180/np.pi*th_vec[Ntst])+'$^o$')

         plt.show(fig)
         ax1.cla()
         ax2.cla()
         plt.close(fig)

    
      # if 1:
      #    # test edge conditions
      #    semiinf  = False
      #    Fbs.test_edge_cons(out,semiinf=semiinf,lhs=True)
      #    Fbs.test_edge_cons(out,semiinf=semiinf,lhs=False)

      # plot energy
      if 0:
         out_plot = Fbs.plot_energy(out,width=width,n_test=0)
