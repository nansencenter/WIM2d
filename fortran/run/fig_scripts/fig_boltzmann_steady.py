import numpy as np
import os
import sys
import struct
import matplotlib.pyplot as plt
# import matplotlib.rcsetup as rc

##
## NB run from 'run' directory !!
##
dd    = os.path.abspath("..")
dirs  = [dd+"/bin",dd+"/run",dd+"/py_funs"]
for n in range(0,len(dirs)):
   dd2   = dirs[n]
   if not(dd2 in sys.path):
      print('adding path : '+dd2)
      sys.path.append(dd2)

import WIM2d_f2py    as Mwim
import run_WIM2d     as Rwim
import fns_get_data  as Fdat
import fns_plot_data as Fplt

# linear algebra
from scipy import linalg as LA
eig   = LA.eig          # 
solve = np.linalg.solve # x=solve(A,b) -> solves A*x=b

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

period   = 10

##############################################
def diag(dd):
   # diagonal matrix with dd as diagonal

   N  = len(dd)

   if abs(dd.imag).sum()!=0.0:
      # if any of dd complex,
      # need to initialise DD as complex
      # (otherwise can't put elements of dd into it)
      DD = 0j*np.ones((N,N))
   else:
      # initialize DD as real
      DD = np.zeros((N,N))

   for n in range(N):
      DD[n,n]  = dd[n]

   return DD
##############################################

##############################################
def linspace(N,a=0.,b=1.):
   v  = np.arange(N)/float(N-1)
   y  = a+v*(b-a)
   return y
##############################################

##############################################
def dirspec_inc(th_vec):

   # incident spectrum:
   # - 2/pi*cos^2(th)
   cc             = np.cos(th_vec)
   D_inc          = 2.0/np.pi*cc**2
   D_inc[cc<0.0]  = 0.0

   return D_inc
##############################################

##############################################
def matrix_isotropic(alp,N):

   Mb = alp*np.ones((N,N))/float(N)
   for n in range(N):
      Mb[n,n]  = Mb[n,n]-alp

   # angles:
   th_vec   = linspace(N+1,0.,2*np.pi)
   th_vec   = th_vec[:-1]  # last element is 2\pi \equiv 0
   C        = diag(np.cos(th_vec))
   # print(C)

   return Mb,th_vec,C
##############################################

##############################################
def get_ft_kernel(alp,N,OPT=0):
   # fourier coefficients of kernel:
   # K=\sum_n{ K_n/2/pi*exp(-1i*n*theta) }
   if OPT is 0:
      # isotropic scattering
      # (K=alp, so K_n = \delta_{n,0})
      K_fou = alp*np.eye(N)[:,0]
   elif OPT is 1:
      # start with isotropic scattering
      # (K=alp, so K_n = \delta_{n,0})
      K_fou = alp*np.eye(N)[:,0]

      # introduce some directionality to scattering
      # K->K*(1+a*cos(\theta))
      a  = .7

      print('Warning!! - not isotropic scattering !!')
      print('Relative amp of cosine = '+str(a)+'\n')

      K_fou[1]    = a/2.0*K_fou[0]
      K_fou[N-1]  = a/2.0*K_fou[0]

   return K_fou
##############################################

##############################################
def solve_boltzmann_ft_semiinf(alp=1.0,N=8,alp_dis=0.0,cg=1.0,f_inc=dirspec_inc):
   # fourier coefficients of kernel:
   # K=\sum_n{ K_n/2/pi*exp(-1i*n*theta) }
   OPT   = 0 # isotropic scattering
   # OPT   = 0 # add some directionality

   K_fou = get_ft_kernel(alp,N,OPT=OPT)

   sn    = cg*(K_fou-alp-alp_dis) # Sn=sn*En : coeffs of Ft of source function

   ##############################################
   # LH matrix
   Lmat           = np.zeros((N,N))
   Lmat[0,1]      = .5
   Lmat[N-1,N-2]   = .5
   for n in range(1,N-1):
      Lmat[n,[n-1,n+1]] = .5
   ##############################################

   ##############################################
   # get eigenvalues & eigenvectors
   Ds       = diag(sn)
   evals,U  = LA.eig(Ds,b=Lmat) # DS*U=\lambda*Lmat*U
   evals    = evals.real        # inv(Lmat)*Ds = real, symmetric 

   # sort evals
   nn       = np.arange(N)
   jz       = nn[evals==0]  # here are the zero eigenvals
   jp       = nn[evals>0]   # here are the positive eigenvals
   jm       = nn[evals<0]   # here are the positive eigenvals
   jkeep    = np.concatenate([jz[0:1],jm])
   
   # if alp_dis==0:
   # 2 zero evals, but e-vecs are lin dep so can keep either
   # (other grows linearly)
   # always discard pos evals (grow exponentially)
   # keep neg evals (grow exponentially)
   M_c2ft   = U[:,jkeep]
   lam      = evals[jkeep]
   ##############################################

   ##############################################
   if 0:
      # test ode solved correctly
      Nc    = len(lam)
      Lmi   = LA.inv(Lmat)
      Mr    = Lmi.dot(Ds)
      npts  = 500
      xx    = linspace(npts,0.0,1.0e3)
      f     = np.zeros((N,npts)) # (fourier modes,position)
      g     = np.zeros((N,npts))
      df    = np.zeros((N,npts))

      # eigen-pair to test
      for nc in [1]:
         e_vec = M_c2ft[:,nc]
         e_val = lam[nc]

         for rx in range(npts):
            x        = xx[rx]
            ff       = e_vec*np.exp(e_val*x)
            f[:,rx]  = ff
            g[:,rx]  = Mr.dot(ff)
            df[:,rx] = e_val*ff


         # fourier mode to test
         r_test   = 3
         print(nc,e_val)
         print(r_test,f[r_test,[0,-1]])
         print(g[r_test,[0,-1]],df[r_test,[0,-1]])

         plt.plot(xx,g[r_test,:].real)
         plt.plot(xx,df[r_test,:].real,'--r')
         plt.show()
   ##############################################

   ##############################################
   # angles:
   M           = int(round(N/2.0))
   dth         = 2*np.pi/N
   shift    = dth/2.0
   th_vec   = dth*nn+shift
   # solution of edge conditions
   # 0: explicit solution
   #    (need to avoid \pm\pi/2)
   ECON_SOLN   = 0

   # indices: go between real index (n>=0)
   # and physical index (n>=-M,n<M)
   n2m   = np.zeros(N)
   m2n   = {0:0}
   for m in range(-M,M):
      if m>=0:
         n  = m
      else:
         n  = N+m
      m2n.update({m:n})
      n2m[n]   = m
   ##############################################

   ##############################################
   # matrices to  go between eigenspace,
   # Fourier space and theta space
   M_ift = 0j*np.zeros((N,N))
   M_ft  = 0j*np.zeros((N,N))# for testing:
   fac   = .5/np.pi
   for n in range(N):
      m           = n2m[n]
      M_ift[:,n]  = fac*np.exp(-1j*m*th_vec)
      M_ft[n,:]   = dth*np.exp(1j*m*th_vec)# for testing:

   M_c2th   = M_ift.dot(M_c2ft)
   ##############################################

   ##############################################
   if 0:
      print(100*'*')
      print('test M_ft/M_ift:')
      m_test   = -1
      print('test n = '+str(m_test))
      f        = fac*np.exp(-1j*m_test*th_vec)
      if 0:
         # test indices mainly
         n  = m2n[m_test]
         g  = fac*np.exp(-1j*n2m[n]*th_vec)
         print('f')
         print(f)
         print('g')
         print(g)
         print('|f-g| = '+str(abs(f-g).sum()))

      if 0:
         # test M_ft
         fn = M_ft.dot(f)
         print('Im[fn]:')
         print(fn.imag)
         print('Re[fn]:')
         print(fn.real)
         n  = m2n[m_test]
         gn = np.eye(N)[:,n] 
         print('Expected fn:')
         print(gn)
         print('|fn-gn| = '+str(abs(fn-gn).sum()))

      if 0:
         print('f')
         print(f)
         n  = m2n[m_test]
         fn = np.eye(N)[:,n]
         g  = M_ift.dot(fn)
         print('g')
         print(g)
         print('|f-g|')
         print('|f-g| = '+str(abs(f-g).sum()))

      print(100*'*'+'\n')
   ##############################################

   ##############################################
   # find coeffs of eigenfunction expansion
   # from incident spectrum:
   if ECON_SOLN is 0:
      # explicit solution
      # (assumes incident wave is even)
      D_inc    = f_inc(th_vec)
      n_inc    = nn[np.cos(th_vec)>0] # N/2 unknowns
      coeffs   = solve(M_c2th[n_inc,:],D_inc[n_inc])
   ##############################################
   
   ##############################################
   En_edge  = M_c2ft.dot(coeffs)
   E_edge   = M_c2th.dot(coeffs)
   Sn_edge  = Ds.dot(En_edge) # source

   out   = {'En_edge': En_edge,'E_edge': E_edge,'Sn_edge': Sn_edge,
            'eig_coeffs':coeffs,'eig_vals':lam,
            'M_c2th':M_c2th,'M_c2ft':M_c2ft,
            'angles':th_vec,'dirspec_inc':D_inc}

   return out
##############################################

N        = 2**8
alp      = 1.0e-5
alp_dis  = 0.0e-5
cg       = 1.0
out      = solve_boltzmann_ft_semiinf(
            alp,N,alp_dis=alp_dis,cg=cg,f_inc=dirspec_inc)

cn       = out['eig_coeffs']
M_c2ft   = out['M_c2ft']
M_c2th   = out['M_c2th']

##############################################
if 1:
   # test edge conditions
   print('Test edge conditions:')

   ang   = out['angles'] 

   if 0:
      print('angles (deg)')
      print(ang*180/np.pi)
      print('cosines')
      print(np.cos(ang))
      print('Im(E_edge)')
      print(out['E_edge'].imag)
      print('Re(E_edge)')
      print(out['E_edge'].real)
      print('Incident spectrum')
      print(out['dirspec_inc'])

   if 0:
      plt.plot(ang*180/np.pi,out['E_edge'].real)
      plt.plot(ang*180/np.pi,out['dirspec_inc'],'--r')
      plt.show()
   elif 1:
      fig   = Fplt.plot_1d(ang*180/np.pi,out['E_edge'].real,
               labs=['Angle,  degrees','$E$'],linestyle='-',color='k')
      Fplt.plot_1d(ang*180/np.pi,out['dirspec_inc'].real,
               labs=None,f=fig,linestyle='--',color='r')

      fname = 'fig_scripts/figs/SSboltzmann-EdgeCons.png'
      plt.savefig(fname,bbox_inches='tight',pad_inches=0.05)
      plt.close()
      fig.clf()
      print('Saving to file : '+fname)
##############################################

##############################################
if 1:
   # plot energy vs x:
   npts  = 5000
   xx    = linspace(npts,0.0,500.0e3)
   E_n   = 0.0j*xx
   lam   = out['eig_vals']

   # Fourier component to plot:
   # (0 is the energy)
   # NB odd n are zero
   n_test   = 0

   for n in range(npts):
      x        = xx[n]
      D_exp    = diag(np.exp(lam*x))
      cn_x     = D_exp.dot(cn)
      E_n[n]   = M_c2ft[n_test,:].dot(cn_x)

   ang      = out['angles'] 
   dth      = 2*np.pi/float(N)
   fwd_mask = np.zeros(N)
   fwd_mask[np.cos(ang)>0] = 1.0
   #
   if 0:
      print('check edge E')
      E_edge   = out['E_edge']
      E0_p     = dth*(fwd_mask*E_edge).sum()
      E0_m     = dth*((1-fwd_mask)*E_edge).sum()
      if n_test==0:
         print('E0 edge')
         print(out['En_edge'][0],E_n[0],E0_p+E0_m) # consistency of eig expansion and fourier coeffs
         print('E0 fwd  = '+str(E0_p))
         print('E0 back = '+str(E0_m))
   else:
      xt = xx[200]
      print('check E, x='+str(xt))
      D_exp    = diag(np.exp(lam*x))
      Mc2th_x  = M_c2th.dot(D_exp)
      Ex       = Mc2th_x.dot(cn)
      Ex_p     = dth*(fwd_mask*Ex).sum()
      Ex_m     = dth*((1-fwd_mask)*Ex).sum()
      print('E fwd  = '+str(Ex_p))
      print('E back = '+str(Ex_m))

      Mc2ft_x  = M_c2ft.dot(D_exp)
      En_x     = Mc2ft_x.dot(cn)
      print('E tot = '+str(En_x[0])+str(Ex_p+Ex_m))

   if 0:
      fig   = plt.figure()
      ax    = fig.add_subplot(1,1,1)
      ax.plot(xx,E_n.real)
      ax.set_yscale('log')
      plt.show()
   elif 1:
      fig   = plt.figure()
      ax    = fig.add_subplot(1,1,1)
      Fplt.plot_1d(xx/1.0e3,E_n.real,f=fig,
               labs=['$x$, km','$E$'],linestyle='-',color='k')

      if alp_dis>0.0:
         ax.set_yscale('log')

      fname = 'fig_scripts/figs/SSboltzmann_E'+str(n_test)+'_profile.png'
      plt.savefig(fname,bbox_inches='tight',pad_inches=0.05)
      plt.close()
      fig.clf()
      print('Saving to file : '+fname)
##############################################
