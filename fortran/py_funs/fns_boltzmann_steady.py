import numpy as np
import os
import sys
import struct
import matplotlib.pyplot as plt
# import matplotlib.rcsetup as rc

import WIM2d_f2py    as Mwim
import run_WIM2d     as Rwim
import fns_get_data  as Fdat
import fns_plot_data as Fplt
import fns_inprod    as Finpr

# linear algebra
from scipy import linalg as LA
eig   = LA.eig          #
solve = np.linalg.solve # x=solve(A,b) -> solves A*x=b

##############################################
def dirspec_inc_spreading(th_vec,inputs=None):

   if inputs is None:
      Hs = 1.
   else:
      Hs    = inputs['Hs']
      mwd   = inputs['mwd']
      th0   = -np.pi/180.*(90.+mwd) # from -90 deg -> to 0 rad
                                    # from   0 deg -> to -pi/2 rad

   # incident spectrum:
   # = 2/pi*cos^2(th)
   cc             = np.cos(th_vec-th0)
   D_inc          = 2.0/np.pi*cc**2
   D_inc[cc<0.0]  = 0.0
      # integral=1, so this corresponds to Hs=4*sqrt(1)=4

   D_inc = D_inc*pow(Hs/4.,2)
   # print(th_vec)
   # print(D_inc)
   # sys.exit('dirspec_inc_spreading')

   return D_inc
##############################################

##############################################
def dirspec_inc_plane(th_vec,inputs):

   if inputs is None:
      Hs    = 1.
      dth   = th_vec[1]-th_vec[0]
   else:
      Hs    = inputs['Hs']
      dth   = inputs['dth']

   # incident spectrum:
   # = A^2/2\delta(\theta) : approximate with constant in neigbouring cells to \theta=0
   A  = Hs/2./np.sqrt(2.)
   C  = A*A/2./(2.*dth)
   #
   D_inc                            = 0.*th_vec
   D_inc[abs(th_vec)<dth]           = C
   D_inc[abs(2*np.pi-th_vec)<dth]   = C
   if D_inc.sum()>C:
      # zero only corresponds to one interval
      # => must double C
      D_inc = 2*D_inc

   return D_inc
##############################################

# ##############################################
# def matrix_isotropic(alp,N):
#
#    Mb = alp*np.ones((N,N))/float(N)
#    for n in range(N):
#       Mb[n,n]  = Mb[n,n]-alp
#
#    # angles:
#    th_vec   = np.linspace(0.,2*np.pi,N+1)
#    th_vec   = th_vec[:-1]  # last element is 2\pi \equiv 0
#    C        = np.diag(np.cos(th_vec))
#    # print(C)
#
#    return Mb,th_vec,C
# ##############################################

###############################################
def read_ft_kernel(filename,Nout):
   f     = open(filename)
   lines = f.readlines()
   f.close()

   nl    = len(lines)
   start = 0
   Kcos  = []
   for n in range(nl):
      ll = lines[n].split()
      if len(ll)>0:
         if ll[0]=='Coefficients':
            start = 1
         elif start:
            Kcos.append(float(ll[0]))

   No2      = len(Kcos)-1
   N        = 2*No2
   K_fou     = np.zeros(Nout)
   K_fou[0]  = .5*Kcos[0]
   for n in range(1,min(int(Nout/2.),No2)):
      K_fou[n]       = .5*Kcos[n]
      K_fou[Nout-n]  = .5*Kcos[n]

   # print(Kcos)
   # print(K_fou)
   return K_fou
###############################################

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

   elif OPT is 2:
      print('Warning!! - not isotropic scattering !!')
      print('Loading kernel from file:')
      # load from file
      fdir  = 'fig_scripts/Boltzmann_kernel/out/'
      #dfil  = fdir+'Kfou_h1_E4_T10_D150.dat'
      dfil  = fdir+'Kfou_h1_E5_T10_D150.dat'
      #dfil  = fdir+'Kfou_h1_E5_T5_D150.dat'
      #dfil  = fdir+'Kfou_test.dat'
      print(dfil+'\n')

      # start with isotropic scattering
      # (K=alp, so K_n = \delta_{n,0})
      K_fou = read_ft_kernel(dfil,N)
      alp   = K_fou[0]

   return K_fou,alp
##############################################

##############################################
def solve_boltzmann_ft_semiinf(alp=1.0,N=8,alp_dis=0.0,cg=1.0,Hs=1.,f_inc=None):
   # steady state solution - semi-finite width
   # \pa_t.E + c_g.\pa_x.E=S -> c_g.\pa_x.E=S

   # option for kernel:
   # OPT   = 0 # isotropic scattering
   # OPT   = 1 # add some directionality
   OPT   = 2 # real coefficients from file

   # get fourier coefficients of kernel:
   # K=\sum_n{ K_n/2/pi*exp(-1i*n*theta) }
   K_fou,alp   = get_ft_kernel(alp,N,OPT=OPT)
   sn          = cg*(K_fou-alp-alp_dis) # Sn=sn*En : coeffs of Ft of source function

   ##############################################
   # LH matrix: c_g*cos(theta)*d/dx
   Lmat           = np.zeros((N,N))
   Lmat[0,1]      = .5*c_g
   Lmat[N-1,N-2]   = .5*c_g
   for n in range(1,N-1):
      Lmat[n,[n-1,n+1]] = .5*c_g
   ##############################################

   ##############################################
   # RH matrix: S
   Ds       = np.diag(sn)

   # get eigenvalues & eigenvectors
   # - solve DS*U=\lambda*Lmat*U
   evals,U  = LA.eig(Ds,b=Lmat)
   evals    = evals.real        # inv(Lmat)*Ds = real, symmetric => eig's should be real

   # sort evals
   nn       = np.arange(N)
   jz       = nn[evals==0]  # here are the zero eigenvals
   jp       = nn[evals>0]   # here are the positive eigenvals
   jm       = nn[evals<0]   # here are the negative eigenvals [decay as required]
   jkeep    = np.concatenate([jz[0:1],jm])

   # if alp_dis==0:
   # 2 zero evals, but e-vecs are lin dep so can keep either
   # (other grows linearly)
   # always discard pos evals (grow exponentially)
   # keep neg evals (grow exponentially)
   M_c2ft   = U[:,jkeep]
   eig_vecs = [M_c2ft]
   lam      = -evals[jkeep]
   if 0:
      U0 = M_c2ft[:,0]

      # solve Ds*V0=Lmat*U0
      v0 = Lmat[1:,:].dot(U0)
      V0 = np.concatenate([[0.],v0/sn[1:]])

   ##############################################

   ##############################################
   if 0:
      # test ode solved correctly
      Nc    = len(lam)
      Lmi   = LA.inv(Lmat)
      Mr    = Lmi.dot(Ds)
      npts  = 500
      xx    = np.linspace(0.0,1.0e3,npts)
      f     = np.zeros((N,npts)) # (fourier modes,position)
      g     = np.zeros((N,npts))
      df    = np.zeros((N,npts))

      # eigen-pair to test
      for nc in [1]:
         e_vec = M_c2ft[:,nc]
         e_val = -lam[nc]

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
   # default incident wave spectrum
   # if f_inc is None:
   #    f_inc    = dirspec_inc_spreading  # spreading around central peak - Hs=1 is default
   if f_inc is None:
      f_inc    = dirspec_inc_plane      # delta function - numerical approx

   if f_inc is dirspec_inc_plane:
      finc_in  = {'Hs':Hs,'dth':dth}
   elif f_inc is dirspec_inc_spreading:
      finc_in  = {'Hs':Hs,'mwd':-90}
   ##############################################

   ##############################################
   # find coeffs of eigenfunction expansion
   # from incident spectrum:
   if ECON_SOLN is 0:
      # explicit solution
      # (assumes incident wave is even)
      D_inc    = f_inc(th_vec,finc_in)
      n_inc    = nn[np.cos(th_vec)>0] # N/2 unknowns
      coeffs   = solve(M_c2th[n_inc,:],D_inc[n_inc])
   ##############################################

   ##############################################
   En_edge  = M_c2ft.dot(coeffs)
   E_edge   = M_c2th.dot(coeffs)
   Sn_edge  = Ds.dot(En_edge) # source

   edge_lhs = {'En_edge'      : En_edge,
               'E_edge'       : E_edge,
               'M_c2th'       : M_c2th,
               'M_c2ft'       : M_c2ft,
               'Sn_edge'      : Sn_edge,
               'dirspec_inc'  : D_inc,
               'which_edge'   : 'lhs'}

   out   = {'edge_lhs'     : edge_lhs,
            'eig_coeffs'   : coeffs,
            'eig_vals'     : lam,
            'eig_vecs'     : eig_vecs,
            'M_En2E'       : M_ift,
            'inputs'       : [alp,N,alp_dis,cg,Hs,f_inc],
            'angles'       : th_vec}

   return out
##############################################

##############################################
def solve_boltzmann_ft(width=1.,alp=1.0,N=8,alp_dis=0.0,cg=1.0,Hs = 1.,f_inc=None):
   # steady state solution - finite width
   # \pa_t.E + c_g.\pa_x.E=S -> c_g.\pa_x.E=S

   # option for kernel:
   OPT   = 0 # isotropic scattering
   # OPT   = 1 # add some directionality
   # OPT   = 2 # real coefficients from file

   # print('alp = '+str(alp))
   # print('width = '+str(width))
   # print(np.exp(-alp*width))
   # sys.exit('solve_boltzmann_ft')

   # get fourier coefficients of kernel:
   # K=\sum_n{ K_n/2/pi*exp(-1i*n*theta) }
   K_fou,alp   = get_ft_kernel(alp,N,OPT=OPT)
   sn          = cg*(K_fou-alp-alp_dis) # Sn=sn*En : coeffs of Ft of source function

   ##############################################
   # LH matrix: c_g*cos(theta)*d/dx
   Lmat           = np.zeros((N,N))
   Lmat[0,1]      = .5*cg
   Lmat[N-1,N-2]  = .5*cg
   for n in range(1,N-1):
      Lmat[n,[n-1,n+1]] = .5*cg
   ##############################################

   ##############################################
   # RH matrix: S
   Ds       = np.diag(sn)

   # get eigenvalues & eigenvectors
   # - solve DS*U=\lambda*Lmat*U
   evals,U  = LA.eig(Ds,b=Lmat)
   evals    = evals.real        # inv(Lmat)*Ds = real, symmetric => eig's should be real

   # sort evals
   nn       = np.arange(N)
   jz       = nn[evals==0]  # here are the zero eigenvals
   jp       = nn[evals>0]   # here are the positive eigenvals
   jm       = nn[evals<0]   # here are the positive eigenvals
   jkeep    = np.concatenate([jz[0:1],jm])

   # N/2-1 corresp to lh edge:
   evl   = evals[jm] # <0, scattered by lh edge
   Ul    = U[:,jm] #

   # N/2-1 corresp to rh edge:
   evr   = evals[jp] # >0, scattered by rh edge
   Ur    = U[:,jp] #

   No2   = int(N/2.)
   expL  = np.exp(-evr*width)# = exp(evl*width)

   M_c2ft_0       = np.zeros((N,N))*0j
   M_c2ft_L       = np.zeros((N,N))*0j
   if len(jz)>0:
      """
      zero is a repeated eigenvalue,
      so need to find coefficient vector of 'x':
      - look for y=x*U0+V0
        - U0 is one of the original eigenvectors: DS*U0=\lambda*Lmat*U0
        - V0 satisfies (DS-\lambda*Lmat)*V0=Lmat*U0,
          so since \lambda=0,
          DS*V0=Lmat*U0
      """
      M     = No2-1
      lam   = np.concatenate([[0.],np.array(evr)])
      U0    = U[:,jz[0]]
      v0    = Lmat[1:,:].dot(U0)
      V0    = np.concatenate([[0.],v0/sn[1:]]) # V0[0] is arbitrary as sn[0]=0

      """
      General soln for E_n is:
      y  = sum_n a_n*Ul_n*exp(evl_n*x)
            + sum_n b_n*Ur_n*exp(evr_n*(x-L))
            + c0*U0+c1*(x*U0+V0)
       => 2+2*(N/2-1)=N unknowns: coeffs=[c0,c1,[a_n],[b_n]]
      """
      # E_n, x=0:
      M_c2ft_0[:,0]  = U0
      M_c2ft_0[:,1]  = V0

      # E_n, x=L:
      M_c2ft_L[:,0]  = U0
      M_c2ft_L[:,1]  = V0+U0*width

      eig_vecs = [U0,V0,Ul,Ur]
   else:
      """
      General soln for E_n is:
      y  = sum_n a_n*Ul_n*exp(evl_n*x)
            + sum_n b_n*Ur_n*exp(evr_n*(x-L))
       => 2*(N/2)=N unknowns: coeffs=[[a_n],[b_n]]
      """
      M        = No2
      lam      = np.array(evr)
      eig_vecs = [Ul,Ur]

   # E_n, x=0:
   n0 = N-2*M
   M_c2ft_0[:,n0:n0+M]  = Ul
   M_c2ft_0[:,n0+M:]    = Ur.dot(np.diag(expL))

   # E_n, x=L:
   M_c2ft_L[:,n0:n0+M]  = Ul.dot(np.diag(expL))
   M_c2ft_L[:,n0+M:]    = Ur
   ##############################################

   ##############################################
   if 0:
      # test ode solved correctly
      Nc    = len(lam)
      Lmi   = LA.inv(Lmat)
      Mr    = Lmi.dot(Ds)
      npts  = 500
      xx    = np.linspace(0.0,1.0e3,npts)
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
   M        = int(round(N/2.0))
   dth      = 2*np.pi/N
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
   ##############################################

   ##############################################
   """
   E_n -> E(\theta)
   """
   M_c2th_0 = M_ift.dot(M_c2ft_0)   # x=0
   M_c2th_L = M_ift.dot(M_c2ft_L)   # x=L
   ##############################################

   ##############################################
   # Default incident wave spectrum
   # if f_inc is None:
   #    f_inc    = dirspec_inc_spreading  # spreading around central peak - Hs=1 is default
   if f_inc is None:
      f_inc    = dirspec_inc_plane      # delta function - numerical approx

   if f_inc is dirspec_inc_plane:
      finc_in  = {'Hs':Hs,'dth':dth}

   elif f_inc is dirspec_inc_spreading:
      finc_in  = {'Hs':Hs,'mwd':-90}
   ##############################################

   ##############################################
   # find coeffs of eigenfunction expansion
   # from incident spectrum:
   if ECON_SOLN is 0:
      # explicit solution
      # (assumes incident wave is even)
      rhs      = np.zeros(N)*0j
      Mlhs     = np.zeros((N,N))*0j
      n_incL   = nn[np.cos(th_vec)<0] # N/2 unknowns

      #
      if f_inc is not None:
         print('solve_boltzmann_ft: finc, input')
         # print(f_inc)
         # print(finc_in)
         # sys.exit()
         #
         n_inc0         = nn[np.cos(th_vec)>0] # N/2 unknowns
         D_inc          = f_inc(th_vec,finc_in)
         rhs[0:No2]     = D_inc[n_inc0]      # waves to right only
         Mlhs[0:No2,:]  = M_c2th_0[n_inc0,:] # waves to right only
         Mlhs[No2:,:]   = M_c2th_L[n_incL,:] # waves to left  only
         coeffs         = solve(Mlhs,rhs)
      else:
         # delta function (plane wave) - analytically TODO looks wrong

         # order lhs correctly for inprod with cos over [-pi/2,pi/2]
         ni0      = nn[np.logical_and(th_vec>0,th_vec<np.pi/2.)]
         ni1      = nn[np.logical_and(th_vec<0,th_vec>-np.pi/2.)]
         if len(ni1)==0:
            ni1      = nn[np.logical_and(th_vec>1.5*np.pi,th_vec<2*np.pi)]
         n_inc0   = np.concatenate([ni1,ni0])
         #
         Mip_out        = Finpr.ipmat_cos(No2)
         A              = Hs/2./np.sqrt(2.)

         # Plane wave to right only: A^2/2*\delta(\theta)
         # 2/\pi*\int_0^\pi.cos(n\theta')*A^2/2*\delta(\theta)d\theta'
         # = A^2/\pi*cos(n*\pi/2)
         rhs[0:No2]     = (A*A/np.pi)*np.cos(np.arange(No2)*np.pi/2)
         Mlhs[0:No2,:]  = Mip_out[0].dot(M_c2th_0[n_inc0,:])   # waves to right only
         Mlhs[No2:,:]   = M_c2th_L[n_incL,:]                   # waves to left  only (rhs=0)
         coeffs         = solve(Mlhs,rhs)
         D_inc          = 0*rhs[0:No2]
   ##############################################

   ##############################################
   En_edge_0   = M_c2ft_0.dot(coeffs)
   E_edge_0    = M_c2th_0.dot(coeffs)
   Sn_edge_0   = Ds.dot(En_edge_0) # source
   #
   En_edge_L   = M_c2ft_L.dot(coeffs)
   E_edge_L    = M_c2th_L.dot(coeffs)
   Sn_edge_L   = Ds.dot(En_edge_L) # source

   edge_lhs = {'En_edge'      : En_edge_0,
               'E_edge'       : E_edge_0,
               'Sn_edge'      : Sn_edge_0,
               'M_c2th'       : M_c2th_0,
               'M_c2ft'       : M_c2ft_0,
               'dirspec_inc'  : D_inc,
               'which_edge'   : 'lhs'}

   edge_rhs = {'En_edge'      : En_edge_L,
               'E_edge'       : E_edge_L,
               'Sn_edge'      : Sn_edge_L,
               'M_c2th'       : M_c2th_L,
               'M_c2ft'       : M_c2ft_L,
               'dirspec_inc'  : 0*D_inc,
               'which_edge'   : 'rhs'}

   out      = {'edge_lhs'     : edge_lhs,
               'edge_rhs'     : edge_rhs,
               'M_En2E'       : M_ift,
               'eig_coeffs'   : coeffs,
               'eig_vals'     : lam,
               'eig_vecs'     : eig_vecs,
               'inputs'       : [alp,N,alp_dis,cg,Hs,f_inc],
               'angles'       : th_vec}

   return out
##############################################

##############################################
def test_edge_cons(out,semiinf=True,lhs=True):
   # test LHS edge conditions (semi-infinite)
   print('Test edge conditions:')

   ang   = out['angles']

   if semiinf:
      figname  = 'fig_scripts/figs/SSboltzmann-EdgeCons-semiinf.png'
      edge     = out['edge_lhs']
   elif lhs:
      figname  = 'fig_scripts/figs/SSboltzmann-EdgeCons-lhs.png'
      edge     = out['edge_lhs']
   else:
      figname  = 'fig_scripts/figs/SSboltzmann-EdgeCons-rhs.png'
      edge     = out['edge_rhs']

   if 1:
      print('angles (deg)')
      print(ang[:10]*180/np.pi)
      print('cosines')
      print(np.cos(ang)[:10])
      print('Im(E_edge),x=0')
      print(edge['E_edge'].imag[:10])
      print('Re(E_edge),x=0')
      print(edge['E_edge'].real[:10])
      print('Incident spectrum')
      print(edge['dirspec_inc'][:10])

   if 0:
      plt.plot(ang*180/np.pi,edge['E_edge'].real)
      plt.plot(ang*180/np.pi,edge['dirspec_inc'],'--r')
      plt.show()
   elif 1:
      fig   = Fplt.plot_1d(ang*180/np.pi,edge['E_edge'].imag,
               labs=['Angle,  degrees','$E$'],linestyle='-',color='b')
      Fplt.plot_1d(ang*180/np.pi,edge['E_edge'].real,
               labs=None,f=fig,linestyle='-',color='k')
      Fplt.plot_1d(ang*180/np.pi,edge['dirspec_inc'].real,
               labs=None,f=fig,linestyle='--',color='r')

      plt.savefig(figname,bbox_inches='tight',pad_inches=0.05)
      plt.close()
      fig.clf()
      print(' ')
      print('Saving to file : '+figname)
      print(' ')
   return
##############################################

##############################################
def calc_energy(out,xx,L=None,n_test=None):

   #############################################
   lam   = out['eig_vals']
   cn    = out['eig_coeffs']
   if L is None:
      # semi-infinite
      semiinf  = True
      No2      = len(cn)
      N        = 2*No2
   else:
      # finite width
      semiinf  = False
      N        = len(cn)
      No2      = int(N/2.)

   npts  = len(xx)
   if n_test is None:
      n_test   = range(N)
   E_n   = np.zeros((npts,len(n_test)))*0j

   print('Reconstructing energy from eigenvalues in calc_energy...')

   #############################################
   if not semiinf:
      if len(out['eig_vecs'])==2:
         # only single eigen-values
         print('only simple eigenvalues')
         print('\n')
         #
         M        = No2
         an       = cn[0:M]
         bn       = cn[M:]
      else:
         # repeated eigen-value
         print('zero is a repeated eigenvalue')
         print('\n')
         #
         lam      = lam[1:]
         M        = No2-1
         an       = cn[2:][0:M]
         bn       = cn[2:][M:]
   else:
      print('semi-infinite case')
      print('\n')
   #############################################

   #############################################
   for m in n_test:

      #############################################
      if semiinf:
         ##########################################
         # no need to distinguish between repeated eigenvalue case
         Ul       = out['eig_vecs'][0][m,:]
         E_n[:,m] = 0*xx
         for n in range(npts):
            x        = xx[n]
            cn_x     = np.exp(-lam*x)*cn
            E_n[n,m] = E_n[n,m] + Ul.dot(cn_x)
         ##########################################
      else:
         ##########################################
         # finite width
         if len(out['eig_vecs'])==2:
            # only single eigen-values
            Ul = out['eig_vecs'][0][m,:]
            Ur = out['eig_vecs'][1][m,:]
            #
            E_n[:,m] = 0.*xx
         else:
            # repeated eigen-value
            U0 = out['eig_vecs'][0][m]
            V0 = out['eig_vecs'][1][m]
            Ul = out['eig_vecs'][2][m,:]
            Ur = out['eig_vecs'][3][m,:]
            #
            E_n[:,m] = cn[0]*U0+cn[1]*(U0*xx+V0)

         for n in range(npts):
            x        = xx[n]
            an_x     = np.exp(-lam*x)*an
            bn_x     = np.exp(-lam*(L-x))*bn
            E_n[n,m] = E_n[n,m] + Ul.dot(an_x) \
                                + Ur.dot(bn_x)
      #############################################

   #############################################
   # finished loop over m - return Fourier coefficients
   return E_n
################################################

##############################################
def plot_energy(out,width=None,n_test=0,Hs=1.,f_inc=None):
   # n_test is the Fourier component to plot:
   # (0 is the energy)
   # NB odd n are zero

   fdir  = 'fig_scripts/figs/profiles/'
   fdir2 = 'fig_scripts/figs/profiles/IsoFrac/'
   if not os.path.exists(fdir):
      os.mkdir(fdir)
   if not os.path.exists(fdir2):
      os.mkdir(fdir2)

   if width is None:
      semiinf     = True
      L           = 500.0e3
      figname     = fdir+'SSboltzmann_E'+str(n_test)+'_profile-semiinf.png'
      figname2    = fdir2+'SSboltzmann_IsoFrac_profile-semiinf.png'
      edge_keys   = ['lhs']
   else:
      semiinf     = False
      L           = width
      figname     = fdir+'SSboltzmann_E'+str(n_test)+'_profile-L'+str(width)+'.png'
      figname2    = fdir2+'SSboltzmann_IsoFrac_profile-L'+str(width)+'.png'
      edge_keys   = ['lhs','rhs']

   # plot energy vs x:
   npts  = 5000
   xx    = np.linspace(0.0,L,npts)
   E_n   = 0.0j*xx
   lam   = out['eig_vals']
   #
   cn       = out['eig_coeffs']
   M_c2ft0  = out['edge_lhs']['M_c2ft']
   M_c2th   = out['edge_lhs']['M_c2th']
   alp      = out['inputs'][0]
   alp_dis  = out['inputs'][2]

   ################################################################################
   # calc Fourier coefficients at x:
   if semiinf:
      # semi-infinite
      E_n      = calc_energy(out,xx,L=None,n_test=[n_test])
      No2      = len(cn) # N/2
      N        = 2*No2
      En_si    = E_n

      En_all   = calc_energy(out,xx,L=None)
      E_coh    = 0*xx # coherent energy: E going to right
      iso_frac = 0*xx # amount of isotropy measure
      dth      = out['angles'][1]-out['angles'][0]
      for n in range(npts):
         E_th        = out['M_En2E'].dot(En_all[n,:])
         E_f         = E_th[abs(out['angles'])<dth].real # two closest intervals around 0
         E_coh[n]    = dth*(E_f.sum())                   # integrate over these intervals
         #
         e0          = abs(E_n[n])*abs(E_n[n])
         ea          = abs(En_all[n,:])*abs(En_all[n,:])
         iso_frac[n] = e0/np.sum(ea)
   else:
      # finite width
      E_n   = calc_energy(out,xx,L=L,n_test=[n_test])
      N     = len(cn)
      No2   = int(N/2.)
      if 1:
         comp_si  = 1
         cg       = out['inputs'][3]
         Hs       = out['inputs'][4]
         f_inc    = out['inputs'][5]
         out_si   = solve_boltzmann_ft_semiinf(alp=alp,N=N,alp_dis=alp_dis,cg=cg,Hs=Hs,f_inc=f_inc)
         En_si    = calc_energy(out_si,xx,L=None,n_test=[n_test]) # semi-infinite profile
      #
      En_all   = calc_energy(out,xx,L)
      E_coh    = 0*xx # coherent energy: E going to right
      iso_frac = 0*xx # amount of isotropy measure
      dth      = out['angles'][1]-out['angles'][0]
      for n in range(npts):
         E_th     = out['M_En2E'].dot(En_all[n,:])
         E_f      = E_th[abs(out['angles'])<dth].real # two closest intervals around 0
         E_coh[n] = dth*(E_f.sum())                   # integrate over these intervals
         #
         e0          = abs(E_n[n])*abs(E_n[n])
         ea          = abs(En_all[n,:])*abs(En_all[n,:])
         iso_frac[n] = e0/np.sum(ea)
   ################################################################################

   ang      = out['angles']
   dth      = 2*np.pi/float(N)
   fwd_mask = np.zeros(N)
   fwd_mask[np.cos(ang)>0] = 1.0
   if 1:
      if f_inc is None:
         #dirspec_plane or analytical delta function
         A        = Hs/2./np.sqrt(2.)
         Ec_tst   = A*A/2.
      else:
         #dirspec_spreading
         Ec_tst   = 1./np.pi*(dth+.5*np.sin(2*dth))*pow(Hs/4.,2)

      print('\n')
      print('Check coherent energy (converted to Hs)')
      print('(analytic integration at x=0)')
      print(4*np.sqrt(E_coh[0]),4*np.sqrt(Ec_tst))
      print('\n')
   #
   if 1:
      for key in edge_keys:
         print('************************************************')
         print('check edge E')
         edge     = out['edge_'+key]
         E_edge   = edge['E_edge']
         E0_p     = dth*(fwd_mask*E_edge).sum()
         E0_m     = dth*((1-fwd_mask)*E_edge).sum()

         if n_test==0:
            print('E0 at edge: '+edge['which_edge'])
            print('E0 fwd     = '+str(E0_p))
            print('E0 back    = '+str(E0_m))
            print('E0 total   = '+str(edge['En_edge'][0]))
            print('E0 total 1 = '+str(E0_p+E0_m))

            if edge['which_edge']=='lhs':
               print('E0 total 2 = '+str(E_n[0]))
               print('x          = '+str(xx[0]))
            else:
               print('E0 total 2 = '+str(E_n[-1]))
               print('x          = '+str(xx[-1]))

         print('************************************************')
         print('\n')

   if 0:
      fig   = plt.figure()
      ax    = fig.add_subplot(1,1,1)
      ax.plot(xx,E_n.real)
      ax.set_yscale('log')
      plt.show()
   elif n_test==0:
      fig   = plt.figure()
      ax    = fig.add_subplot(1,1,1)
      if (not semiinf) and (comp_si==1):
         Fplt.plot_1d(xx/1.0e3,4*np.sqrt(En_si.real),f=fig,
               linestyle='-',color='g')
      Fplt.plot_1d(xx/1.0e3,4*np.sqrt(E_n.real),f=fig,
               labs=['$x$, km','$H_s$, m'],linestyle='-',color='k')
      Fplt.plot_1d(xx/1.0e3,4*np.sqrt(E_coh),f=fig,
               linestyle='--',color='r')

      # check vs exponential decay
      Fplt.plot_1d(xx/1.0e3,Hs*np.exp(-alp/2.*xx),f=fig,
               linestyle='--',color='b')

      # check vs equi-partition of energy
      E_eq  = En_si[-1].real/float(N)
      Fplt.plot_1d(xx/1.0e3,4*np.sqrt(E_eq)+0*xx,f=fig,
               linestyle='--',color='g')

      if lam[0]>0.0:
         # dissipation so plot y with log scale (exponential decay)
         ax.set_yscale('log')

      plt.savefig(figname,bbox_inches='tight',pad_inches=0.05)
      plt.close()
      fig.clf()
      print(' ')
      print('Saving to file : '+figname)
      print(' ')

      # plot amount of isotropy:
      fig   = plt.figure()
      Fplt.plot_1d(xx/1.0e3,iso_frac,f=fig,
               labs=['$x$, km','Isotropic fraction'],linestyle='-',color='k')
      Fplt.plot_1d(xx/1.0e3,1-np.exp(-alp*xx),f=fig,
               linestyle='--',color='r')
      plt.savefig(figname2,bbox_inches='tight',pad_inches=0.05)
      plt.close()
      fig.clf()
      print(' ')
      print('Saving to file : '+figname2)
      print(' ')
   elif 1:
      fig   = plt.figure()
      ax    = fig.add_subplot(1,1,1)
      Fplt.plot_1d(xx/1.0e3,E_n.real,f=fig,
               labs=['$x$, km','$E$'],linestyle='-',color='k')

      if alp_dis>0.0:
         ax.set_yscale('log')

      plt.savefig(figname,bbox_inches='tight',pad_inches=0.05)
      plt.close()
      fig.clf()
      print(' ')
      print('Saving to file : '+figname)
      print(' ')

   out   = {'E_n':E_n,'x':xx,'n':n_test,'E_coh':E_coh}

   #######################################################
   # save data to file
   if n_test==0:
      Hs       = 4*np.sqrt(E_n.real)
      Hs_f     = 4*np.sqrt(E_coh)

      ddir  = 'fig_scripts/figs/profiles/datfiles'
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
   #######################################################

   return out
##############################################
