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
dirs  = [dd+"/bin",dd+"/run",dd+"/misc_py"]
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
   DD = np.zeros((N,N))

   if abs(dd).imag.sum()!=0.0:
      #test if any complex
      DD = 0j*np.ones((N,N))

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
def solve_isotropic_ft(alp=1.0,N=8,alp_dis=0.0,cg=1.0,f_inc=dirspec_inc):
   # fourier coefficients of kernel:
   # K=\sum_n{ K_n/2/pi*exp(-1i*n*theta) }
   K_fou = alp*np.eye(N)[:,0]
   if 1:
      # introduce some directionality
      # K->K*(1+a*cos(\theta))
      a           = .7
      K_fou[1]    = a/2.0*K_fou[0]
      K_fou[N-1]  = a/2.0*K_fou[0]

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
   # get evals
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
   # 2 zero evals, but e-vecs are lin dep so keep 1
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
   M        = int(round(N/2.0))
   dth      = 2*np.pi/N
   shift    = dth/2.0
   th_vec   = dth*nn+shift

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
   D_inc    = f_inc(th_vec)
   n_inc    = nn[np.cos(th_vec)>0]
   coeffs   = solve(M_c2th[n_inc,:],D_inc[n_inc])
   ##############################################
   
   ##############################################
   En_edge  = M_c2ft.dot(coeffs)
   E_edge   = M_c2th.dot(coeffs)
   Sn_edge  = Ds.dot(En_edge)

   out   = {'En_edge': En_edge,'E_edge': E_edge,'Sn_edge': Sn_edge,
            'eig_coeffs':coeffs,'eig_vals':lam,
            'M_c2th':M_c2th,'M_c2ft':M_c2ft,
            'angles':th_vec,'dirspec_inc':D_inc}
   return out
##############################################

alp      = 1.0
N        = 2**8
alp_dis  = 0.0e-5
cg       = 1.0
out      = solve_isotropic_ft(
         alp,N,alp_dis=alp_dis,cg=cg,f_inc=dirspec_inc)

cn       = out['eig_coeffs']
M_c2ft   = out['M_c2ft']
M_c2th   = out['M_c2th']

##############################################
if 0:
   # test edge conditions
   print('Test edge conditions:')
   ang   = out['angles'] 
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

   plt.plot(ang*180/np.pi,out['E_edge'].real)
   plt.plot(ang*180/np.pi,out['dirspec_inc'],'--r')
   plt.show()
##############################################

##############################################
if 1:
   # plot energy vs x:
   npts  = 500
   xx    = linspace(npts,0.0,1.0e3)
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

   fig   = plt.figure()
   ax    = fig.add_subplot(1,1,1)
   ax.plot(xx,E_n.real)
   ax.set_yscale('log')
   plt.show()
##############################################



# return
# 
# #########################################################
# v        = np.arange(0,ny)/float(ny-1) # ny points between 0,1
# Tp0      = 8.0
# Tp1      = 18.0
# Tp_vec   = Tp0+(Tp1-Tp0)*v
# 
# if 0:
#    # do a new calculation:
# 
#    ######################################################
#    # set input ice conc/thickness
#    ICE_MASK = np.zeros((nx,ny))
#    ICE_MASK[np.logical_and(X>0.7*xmin,LANDMASK<1)]  = 1 # i>=24
#    icec  = 0.75*ICE_MASK
#    iceh  = 2.0*ICE_MASK
#    dfloe = 250*ICE_MASK
# 
#    # set input wave fields
#    WAVE_MASK   = np.zeros((nx,ny))
#    WAVE_MASK[X<xmin*0.8]   = 1   # i<=15
#    Hs    = 2.0*WAVE_MASK
#    mwd   = -90.0*WAVE_MASK
# 
#    if 0:
#       Tp = 12.0*WAVE_MASK
#    else:
#       Tp = np.zeros((nx,ny))
#       for j in range(0,ny):
#          Tp[:,j]  = Tp_vec[j]*WAVE_MASK[:,j]
#    ######################################################
# 
#    # set in_fields:
#    in_fields   = {'icec'   : icec,
#                   'iceh'   : iceh,
#                   'dfloe'  : dfloe,
#                   'Hs'     : Hs,
#                   'Tp'     : Tp,
#                   'mwd'    : mwd}
# 
#    # parameters for advection/attenuation
#    SOLVER      = 1
#    ADV_DIM     = 1
#    int_prams   = np.array([SOLVER,ADV_DIM])
# 
#    # do calculation in fortran:
#    out_fields,outdir = Rwim.do_run(RUN_OPT=2,
#                                    in_fields=in_fields,
#                                    int_prams=int_prams)
#    ######################################################
# else:
#    # load results of previous run:
#    out_fields,outdir = Rwim.do_run(RUN_OPT=3)
# #########################################################
# 
# #########################################################
# # calculate maxima of tau_x and get Wmiz:
# taux_max = np.zeros(ny)
# Wmiz     = np.zeros(ny)
# for j in range(0,ny):
#    Dmax  = out_fields['dfloe'][:,j]
#    taux  = out_fields['taux'] [:,j]
#    ##
#    taux_max[j] = taux.max() # max stress in Pa
#    ##
#    miz   = np.zeros(nx)
#    miz[np.logical_and(Dmax>0.0,Dmax<250.0)]  = 1.0
#    Wmiz[j]  = sum(miz*dx/1.0e3)  # MIZ width in km
# #########################################################
# 
# #########################################################
# def plot_diagnostics(Tp_vec,Wmiz,taux_max):
#    # plot MIZ width and taux max
# 
#    figdir   = 'out_io/figs_diag'
# 
#    # MIZ width:
#    fig   = figdir+'/Wmiz.png'
#    f     = Fplt.plot_1d(Tp_vec,Wmiz,['$T_p$, s','$W_{MIZ}$, km'])
#    plt.savefig(fig,bbox_inches='tight',pad_inches=0.05)
#    plt.close()
#    f.clf()
# 
#    # tau_x
#    fig   = figdir+'/taux.png'
#    f     = Fplt.plot_1d(Tp_vec,taux_max,['$T_p$, s','stress ($x$ dir), Pa'])
#    plt.savefig(fig,bbox_inches='tight',pad_inches=0.05)
#    plt.close()
#    f.clf()
# #########################################################
# 
# if 1:
#    print("Tp (s) :")
#    print(Tp_vec)
# 
#    print(" ")
#    print("MIZ midth (km) :")
#    print(Wmiz)
# 
#    print(" ")
#    print("max tau_x (Pa) :")
#    print(taux_max)
# 
#    print(" ")
#    print("Plotting diagnostics...")
#    print(" ")
#    plot_diagnostics(Tp_vec,Wmiz,taux_max)
# #########################################################
# 
# if 1:
#    # plot results from binaries:
#    print(" ")
#    print("Getting inputs/outputs from binaries...")
#    print(" ")
# 
#    ## look at initial fields:
#    print("Plotting initial conditions...")
#    grid_prams              = Fdat.fn_check_grid(outdir) # load grid from binaries
#    ice_fields,wave_fields  = Fdat.fn_check_init(outdir) # load initial conditions from binaries
#    ##
#    figdir   = outdir+'/figs/'
#    Fplt.fn_plot_init(grid_prams,ice_fields,wave_fields,figdir) # plot initial conditions
#    print("Plots in "+figdir+"/init")
#    print(" ")
# 
#    ## look at results:
#    print("Plotting results...")
#    out_fields  = Fdat.fn_check_out_bin(outdir)
#    # Fplt.fn_plot_final(grid_prams,out_fields,figdir)
#    Fplt.fn_plot_final_V1d(grid_prams,Tp_vec,out_fields,figdir)
#    print("Plots in "+figdir+"/final")
# elif 0:
#    # plot results from outputs:
# 
#    ## look at initial fields:
#    print("Plotting initial conditions...")
#    figdir   = outdir+'/figs/'
#    Fplt.fn_plot_init(grid_prams,ice_fields,wave_fields,figdir) # plot initial conditions
#    print("Plots in "+figdir+"/init")
#    print(" ")
# 
#    ## look at results:
#    print("Plotting results...")
#    # Fplt.fn_plot_final(grid_prams,out_fields,figdir)
#    Fplt.fn_plot_final_V1d(grid_prams,Tp_vec,out_fields,figdir)
#    print("Plots in "+figdir+"/final")
