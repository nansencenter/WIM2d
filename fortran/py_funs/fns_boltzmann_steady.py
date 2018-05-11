import numpy as np
import os
import sys
import struct
import matplotlib.pyplot as plt

# linear algebra
from scipy import linalg as LA
eig    = LA.eig             #
solve = np.linalg.solve # x=solve(A,b) -> solves A*x=b

import WIM2d_f2py     as Mwim
import run_WIM2d      as Rwim
import fns_get_data  as Fdat
import fns_plot_data as Fplt
import fns_inprod     as Finpr


def calc_M_En2E(th_vec):
    # transfer matrix to get between E_n and E(\theta)
    N   = len(th_vec)
    M   = int(round(N/2.0))
    n2m = np.zeros(N)
    for m in range(-M,M):
        if m>=0:
            n = m
        else:
            n = N+m
        n2m[n] = m

    M_ift = 0j*np.zeros((N,N))
    fac   = .5/np.pi

    for n in range(N):
        m          = n2m[n]
        M_ift[:,n] = fac*np.exp(-1j*m*th_vec)

    return M_ift

def calc_M_E2En(th_vec):
    # transfer matrix to get between E_n and E(\theta)
    N      = len(th_vec)
    M      = int(round(N/2.0))
    n2m    = np.zeros(N)
    for m in range(-M,M):
        if m>=0:
            n = m
        else:
            n = N+m
        n2m[n] = m

    M_ft = 0j*np.zeros((N,N))
    fac  = 2.*np.pi/N

    for n in range(N):
        m         = n2m[n]
        M_ft[n,:] = fac*np.exp(1j*m*th_vec)

    return M_ft

def dirspec_inc_spreading(th_vec,inputs=None):

    if inputs is None:
        Hs = 1.
    else:
        Hs  = inputs['Hs']
        mwd = inputs['mwd']
        th0 = -np.pi/180.*(90.+mwd) # from -90 deg -> to 0 rad
                                    # from   0 deg -> to -pi/2 rad

    # incident spectrum:
    # = 2/pi*cos^2(th)
    cc            = np.cos(th_vec-th0)
    D_inc         = 2.0/np.pi*cc**2
    D_inc[cc<0.0] = 0.0
        # integral=1, so this corresponds to Hs=4*sqrt(1)=4

    D_inc = D_inc*pow(Hs/4.,2)
    # print(th_vec)
    # print(D_inc)
    # sys.exit('dirspec_inc_spreading')

    return D_inc

def dirspec_inc_plane(th_vec,inputs):

    if inputs is None:
        Hs  = 1.
        dth = th_vec[1]-th_vec[0]
    else:
        Hs  = inputs['Hs']
        dth = inputs['dth']

    # incident spectrum:
    # = A^2/2\delta(\theta) : approximate with constant in neigbouring cells to \theta=0
    A  = Hs/2./np.sqrt(2.)
    C  = A*A/2./(2.*dth)
    #
    D_inc                          = 0.*th_vec
    D_inc[abs(th_vec)<dth]         = C
    D_inc[abs(2*np.pi-th_vec)<dth] = C
    if D_inc.sum()>C:
        # zero only corresponds to one interval
        # => must double C
        D_inc = 2*D_inc

    return D_inc

# def matrix_isotropic(alp,N):
#
#     Mb = alp*np.ones((N,N))/float(N)
#     for n in range(N):
#         Mb[n,n]  = Mb[n,n]-alp
#
#     # angles:
#     th_vec    = np.linspace(0.,2*np.pi,N+1)
#     th_vec    = th_vec[:-1]  # last element is 2\pi \equiv 0
#     C          = np.diag(np.cos(th_vec))
#     # print(C)
#
#     return Mb,th_vec,C

def solve_sys(Rmat,width=None,Lmat=None):
    # find and sort eigenvalues
    # evr,Ur,evl,Ul,evz,Uz = sort_evals(evals,U)

    # get eigenvalues & eigenvectors
    # - solve Rmat*U=\lambda*Lmat*U
    evals,U = LA.eig(Rmat,b=Lmat)
    evals   = evals.real          # inv(Lmat)*Ds = real, symmetric => eig's should be real

    Ztol      = 1.e-14
    sort_list = np.array(sorted([(np.abs(ev),i) for i,ev in enumerate(evals)]))
    indx      = np.array(sort_list[:,1],dtype='int')
    evals     = evals[indx]
    U         = U[:,indx]
    out       = {'Rmat':Rmat,'Lmat':Lmat}

    # sort evals:
    # - negative eigenvals
    indx = np.logical_and(evals<0,abs(evals)>=Ztol)
    evl  = evals[indx]
    Ul   = U     [:,indx]

    # - 'zero' eigenvals
    if width is not None:
        out.update({'evals_neg':evl,'evecs_neg':Ul})

        # - positive eigenvals
        indx = np.logical_and(evals>0,abs(evals)>=Ztol)
        evr  = evals[indx]
        Ur   = U     [:,indx]
        out.update({'evals_pos':evr,'evecs_pos':Ur})

        # want both 'zero' evals
        evz  = evals[abs(evals)<Ztol]
        Uz   = U[:,abs(evals)<Ztol]
        # evz[:] = 0.

        if 0:
            print(evals)
            print(evr)
            print(evl)
            print(evz)
            sys.exit()
    else:
        # semi-infinite ice sheet
        # - just take the negative one & combine with evl
        indx  = np.logical_and(evals<=0.,abs(evals)<Ztol)
        if np.any(indx):
            evz = np.array([evals[indx[0]]])
            Uz  = np.array([U    [:,indx[0]]]).transpose()
            if 0:
                print(evals)
                print(evl)
                print(evz)
                print(Uz.shape)
                print(Ul.shape)
                sys.exit()
            # evl    = np.concatenate([0.,evl])
            evl = np.concatenate([evz,evl])
            Ul  = np.concatenate([Uz,Ul],1)
            if 0:
                print(evals)
                print(evl)
                print(evz)
                sys.exit()
        elif 0:
            print(evals)
            print(evl)
            sys.exit()

        out.update({'evals_neg':evl,'evecs_neg':Ul})

    if width is None:
        out.update({'M_c2th_0':Ul})
    else:
        N    = len(evals)
        M    = int(N/2.)
        expL = np.exp(evl*width)# = exp(-evr*width)

        # make matrices to convert coefficents to values at edges
        M_c2th_0 = np.zeros((N,N))*0j
        M_c2th_L = np.zeros((N,N))*0j
        if len(evz)>0:
            """
            zero is a repeated eigenvalue,
            so need to find coefficient vector of 'x':
            - look for y=x*U0+V0
              - U0 is one of the original eigenvectors: Rmat*U0=\lambda*Lmat*U0
              - V0 satisfies (Rmat-\lambda*Lmat)*V0=Lmat*U0,
                 so since \lambda=0,
                 Rmat*V0=Lmat*U0
            """
            M  = M-1
            U0 = Uz[:,0]
            V0 = np.linalg.lstsq(Rmat,Lmat.dot(U0))[0]
            Uz[:,1]  = V0
            out.update({'evals_zero':evz,'evecs_zero':Uz})

            """
            General soln for E_n is:
            y  = sum_n a_n*Ul_n*exp(evl_n*x)
                    + sum_n b_n*Ur_n*exp(evr_n*(x-L))
                    + c0*U0+c1*(x*U0+V0)
             => 2+2*(N/2-1)=N unknowns: coeffs=[c0,c1,[a_n],[b_n]]
            """
            # E_n, x=0:
            M_c2th_0[:,0] = U0
            M_c2th_0[:,1] = V0

            # E_n, x=L:
            M_c2th_L[:,0] = U0
            M_c2th_L[:,1] = V0+U0*width

        """
        If no zero eigenvalues,
        then general soln for E_n is:
        y  = sum_n a_n*Ul_n*exp(evl_n*x)
                + sum_n b_n*Ur_n*exp(evr_n*(x-L))
         => 2*(N/2)=N unknowns: coeffs=[[a_n],[b_n]]
        NB if there are some zero eigenvalues,
        rest of solution proceeds in the same way
        as for this case
        """

        # E_n, x=0:
        n0                  = N-2*M
        M_c2th_0[:,n0:n0+M] = Ul
        M_c2th_0[:,n0+M:]   = Ur.dot(np.diag(expL))
        out.update({'M_c2th_0':M_c2th_0})

        # E_n, x=L:
        M_c2th_L[:,n0:n0+M] = Ul.dot(np.diag(expL))
        M_c2th_L[:,n0+M:]   = Ur
        out.update({'M_c2th_L':M_c2th_L})

    return out

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
    K_fou    = np.zeros(Nout)
    K_fou[0] = .5*Kcos[0]
    for n in range(1,min(int(Nout/2.),No2)):
        K_fou[n]       = .5*Kcos[n]
        K_fou[Nout-n]  = .5*Kcos[n]

    # print(Kcos)
    # print(K_fou)
    return K_fou

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

        K_fou[1]   = a/2.0*K_fou[0]
        K_fou[N-1] = a/2.0*K_fou[0]

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

def solve_boltzmann(width=None,alp=1.0,N=8,alp_dis=0.0,cg=1.0,Hs = 1.,f_inc=None):
    # steady state solution - finite width
    # \pa_t.E + c_g.\pa_x.E=S -> c_g.\pa_x.E=S

    # option for kernel:
    OPT = 0 # isotropic scattering
    # OPT    = 1 # add some directionality
    # OPT    = 2 # real coefficients from file

    # print('alp = '+str(alp))
    # print('width = '+str(width))
    # print(np.exp(-alp*width))
    # sys.exit('solve_boltzmann_ft')

    # angles:
    nn     = np.arange(N)
    M      = int(round(N/2.0))
    dth    = 2*np.pi/N
    shift  = dth/2.0
    th_vec = dth*nn+shift

    if 1:
        # LH matrix: c_g*cos(theta)*d/dx
        Lmat  = np.diag(cg*np.cos(th_vec))

        # RH matrix: K
        qtot  = cg*(alp+alp_dis)

        print('Inputs:')
        print(alp_dis,qtot,alp*cg)
        print('\n')
        if OPT==0:
            Rmat = -qtot*np.eye(N)+(alp*cg)/float(N)*np.ones((N,N))
    else:
        # inv of LH matrix: c_g*cos(theta)*d/dx
        Lmi  = np.diag((1./cg)/np.cos(th_vec))
        Lmat = np.eye(N)

        # RH matrix: K
        qtot = cg*(alp+alp_dis)

        print('Inputs:')
        print(alp_dis,qtot,alp*cg)
        print('\n')
        if OPT==0:
            Rmat = Lmi.dot(-qtot*np.eye(N)+(alp*cg)/float(N)*np.ones((N,N)))

    # get eigenvalues & eigenvectors
    # - solve Rmat*U=\lambda*Lmat*U
    soln = solve_sys(Rmat,width=width,Lmat=Lmat)

    if 0:
        # test ode solved correctly
        npts = 500
        xx   = np.linspace(0.0,1.0e3,npts)
        f    = np.zeros((N,npts)) # (fourier modes,position)
        LHS  = np.zeros((N,npts))
        RHS  = np.zeros((N,npts))

        # eigen-pair to test
        for nc in [0]:

            print(90*'#')
            if 0:
                print('testing negative eval')
                print(90*'#')
                e_val = soln['evals_neg'][nc]
                e_vec = soln['evecs_neg'][:,nc]

                for rx in range(npts):
                    x         = xx[rx]
                    ff        = e_vec*np.exp(e_val*x)
                    f[:,rx]   = ff
                    RHS[:,rx] = Rmat.dot(ff)
                    LHS[:,rx] = Lmat.dot(e_val*ff) # Lmat*df
            elif 1 and (width is not None):
                print('testing positive eval')
                print(90*'#')
                e_val = soln['evals_pos'][nc]
                e_vec = soln['evecs_pos'][:,nc]

                for rx in range(npts):
                    x         = xx[rx]
                    ff        = e_vec*np.exp(e_val*(x-width))
                    f[:,rx]   = ff
                    RHS[:,rx] = Rmat.dot(ff)
                    LHS[:,rx] = Lmat.dot(e_val*ff) # Lmat*df

            elif nc==0 and 'evals_zero' in soln.keys():
                print('testing first zero evec')
                print(90*'#')
                e_val = soln['evals_zero'][nc]
                e_vec = soln['evecs_zero'][:,nc]

                for rx in range(npts):
                    x         = xx[rx]
                    ff        = e_vec # constant
                    f[:,rx]   = ff
                    RHS[:,rx] = Rmat.dot(ff)
                    LHS[:,rx] = Lmat.dot(e_val*ff) # Lmat*df

            elif nc==1 and 'evals_zero' in soln.keys():
                print('testing second zero evec')
                print(90*'#')
                e_val = soln['evals_zero'][nc]
                U0    = soln['evecs_zero'][:,0]
                V0    = soln['evecs_zero'][:,1]

                for rx in range(npts):
                    x         = xx[rx]
                    ff        = U0*x+V0
                    f[:,rx]   = ff
                    RHS[:,rx] = Rmat.dot(ff)
                    LHS[:,rx] = Lmat.dot(U0)

            # angle to test
            r_test = 12
            print(nc,e_val)
            print(r_test,f[r_test,[0,-1]])
            print(RHS[r_test,[0,-1]],RHS[r_test,[0,-1]])

            plt.plot(xx,LHS[r_test,:].real)
            plt.plot(xx,RHS[r_test,:].real,'--r')
            plt.show()

        sys.exit()

    # Default incident wave spectrum
    # if f_inc is None:
    #     f_inc     = dirspec_inc_spreading  # spreading around central peak - Hs=1 is default
    if f_inc is None:
        f_inc = dirspec_inc_plane # delta function - numerical approx

    if f_inc is dirspec_inc_plane:
        finc_in = {'Hs':Hs,'dth':dth}

    elif f_inc is dirspec_inc_spreading:
        finc_in = {'Hs':Hs,'mwd':-90}

    # find coeffs of eigenfunction expansion
    # from incident spectrum:

    print('solve_boltzmann: finc, input')
    # print(f_inc)
    # print(finc_in)
    # sys.exit()
    #
    if width is None:
        M_c2th_0 = soln['M_c2th_0']
        soln.pop('M_c2th_0')
        #
        n_inc0 = nn[np.cos(th_vec)>0] # N/2 unknowns
        Mlhs   = M_c2th_0[n_inc0,:] # waves to right only
        #
        D_inc  = f_inc(th_vec,finc_in)
        rhs    = D_inc[n_inc0]        # waves to right only
        coeffs = solve(Mlhs,rhs)
    else:
        M_c2th_0 = soln['M_c2th_0']
        M_c2th_L = soln['M_c2th_L']
        soln.pop('M_c2th_0')
        soln.pop('M_c2th_L')
        #
        n_inc0 = nn[np.cos(th_vec)>0] # N/2 unknowns
        n_incL = nn[np.cos(th_vec)<0] # N/2 unknowns
        M      = len(n_inc0)
        #
        Mlhs        = np.zeros((N,N))*0j
        Mlhs[0:M,:] = M_c2th_0[n_inc0,:] # waves to right only
        Mlhs[M:,:]  = M_c2th_L[n_incL,:] # waves to left  only
        #
        rhs      = np.zeros(N)*0j
        D_inc    = f_inc(th_vec,finc_in)
        rhs[0:M] = D_inc[n_inc0]        # waves to right only
        coeffs   = solve(Mlhs,rhs)

    # outputs
    E_edge_0 = M_c2th_0.dot(coeffs)
    M_E2En   = calc_M_E2En(th_vec)

    edge_lhs = {'E_edge'      : E_edge_0,
                'M_c2th'      : M_c2th_0,
                'dirspec_inc' : D_inc,
                'which_edge'  : 'lhs'}

    inputs    = {'alp_scat'     : alp,
                 'Ndir'         : N,
                 'alp_dis'      : alp_dis,
                 'group_vel'    : cg,
                 'Hs_inc'       : Hs,
                 'f_dirspec_inc': f_inc}

    if width is None:
        out = {'edge_lhs'   : edge_lhs,
               'solution'   : soln,
               'eig_coeffs' : coeffs,
               'M_E2En'     : M_E2En,
               'inputs'     : inputs,
               'angles'     : th_vec}
    else:
        E_edge_L  = M_c2th_L.dot(coeffs)
        inputs.update({'width':width})

        edge_rhs = {'E_edge'      : E_edge_L,
                    'M_c2th'      : M_c2th_L,
                    'dirspec_inc' : 0*D_inc,
                    'which_edge'  : 'rhs'}

        out        = {'edge_lhs'     : edge_lhs,
                        'edge_rhs'   : edge_rhs,
                        'solution'   : soln,
                        'eig_coeffs' : coeffs,
                        'M_E2En'     : M_E2En,
                        'inputs'     : inputs,
                        'angles'     : th_vec}

    return out

def test_edge_cons(out,semiinf=True,lhs=True):
    # test LHS edge conditions (semi-infinite)
    print('Test edge conditions:')

    th_vec = out['angles']
    print(len(th_vec))

    if semiinf:
        figname = 'fig_scripts/figs/SSboltzmann-EdgeCons-semiinf.png'
        edge    = out['edge_lhs']
    elif lhs:
        figname = 'fig_scripts/figs/SSboltzmann-EdgeCons-lhs.png'
        edge    = out['edge_lhs']
    else:
        figname = 'fig_scripts/figs/SSboltzmann-EdgeCons-rhs.png'
        edge    = out['edge_rhs']

    if 1:
        print('angles (deg)')
        print(th_vec[:10]*180/np.pi)
        print('cosines')
        print(np.cos(th_vec)[:10])
        print('Im(E_edge),x=0')
        print(edge['E_edge'].imag[:10])
        print('Re(E_edge),x=0')
        print(edge['E_edge'].real[:10])
        print('Incident spectrum')
        print(edge['dirspec_inc'][:10])

    if 0:
        plt.plot(th_vec*180/np.pi,edge['E_edge'].real)
        plt.plot(th_vec*180/np.pi,edge['dirspec_inc'],'--r')
        plt.show()
    elif 1:
        pobj  = Fplt.plot_1d(th_vec*180/np.pi,edge['E_edge'].imag,
                    labs=['Angle,  degrees','$E$'],linestyle='-',color='b')
        Fplt.plot_1d(th_vec*180/np.pi,edge['E_edge'].real,
                    labs=None,pobj=pobj,linestyle='-',color='k')
        Fplt.plot_1d(th_vec*180/np.pi,edge['dirspec_inc'].real,
                    labs=None,pobj=pobj,linestyle='--',color='r')

        fig,ax,line = pobj
        fig.savefig(figname,bbox_inches='tight',pad_inches=0.05)
        ax.cla()
        plt.close(fig)
        print(' ')
        print('Saving to file : '+figname)
        print(' ')
    return

def calc_expansion(out,xx,L=None,n_test=None):
    # *calculates eigenfunction expansion "Y" at x in "xx"
    #  using eigenvalues from  out["eig_vals"]
    #  and    eigenvectors from out["eig_vals"]
    # *if "out" comes from solve_boltzmann_ft,
    #  "Y" gives the Fourier series coefficents for the
    #  directional energy spectrum E(th)
    # *if "out" comes from solve_boltzmann,
    #  "Y" gives the directional energy spectrum E(th) itself

    soln  = out['solution']
    lam   = soln['evals_neg']
    Zeros = ('evals_zero' in soln.keys())
    cn    = out['eig_coeffs']
    if len(soln['evecs_neg'].shape)==1:
        Nrows = 1
    else:
        Nrows,Nl = soln['evecs_neg'].shape

    if L is None:
        # semi-infinite
        semiinf = True
        No2     = len(cn)
        N       = 2*No2
        Ul      = soln['evecs_neg']

    else:
        # finite width
        semiinf = False
        N       = len(cn)
        No2     = int(N/2.)

    npts  = len(xx)
    if n_test is None:
        n_test = range(Nrows)
    E_n = np.zeros((npts,len(n_test)))*0j

    print('Reconstructing energy from eigenvalues in calc_expansion...')

    if not semiinf:
        if not Zeros:
            # only single eigen-values
            print('only simple eigenvalues')
            print('\n')
            #
            M  = No2
            an = cn[0:M]
            bn = cn[M:]
        else:
            # repeated eigen-value
            print('zero is a repeated eigenvalue')
            print('\n')
            #
            M  = No2-1
            an = cn[2:][0:M]
            bn = cn[2:][M:]
    else:
        print('semi-infinite case')
        print('\n')

    for m in n_test:
        # loop over elements in eigenvectors
        # eg this could be th_vec or Fourier index

        if semiinf:
            # no need to distinguish between repeated eigenvalue case
            E_n[:,m] = 0*xx
            for n in range(npts):
                x    = xx[n]
                cn_x = np.exp(lam*x)*cn # element by element multiplication (lam[i] corresponds to cn[i])
                E_n[n,m] = E_n[n,m] + Ul[m,:].dot(cn_x)
        else:
            # finite width
            if not Zeros:
                # only single eigen-values
                Ul = soln['evecs_neg'][m,:]
                Ur = soln['evecs_pos'][m,:]
                #
                E_n[:,m] = 0.*xx
            else:
                # repeated eigen-value
                U0 = soln['evecs_zero'][m,0]
                V0 = soln['evecs_zero'][m,1]
                Ul = soln['evecs_neg' ][m,:]
                Ur = soln['evecs_pos' ][m,:]
                #
                E_n[:,m] = cn[0]*U0+cn[1]*(U0*xx+V0)

            for n in range(npts):
                x    = xx[n]
                an_x = np.exp(lam*x)*an
                bn_x = np.exp(lam*(L-x))*bn
                E_n[n,m] = E_n[n,m] + Ul.dot(an_x) \
                                          + Ur.dot(bn_x)

    # finished loop over m - return expansion
    return E_n # NB this is E if out comes from solve_boltzmann

def calc_expansion_deriv(out,xx,L=None,n_test=None):
    # *calculates deriv of eigenfunction expansion "Y" at x in "xx"
    #  using eigenvalues from  out["eig_vals"]
    #  and    eigenvectors from out["eig_vals"]
    # *if "out" comes from solve_boltzmann_ft,
    #  "Y" gives the Fourier series coefficents for the
    #  directional energy spectrum E(th)
    # *if "out" comes from solve_boltzmann,
    #  "Y" gives the directional energy spectrum E(th) itself

    soln  = out['solution']
    lam   = soln['evals_neg']
    Zeros = ('evals_zero' in soln.keys())
    cn    = out['eig_coeffs']
    if len(soln['evecs_neg'].shape)==1:
        Nrows = 1
    else:
        Nrows,Nl = soln['evecs_neg'].shape

    if L is None:
        # semi-infinite
        semiinf = True
        No2     = len(cn)
        N       = 2*No2
        Ul      = soln['evecs_neg']

    else:
        # finite width
        semiinf = False
        N       = len(cn)
        No2     = int(N/2.)

    npts  = len(xx)
    if n_test is None:
        n_test = range(Nrows)
    E_n = np.zeros((npts,len(n_test)))*0j

    print('Reconstructing energy from eigenvalues in calc_expansion_deriv...')

    if not semiinf:
        if not Zeros:
            # only single eigen-values
            print('only simple eigenvalues')
            print('\n')
            #
            M  = No2
            an = cn[0:M]
            bn = cn[M:]
        else:
            # repeated eigen-value
            print('zero is a repeated eigenvalue')
            print('\n')
            #
            M  = No2-1
            an = cn[2:][0:M]
            bn = cn[2:][M:]
    else:
        print('semi-infinite case')
        print('\n')

    for m in n_test:
        # loop over elements in eigenvectors
        # eg this could be th_vec or Fourier index

        if semiinf:
            # no need to distinguish between repeated eigenvalue case
            E_n[:,m] = 0*xx
            for n in range(npts):
                x    = xx[n]
                cn_x = np.exp(lam*x)*cn # element by element multiplication (lam[i] corresponds to cn[i])
                E_n[n,m] = E_n[n,m] + Ul[m,:].dot(cn_x*lam)
        else:
            # finite width
            if not Zeros:
                # only single eigen-values
                Ul = soln['evecs_neg'][m,:]
                Ur = soln['evecs_pos'][m,:]
                #
                E_n[:,m] = 0.*xx
            else:
                # repeated eigen-value
                U0 = soln['evecs_zero'][m,0]
                V0 = soln['evecs_zero'][m,1]
                Ul = soln['evecs_neg' ][m,:]
                Ur = soln['evecs_pos' ][m,:]
                #
                E_n[:,m] = cn[1]*U0

            for n in range(npts):
                x    = xx[n]
                an_x = np.exp(lam*x)*an
                bn_x = np.exp(lam*(L-x))*bn
                E_n[n,m] = E_n[n,m] + Ul.dot(an_x*lam) \
                                          - Ur.dot(bn_x*lam)

    # finished loop over m - return expansion
    return E_n # NB this is E if out comes from solve_boltzmann

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
        semiinf   = True
        L         = 500.0e3 # plotting limit
        figname   = fdir+'SSboltzmann_E'+str(n_test)+'_profile-semiinf.png'
        figname2  = fdir2+'SSboltzmann_IsoFrac_profile-semiinf.png'
        edge_keys = ['lhs']
    else:
        semiinf   = False
        L         = width
        figname   = fdir+'SSboltzmann_E'+str(n_test)+'_profile-L'+str(width)+'.png'
        figname2  = fdir2+'SSboltzmann_IsoFrac_profile-L'+str(width)+'.png'
        edge_keys = ['lhs','rhs']

    # plot energy vs x:
    npts = 5000
    xx   = np.linspace(0.0,L,npts)
    E_n  = 0.0j*xx
    #
    cn      = out['eig_coeffs']
    M_c2th  = out['edge_lhs']['M_c2th']
    alp     = out['inputs'][0]
    alp_dis = out['inputs'][2]
    th_vec  = out['angles']

    if 'M_E2En' in out.keys():
        # have eigenvectors in position space;
        # transform n=0 into Fourier space
        inp0    = {}
        inp_all = {}
        for key in out.keys():
            if type(out[key])==type([]):
                print(key+' (list)')
                inp0.update({key      : list(out[key])})
                inp_all.update({key  : list(out[key])})
            else:
                print(key+' (np array)')
                inp0.update({key:out[key].copy()})
                inp_all.update({key:out[key].copy()})

        Mt   = out['M_E2En']
        soln = out['solution']
        for key in soln.keys():
            if 'evecs' in key:
                print('Taking FT of: '+key)
                print(Mt.shape,inp_all['solution'][key].shape,inp0['solution'][key].shape)
                EV = soln[key]
                inp0['solution'][key] = np.array([Mt[n_test,:].dot(EV)]) # convert back to rank 2
                print(Mt.shape,inp_all['solution'][key].shape,inp0['solution'][key].shape)
                inp_all['solution'][key] = Mt.dot(EV)
    else:
        inp0    = out.copy()
        inp_all = out
        for key in out['solution'].keys():
            if 'evecs' in key:
                print('Taking '+str(n_test)+'-th row of: '+key)
                inp0['solution'][key] = np.array([inp0['solution'][key][n_test,:]]) # convert back to rank 2

    # calc Fourier coefficients at x:
    if semiinf:
        # semi-infinite
        E_n   = calc_expansion(out,xx,L=None,n_test=[n_test])
        No2   = len(cn) # N/2
        N     = 2*No2
        En_si = E_n

        En_all   = calc_expansion(out,xx,L=None)
        E_coh    = 0*xx # coherent energy: E going to right
        iso_frac = 0*xx # amount of isotropy measure
        dth      = th_vec[1]-th_vec[0]
        for n in range(npts):
            E_th     = out['M_En2E'].dot(En_all[n,:])
            E_f      = E_th[abs(th_vec)<dth].real # two closest intervals around 0
            E_coh[n] = dth*(E_f.sum()) # integrate over these intervals
            #
            e0          = abs(E_n[n])*abs(E_n[n])
            ea          = abs(En_all[n,:])*abs(En_all[n,:])
            iso_frac[n] = e0/np.sum(ea)
    else:
        # finite width
        E_n = calc_expansion(inp0,xx,L=L)
        N   = len(cn)
        No2 = int(N/2.)
        if 0:
            comp_si= 1
            cg     = out['inputs'][3]
            Hs     = out['inputs'][4]
            f_inc  = out['inputs'][5]
            out_si = solve_boltzmann_ft_semiinf(alp=alp,N=N,alp_dis=alp_dis,cg=cg,Hs=Hs,f_inc=f_inc)
            En_si  = calc_expansion(out_si,xx,L=None,n_test=[n_test]) # semi-infinite profile
        #
        En_all   = calc_expansion(inp_all,xx,L)
        E_coh    = 0*xx # coherent energy: E going to right
        iso_frac = 0*xx # amount of isotropy measure
        dth      = th_vec[1]-th_vec[0]
        M_ift    = calc_M_En2E(th_vec)

        for n in range(npts):
            E_th     = M_ift.dot(En_all[n,:])
            E_f      = E_th[abs(th_vec)<dth].real # two closest intervals around 0
            E_coh[n] = dth*(E_f.sum()) # integrate over these intervals
            #
            e0 = abs(E_n[n])*abs(E_n[n])
            ea = abs(En_all[n,:])*abs(En_all[n,:])
            iso_frac[n] = e0/np.sum(ea)

    th_vec = out['angles']
    dth    = 2*np.pi/float(N)
    fwd_mask = np.zeros(N)
    fwd_mask[np.cos(th_vec)>0] = 1.0
    if 1:
        if f_inc is None:
            #dirspec_plane or analytical delta function
            A      = Hs/2./np.sqrt(2.)
            Ec_tst = A*A/2.
        else:
            #dirspec_spreading
            Ec_tst = 1./np.pi*(dth+.5*np.sin(2*dth))*pow(Hs/4.,2)

        print('\n')
        print('Check coherent energy (converted to Hs)')
        print('(analytic integration at x=0)')
        print(4*np.sqrt(E_coh[0]),4*np.sqrt(Ec_tst))
        print('\n')
    #
    if 0:
        for key in edge_keys:
            print('************************************************')
            print('check edge E')
            edge   = out['edge_'+key]
            E_edge = edge['E_edge']
            E0_p   = dth*(fwd_mask*E_edge).sum()
            E0_m   = dth*((1-fwd_mask)*E_edge).sum()

            if n_test==0:
                print('E0 at edge: '+edge['which_edge'])
                print('E0 fwd    = '+str(E0_p))
                print('E0 back   = '+str(E0_m))
                # print('E0 total    = '+str(edge['En_edge'][0]))
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
        fig = plt.figure()
        ax  = fig.add_subplot(1,1,1)
        ax.plot(xx,E_n.real)
        ax.set_yscale('log')
        plt.show()
    elif 0:#n_test==0:
        fig  = plt.figure()
        ax   = fig.add_subplot(1,1,1)
        pobj = [fig,ax]
        if (not semiinf) and (comp_si==1):
            Fplt.plot_1d(xx/1.0e3,4*np.sqrt(En_si.real),pobj=pobj,
                    plot_steps=False,linestyle='-',color='g')
        Fplt.plot_1d(xx/1.0e3,4*np.sqrt(E_n.real),pobj=pobj,
                    plot_steps=False,labs=['$x$, km','$H_s$, m'],linestyle='-',color='k')
        Fplt.plot_1d(xx/1.0e3,4*np.sqrt(E_coh),pobj=pobj,
                    plot_steps=False,linestyle='--',color='r')

        # check vs exponential decay
        Fplt.plot_1d(xx/1.0e3,Hs*np.exp(-alp/2.*xx),pobj=pobj,
                    plot_steps=False,linestyle='--',color='b')

        # check vs equi-partition of energy
        E_eq  = En_si[-1].real/float(N)
        Fplt.plot_1d(xx/1.0e3,4*np.sqrt(E_eq)+0*xx,pobj=pobj,
                    plot_steps=False,linestyle='--',color='g')

        if lam[0]>0.0:
            # dissipation so plot y with log scale (exponential decay)
            ax.set_yscale('log')

        plt.savefig(figname,bbox_inches='tight',pad_inches=0.05)
        ax.cla()
        plt.close(fig)
        print(' ')
        print('Saving to file : '+figname)
        print(' ')

        # plot amount of isotropy:
        fig = plt.figure()
        Fplt.plot_1d(xx/1.0e3,iso_frac,pobj=pobj,
                    plot_steps=False,labs=['$x$, km','Isotropic fraction'],linestyle='-',color='k')
        Fplt.plot_1d(xx/1.0e3,1-np.exp(-alp*xx),pobj=pobj,
                    plot_steps=False,linestyle='--',color='r')
        plt.savefig(figname2,bbox_inches='tight',pad_inches=0.05)
        plt.close()
        fig.clf()
        print(' ')
        print('Saving to file : '+figname2)
        print(' ')
    elif 1:
        # simple plot
        fig = plt.figure()
        ax  = fig.add_subplot(1,1,1)
        if 0:
            # of E
            Fplt.plot_1d(xx/1.0e3,E_n.real,pobj=[fig,ax],
                    plot_steps=False,labs=['$x$, km','$E$, m$^2$'],linestyle='-',color='k')
        else:
            # of Hs
            Fplt.plot_1d(xx/1.0e3,4*np.sqrt(E_n.real),pobj=[fig,ax],
                    plot_steps=False,labs=['$x$, km','$H_s$, m'],linestyle='-',color='k')

        if alp_dis>0.0:
            ax.set_yscale('log')

        fig.savefig(figname,bbox_inches='tight',pad_inches=0.05)
        ax.cla()
        plt.close(fig)
        print(' ')
        print('Saving to file : '+figname)
        print(' ')

    out = {'E_n':E_n,'x':xx,'n':n_test,'E_coh':E_coh}

    # save data to file
    if n_test==0:
        Hs   = 4*np.sqrt(E_n.real)
        Hs_f = 4*np.sqrt(E_coh)

        ddir = 'fig_scripts/figs/profiles/datfiles'
        if not os.path.exists(ddir):
            os.mkdir(ddir)

        if width is None:
            out_file = ddir+'/steady_semiinf.dat'
        else:
            out_file = ddir+'/steady_L'+str(width)+'.dat'

        print('\n')
        print('Saving data points to '+out_file)
        of1    = open(out_file,'w')
        of1.write('x, m         Hs, m         Hs (coherent), m\n')
        for n in range(len(xx)):
            of1.write('%f    %f     %f\n'%(xx[n],Hs[n],Hs_f[n]))

        of1.close()

    return out

def boltz_main_outputs(N=2**5,h=2.,period=12.,youngs=5.e9,\
          visc_rp=0,width=100.,Hs_inc=3.,conc=.7,dmean=100.,ax=None):

    # N            : no of directions
    # h            : # thickness [m]
    # period     : # period [s]
    # youngs     : # period [s]
    # visc_rp    : # R-P damping parameter [m^s/s]
    # Hs          : # sig wave height [m]
    # conc        : # concentration [0-1]
    # dmean      : # mean floe size [m]

    om        = 2*np.pi/period
    gravity   = 9.81
    atten_in  = np.array([h,om,youngs,visc_rp])
    atten_out = Mwim.atten_youngs(atten_in)
    # print(atten_in)
    # print(atten_out)

    alp     = conc/dmean*atten_out[4]  # scattering "attenuation" [m^{-1}]
    alp_dis = 2*conc*atten_out[0]        # damping [m^{-1}]
    kwtr    = atten_out[2]
    cp      = om/kwtr # phase vel (open water) [m/s]
    cg      = cp/2.    # group vel (open water, inf depth relation) [m/s]

    print('\n')
    print('width = '+str(width)+'m')
    print('\n')

    out    = Fbs.solve_boltzmann(width=width,
                                          alp=alp,N=N,alp_dis=alp_dis,cg=cg,Hs=Hs_inc,
                                          f_inc=dirspec_inc_spreading)

    th_vec  = out['angles']
    dth      = 2*np.pi/N

    xx     = np.linspace(0,width,num=800)
    E_th  = Fbs.calc_expansion(out,xx,L=width)
    #
    E0        = dth*E_th.sum(1) # integral over directions (freq spec)
    Hs_steady = 4*np.sqrt(E0.real)
    #
    S_th        = E_th.dot(out['solution']['Rmat'].transpose()) # source term
    S_cos       = S_th.dot(dtheta*np.cos(out['angles'])) # integral with cos over directions
    rhow        = 1025 # kg/m^3
    gravity     = 9.81 # m/s^2
    taux_steady = -(rhow*gravity/cp)*S_cos.real

    #test plot if desired
    if ax is not None:
        ax.plot(xx/1.e3,Hs_out)

    return xx,Hs_out
