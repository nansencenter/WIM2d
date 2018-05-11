import numpy  as np
from matplotlib import pyplot as plt

def ipmat_cos(N,theta_in=None):
    # for \theta\in[0,\pi],
    # expansion of f(\theta)=f_0/2+\sum_{n=1}^{N-1}f_n*cos(n\theta)
    # f_m = 2/\pi\int_0^\pi f(\theta)*cos(m\theta)d\theta
    #      ~ 2/\pi*dth*\sum_{n=0}^{N-1}f(\theta_n)*cos(m\theta_n);
    # in matrix form:
    # fvec    = M_ipcos*f(theta_vec)
    #
    # if theta_in is not None:
    # - do interpolation at values of theta_in

    if theta_in is None:
        M_ipcos  = np.zeros((N,N))
        dth        = np.pi/N
        th         = dth*np.arange(N)+.5*dth # need to be in the centre of the interval
        for n in range(N):
            fac                = (2./np.pi)*dth
            M_ipcos[n,:]    = fac*np.cos(n*th)

        # output is inner product matrix
        # and points where function values are needed
        out    = [M_ipcos,th]
    else:
        M      = len(theta_in)
        M_cos = np.zeros((M,N))
        dth    = np.pi/N
        th     = dth*np.arange(N)
        for n in range(N):
            if n==0:
                fac    = 0.5
            else:
                fac    = 1.
            M_cos[:,n]  = fac*np.cos(n*theta_in)

        # output is matrix to do interpolation
        out    = [M_cos]

    return out

def ipmat_exp(N,theta_in=None):
    # for \theta\in[0,\pi],
    # expansion of f(\theta)=1/(2\pi)\sum_{n=0}^{N-1}f_n*exp(-1j*n\theta)
    # f_m = \int_0^{2\pi} f(\theta)*exp(1j*m\theta)d\theta
    #      ~ dth*\sum_{n=0}^{N-1}f(\theta_n)*exp(1j*m\theta_n);
    # in matrix form:
    # fvec    = M_ip*f(theta_vec)
    #
    # if theta_in is not None:
    # - do interpolation at values of theta_in

    if theta_in is None:
        M_ip  = np.zeros((N,N))*0j
        dth    = 2*np.pi/N
        th     = dth*np.arange(N)
        for n in range(N):
            M_ip[n,:]    = dth*np.exp(1j*n*th)

        # output is inner product matrix
        # and points where function values are needed
        out    = [M_ip,th]
    else:
        M      = len(theta_in)
        M_exp = np.zeros((M,N))*0j
        fac    = .5/np.pi
        for n in range(N):
            M_exp[:,n]  = fac*np.exp(-1j*n*theta_in)

        # output is matrix to do interpolation
        out    = [M_exp]

    return out

def test_ipmat_cos(nc=4):

    N  = np.max((10,4*nc))
    # test cosine expansion of incident wave spectrum:
    out        = ipmat_cos(N)
    M_ipcos  = out[0]
    th_in    = out[1]
    if nc==0:
        f  = .5+0*th_in 
    else:
        f  = np.cos(nc*th_in)

    fc = M_ipcos.dot(f)
    n0 = np.max((0,nc-4))
    n1 = np.min((N,nc+4))
    print(fc[n0:n1])
    #
    out2  = ipmat_cos(N,theta_in=th_in)
    g      = out2[0].dot(fc)
    #
    fig        = plt.figure()
    ax         = fig.add_subplot(111)
    ax.plot(th_in,f,'-b')
    ax.plot(th_in,g,'--r')
    plt.show()
    fig.clf()

def test_ipmat_exp(nc=4):

    N  = np.max((10,4*abs(nc)))
    if nc<0:
        nc = N+nc

    # test cosine expansion of incident wave spectrum:
    out    = ipmat_exp(N)
    M_ip  = out[0]
    th_in = out[1]
    #
    f      = .5/np.pi*np.exp(-1j*nc*th_in)
    fc = M_ip.dot(f)
    #
    n0 = np.max((0,nc-4))
    n1 = np.min((N,nc+4))
    print('all coeffs')
    print(fc)
    print('this should be 1')
    print(fc[nc])
    #
    out2  = ipmat_exp(N,theta_in=th_in)
    g      = out2[0].dot(fc)
    #
    fig        = plt.figure()
    ax         = fig.add_subplot(1,2,1)
    ax.plot(th_in,f.real,'-b')
    ax.plot(th_in,g.real,'--r')
    ax2    = fig.add_subplot(1,2,2)
    ax2.plot(th_in,f.imag,'-b')
    ax2.plot(th_in,g.imag,'--r')
    plt.show()
    fig.clf()
