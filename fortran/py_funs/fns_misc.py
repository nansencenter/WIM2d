import numpy as np

def SDF_Bretschneider(omega,T_peak,H_sig):

    PI   = np.pi
    T    = 2*PI/omega
    om_m = 2*PI/T_peak
    ##
    f1 = 5.0/16.0*pow(H_sig,2)*pow(om_m,4)
    f2 = 1.0/pow(omega,5)
    f3 = exp(-1.25*pow(T/T_peak,4))
    ##
    S  = f1*f2*f3

    return S

def wt_simpsons(nw,dom=1.):

    wt0 = 4.*np.ones(nw)
    wt0[np.arange(1,nw,2,dtype='int')] = 2.
    wt0[0]  = 1.
    wt0[-1] = 1.

    wt_om = (dom/3.)*wt0
    return wt_om

def theta_dirfrac(th1_,dtheta_,mwd_):
    ## chi=pi/180*(theta-mwd) 
    ##  (convert to radians and centre about mwd, the mean wave direction)
    ## chi1=pi/180*(th1-mwd)
    ## chi2=pi/180*(th2-mwd)
    ## D(chi)=2/pi*cos^2(chi) if chi\in[-pi/2,pi/2] 
    ## D(chi)=0 if chi\in[-pi/2,pi/2] (no backwards waves)
    ## theta_dirfrac  = \int_{chi1}^{chi2}D(chi)dchi
    ##                     = 1 if -chi1=chi2=pi/2
    ##                     = the fraction of the energy going in the directions in [chi1,chi2]

    ##get mwd inside [th1,th1+360)
    mwd = theta_in_range(mwd_,th1_)  #>th1_
    th2 = th1_+dtheta_
    if ((mwd>th2) and (mwd-th2)>abs(mwd-360-th1_)):
        mwd = mwd-360

    th1 = max(mwd-90,th1_)
    th2 = min(mwd+90,th2)
    th2 = max(th1,th2) #make th2>=th1

    chi1 = np.pi*(th1-mwd)/180.
    chi2 = np.pi*(th2-mwd)/180.

    theta_dirfrac = 2*(chi2-chi1)+np.sin(2*chi2)-np.sin(2*chi1)
    theta_dirfrac = theta_dirfrac/2/np.pi

    return theta_dirfrac


def theta_in_range(th_,th1):
    ## th = theta_in_range(th_)
    ## where th1<=th<th1+360.

    th2 = th1+360.0
    if (th_<th1):
        dth   = th1-th_
        njump = np.ceil(dth/360.)
        th    = th_+njump*360.
    elif (th_>th2):
        dth   = th_-th2
        njump = np.ceil(dth/360.)
        th    = th_-njump*360.
    elif (th_==th2):
        th = th1
    else:
        th = th_

    return th


def spectrum_integrals(sdf_dir,freq_vec,wavdir):

    nx,ny,ndir,nw = sdf_dir.shape
    #
    PI   = np.pi
    CDR  = (PI/180.) # degrees to radians
    dirs = -CDR*(wavdir+90.)
    #
    dth  = 2*PI/ndir
    mom0 = np.zeros((nx,ny))
    mom2 = np.zeros((nx,ny))
    Mom0 = np.zeros((nx,ny))

    # calculate the spectral moments:
    if nw==1:
        om = 2*PI*freq_vec[0]
        for wth in range(ndir):
            mom0 = mom0+dth*sdf_dir[:,:,wth,0]
            mom2 = mom2+dth*om*om*sdf_dir[:,:,wth,0]

    else:
        dom   = abs(freq_vec[1]-freq_vec[0])
        wt_om = wt_simpsons(nw,dom)
        for w in range(nw):
            om = 2*PI*freq_vec[w]
            for wth in range(ndir):
                th   = dirs[wth]
                mom0 = mom0+abs( wt_om[w]*dth*sdf_dir[:,:,wth,w]         )
                mom2 = mom2+abs( wt_om[w]*dth*om*om*sdf_dir[:,:,wth,w] )
                Mom0 = Mom0+abs( wt_om[w]*dth*th*sdf_dir[:,:,wth,w]     )
    
    # calculate important quantities from the moments:
    Hs  = 4.*np.sqrt(mom0)
    Tp  = np.zeros((nx,ny))
    mwd = np.zeros((nx,ny))
    for i in range(nx):
        for j in range(ny):
            if mom0[i,j]>0:
                Tp[i,j]  = 2*PI*np.sqrt(mom0[i,j]/mom2[i,j])
                mwd[i,j] = -1./CDR*(Mom0[i,j]/mom0[i,j]) - 90.

    return Hs,Tp,mwd
