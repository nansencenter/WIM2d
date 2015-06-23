###############################################################
def SDF_Bretschneider(omega,T_peak,H_sig):
   import numpy as np

   PI    = np.pi
   T     = 2*PI/omega
   om_m  = 2*PI/T_peak
   ##
   f1 = 5.0/16.0*pow(H_sig,2)*pow(om_m,4)
   f2 = 1.0/pow(omega,5)
   f3 = exp(-1.25*pow(T/T_peak,4))
   ##
   S  = f1*f2*f3

   return S
###############################################################

###############################################################
def wt_simpsons(nw,dom=1.):

   wt0                                 = 4.*np.ones(nw)
   wt0[np.arange(1,nw,2,dtype='int')]  = 2.
   wt0[0]                              = 1.
   wt0[-1]                             = 1.

   wt_om = (dom/3.)*wt0
   return wt_om
###############################################################

###############################################################
def spectrum_integrals(sdf_dir,freq_vec,wavdir):
   import numpy as np

   nx,ny,ndir,nw  = sdf_dir.shape
   PI             = np.pi
   #
   CDR            = (PI/180.) # degrees to radians
   dirs           = -CDR*(wavdir+90.)
   #
   dth   = 2*PI/ndir
   mom0  = np.zeros((nx,ny))
   mom2  = np.zeros((nx,ny))
   Mom0  = np.zeros((nx,ny))

   # calculate the spectral moments:
   if nw==1:

      om = 2*PI*freq_vec[0]
      for wth in range(ndir):
         mom0  = mom0+dth*sdf_dir[:,:,wth,0]
         mom2  = mom2+dth*om*om*sdf_dir[:,:,wth,0]

   else:

      dom   = abs(freq_vec[1]-freq_vec[0])
      wt_om = wt_simpsons(nw,dom)
      for w in range(nw):
         om = 2*PI*freq_vec[w]
         for wth in range(ndir):
            th    = dirs[wth]
            mom0  = mom0+abs( wt_om[w]*dth*sdf_dir[:,:,wth,w]       )
            mom2  = mom2+abs( wt_om[w]*dth*om*om*sdf_dir[:,:,wth,w] )
            Mom0  = Mom0+abs( wt_om[w]*dth*th*sdf_dir[:,:,wth,w]    )
   
   # calculate important quantities from the moments:
   Hs    = 4.*np.sqrt(mom0)
   Tp    = np.zeros((nx,ny))
   mwd   = np.zeros((nx,ny))
   for i in range(nx):
      for j in range(ny):
         if mom0[i,j]>0:
            Tp[i,j]  = 2*PI*np.sqrt(mom0[i,j]/mom2[i,j])
            mwd[i,j] = -1./CDR*(Mom0[i,j]/mom0[i,j]) - 90.

   return Hs,Tp,mwd
###############################################################
