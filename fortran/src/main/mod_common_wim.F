      module mod_common_wim
         !!global variables: ice conc/thickness/floe size

      use mod_gridinfo

      implicit none

      integer,parameter :: itest       = 24
      integer,parameter :: jtest       = 1
      integer,parameter :: nmesh_vars  = 6 !no of fields taken in with info about mesh

      real,dimension(ii,jj)   :: icec,iceh,dfloe,ICE_MASK
      real,dimension(ii,jj)   :: WTR_MASK
      real,dimension(idm,jdm) :: Hs,Tp,mwd,WAVE_MASK  !set in waves_init
      real,dimension(idm,jdm) :: tau_x,tau_y
      !!
      real,parameter          :: dfloe_pack_init   = 300.0
      
      !! set number of wave directions and frequencies in wave_info.h
      !! - n_wavdir,n_wave_freq (should be odd - Simpson's rule)
      include "wave_info.h"

      real,dimension(n_wave_freq)   :: freq_vec !set in waves_init
      real,dimension(n_wavdir)      :: wavdir   !set in waves_init

      !!directional spectra (incident,spectrum that is propagated)
      real,dimension(ii,jj,n_wavdir,n_wave_freq)   :: sdf_inc  !set in waves_init
      real,dimension(ii,jj,n_wavdir,n_wave_freq)   :: sdf_dir  !propagated in p_WIM2d

      end module mod_common_wim