      module mod_common_wim
         !!global variables: ice conc/thickness/floe size

      use mod_gridinfo

      implicit none

      integer,parameter :: itest = 24
      integer,parameter :: jtest = 1

      real,dimension(ii,jj)   :: icec,iceh,dfloe,ICE_MASK
      real,dimension(ii,jj)   :: WTR_MASK
      real,dimension(idm,jdm) :: Hs,Tp,mwd,WAVE_MASK  !set in waves_init
      real,dimension(idm,jdm) :: tau_x,tau_y
      !!
      real,parameter          :: unif_c            = .75
      real,parameter          :: unif_h            = 2.0
      real,parameter          :: dfloe_pack_init   = 250.0
      
      !!wave directions and frequencies
      integer,parameter             :: n_wavdir    = 16
      integer,parameter             :: n_wave_freq = 25!should be odd (Simpson's rule)
      real,dimension(n_wave_freq)   :: freq_vec !set in waves_init
      real,dimension(n_wavdir)      :: wavdir   !set in waves_init

      !!directional spectra (incident,spectrum that is propagated)
      real,dimension(ii,jj,n_wavdir,n_wave_freq)   :: sdf_inc  !set in waves_init
      real,dimension(ii,jj,n_wavdir,n_wave_freq)   :: sdf_dir  !propagated in p_WIM2d

      end module mod_common_wim
