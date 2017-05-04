      module mod_common_wim
         !!global variables: ice conc/thickness/floe size

      use mod_gridinfo

      implicit none

      include "common_wim.h"

      integer           :: itest,jtest
      integer,parameter :: nmesh_vars  = 6 !no of fields taken in with info about mesh

      real,dimension(ii,jj)   :: icec,iceh,dfloe,ICE_MASK
      real,dimension(ii,jj)   :: WTR_MASK
      real,dimension(idm,jdm) :: Hs,Tp,mwd,WAVE_MASK  !set in waves_init
      real,dimension(idm,jdm) :: tau_x,tau_y
      
      !! set number of wave directions and frequencies in wave_info.h
      !! - n_wavdir,n_wave_freq (should be odd - Simpson's rule)
      include "wave_info.h"

      real,dimension(n_wave_freq)   :: freq_vec !set in waves_init
      real,dimension(n_wavdir)      :: wavdir   !set in waves_init

      !!directional spectra (incident,spectrum that is propagated)
      real,dimension(ii,jj,n_wavdir,n_wave_freq)   :: sdf_inc  !set in waves_init
      real,dimension(ii,jj,n_wavdir,n_wave_freq)   :: sdf_dir  !propagated in p_WIM2d

      character(len=80) :: outdir

      integer  :: SCATMOD,ADV_DIM,ADV_OPT,STEADY
      integer  :: DO_CHECK_PROG,DO_CHECK_INIT,DO_CHECK_FINAL
      integer  :: BRK_OPT  !breaking criterion to use
                           ! 0: no breaking; 1: Vernon; 2: Marchenko; 3: Mohr-Coulomb
      integer  :: DO_ATTEN !can turn off attenuation to test advection
      integer  :: FSD_OPT  ! which FSD parameterisation to use (0,1)
                           ! 0: RG (Dumont et al, 2011);1: smooth FSD

      logical  :: INTERP_MESH          = .false.
      logical  :: TAKE_MAX_WAVES       = .false. ! this is independent of inputs on mesh
      logical  :: NEED_TAKE_MAX_WAVES  = .false. ! even if TAKE_MAX_WAVES is true,
                                                 !  this option will not be needed if
                                                 !  interp_mesh is false (eg if no ice present)
#if 0
      real, dimension(:,:), allocatable ::
     &   mom0_max,mom2_max,var_strain_max
#else
      real, dimension(ii,jj) ::
     &   mom0_max,mom2_max,var_strain_max
#endif

      real     :: CFL

      end module mod_common_wim
