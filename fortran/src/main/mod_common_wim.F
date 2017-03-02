      module mod_common_wim
         !!global variables: ice conc/thickness/floe size

      use mod_gridinfo

      implicit none

      integer           :: itest,jtest
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


      type wim_params
         ! Model parameters
         real :: young             ! Young's modulus [Pa]
         real :: drag_rp           ! Robinson-Palmer drag coefficient [Pa.s/m]
         real :: visc_ws  = 0.     ! viscosity for Wang & Shen (2010) viscoelastic rheology
                                   ! (Mosig et al, 2016, approximation) [Pa.s]
         real :: rhowtr   = 1025   ! Ice density      [kg/m^3]
         real :: rhoice   = 922.5  ! Ice density      [kg/m^3]
         real :: gravity  = 9.81   ! Gravity          [m/s^2]
         real :: poisson  = .3     ! Poisson's ratio

         ! Mohr-Coulomb failure criterion (BRK_OPT=3)
         real :: cohesion = 1.e5   ! Cohesion         [Pa]
         real :: friction = .7     ! Tangent of friction angle

         !breaking criterion
         real  :: epsc     !! breaking strain [-]
         real  :: stress_c !! breaking stress [Pa]
         real  :: vbf      !! brine volume fraction [-]
         real  :: sigma_c  !! flexural strength [Pa]
         real  :: flex_rig_coeff
            !! flex rigidity = flex_rig_coeff*h^3
            !! critical stress (for Marchenko criterion)

         ! Parameters for the floe size distribution
         ! TODO add to infile_nonIO.txt and param_vec (input from mex
         ! function)
         real  :: Dmin        = 20   ! min floe size [m]
         real  :: xi          = 2    ! [-]
         real  :: fragility   = 0.9  ! [-]
         real  :: Dthresh     = 200. ! [m]
         real  :: cice_min    = 0.05 ! used to set ICE_MASK (=1 if cice>cice_min)
      end type wim_params

      type(wim_params) :: wim_info

      end module mod_common_wim
