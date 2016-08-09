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

      type wim_params
         ! Model parameters
         real :: young             ! Young's modulus
         real :: drag_rp           ! Robinson-Palmer drag coefficient
         real :: visc_ws  = 0.     ! viscosity for Mosig et al viscoelastic rheology
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
         real  :: Dmin        = 20   ! min floe size [m]
         real  :: xi          = 2    ! [-]
         real  :: fragility   = 0.9  ! [-]
      end type wim_params

      type(wim_params) :: wim_info

      end module mod_common_wim
