      module mod_common_wim

      implicit none
    
      integer           :: itest,jtest
      character(len=80) :: outdir
      integer,parameter :: nmesh_vars  = 6 !no of fields taken in with info about mesh

      integer  :: AUTO_GRID = 0 ! generate grid arrays automatically (0/1)
      integer  :: SCATMOD       ! scattering model to use
      integer  :: ADV_DIM       ! dimension of advection scheme (1/2)
      integer  :: ADV_OPT       ! boundary conditions
      integer  :: STEADY        ! top up waves each time step
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
      
      real           :: CFL
      real,parameter :: dfloe_pack_init   = 300.0

      type wim_params
         ! Model parameters
         real :: young    = 5.49e9 ! Young's modulus [Pa]
         real :: drag_rp  = 13.    ! Robinson-Palmer drag coefficient [Pa.s/m]
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
         real  :: epsc     = 4.9935e-05 !! breaking strain [-]
         real  :: stress_c = 3.0126e+05 !! breaking stress [Pa]
         real  :: vbf      = 1.         !! brine volume fraction [-]
         real  :: sigma_c  = 2.7414e+05 !! flexural strength [Pa]
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

      type atten_outputs
         real   :: atten_nond   ! non-dimensional attenuation (per floe)
         real   :: damping      ! damping (can combine Robinson & Palmer damping or Wang & Shen viscoelasticity)
         real   :: kice         ! wave number in ice
         real   :: kwtr         ! wave number in water
         real   :: int_adm      ! intrinsic admission
         real   :: modT         ! |T|
         real   :: argR         ! arg(R)
         real   :: argT         ! arg(T)
      end type atten_outputs

      type(wim_params) :: wim_info

      end module mod_common_wim
