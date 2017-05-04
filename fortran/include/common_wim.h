      real,parameter    :: dfloe_pack_init   = 300.0
      
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

      type(wim_params) :: wim_info
