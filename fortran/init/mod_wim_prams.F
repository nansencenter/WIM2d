      module mod_wim_prams

      implicit none

      !! Model Parameters
      real,parameter :: rhowtr   = 1025   ! Ice density      [kg/m^3]
      real,parameter :: rhoice   = 922.5  ! Ice density      [kg/m^3]
      real,parameter :: gravity  = 9.81   ! Gravity          [m/s^2]
      real,parameter :: poisson  = .3     ! Poisson's ratio
      real,parameter :: cohesion = 40.e3  ! Cohesion         [Pa]

      ! Parameters for the floe size distribution
      real,parameter :: Dmin        = 20   ! min floe size [m]
      real,parameter :: xi          = 2    ! [-]
      real,parameter :: fragility   = 0.9  ! [-]

      !! parameters that still need to be set
      !! in set_wim_prams
      real  :: visc_rp  !! Robinson-Palmer viscosity coefficient [Pa/(m/s)] (~13)
      real  :: young    !! young's modulus  [Pa]
      real  :: epsc     !! breaking strain [-]
      real  :: stress_c !! breaking stress [Pa]
      real  :: vbf      !! brine volume fraction [-]
      real  :: sigma_c  !! flexural strength [Pa]

      contains

      subroutine set_wim_prams(BRK_OPT,wim_prams_in)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      implicit none

      integer,intent(in)         :: BRK_OPT
      real,optional,intent(in)   :: wim_prams_in(2)

      integer,parameter    :: YOUNG_OPT   = 1

      real  :: flex_rig_coeff
         !! flex strength
         !! critical stress (for Marcheko criterion)
      real  :: vb,salt,temp
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!Brine vol fraction -> Young's modulus and flexural strength
      if (.true.) then
         !!just set
         vbf   = .1        !! brine volume fraction 
         vb    = vbf*1.0   !! brine volume [ppt]
      else
         !!F&G (1967)
         if (.true.) then !!Dany's estimate
            salt  = 5.0   ! Ice salinity    [psu]
            temp  = -10.0 ! Ice temperature [oC]
         else!!vbf=.1
            temp  = -10.0               ![oC]
            salt  = 18.346940647647006  ![psu]
         end if
         ! Brine volume (Frankenstein and Gardner 1967)
         vb    = salt*(49.185/abs(temp) + 0.532)  ! [ppt]
         vbf   = vb/1.0d3                         ! brine volume fraction 
      end if

      if (present(wim_prams_in)) then
         young    = wim_prams_in(1)
         visc_rp  = wim_prams_in(2)
      else
         !!Robinson-Palmer viscosity parameter
         visc_rp  = 13.

         !!Young's modulus
         if (YOUNG_OPT.eq.1) then
            !!just set Young's modulus
            young = 2.0e9     !!lower ~ Marchenko
            !young = 5.49e9    !!higher ~ Vernon
         else
            young = 10.0e9*(1-3.51*vbf)-1.0e9   ! Young's modulus [Pa]
         end if
      end if

      !!Flexural strength (Timco and O'Brien 1994)
      sigma_c  = 1.76d6*exp(-5.88*sqrt(vbf)) ! [Pa]

      !!Breaking criterion option
      if ((BRK_OPT.eq.0).or.(BRK_OPT.eq.1)) then !!beam test
         !! not used if BRK_OPT==0 (no breaking)
         epsc     = sigma_c/young
         stress_c = young*epsc/(1-poisson**2)
      elseif (BRK_OPT.eq.2) then!!Marchenko's stress criterion (plate)
          !! - convert to strain criterion
         stress_c = 2.6*sigma_c!!=E/(1-nu^2)*strain_c = thin plate (plane stress: \sigma_33=0)
         epsc     = (1-poisson**2)*stress_c/young
      elseif (BRK_OPT.eq.3) then!!Mohr-Coulomb criterion (plate)
          !! - convert to strain criterion
         stress_c = cohesion!!=E/(1-nu^2)*strain_c = thin plate (plane stress: \sigma_33=0)
         epsc     = (1-poisson**2)*stress_c/young
      else
         print*,'Unknown value for BRK_OPT: ',BRK_OPT
         stop
      end if

      !! flex rigidity = flex_rig_coeff*h^3
      flex_rig_coeff = young/12.0/(1-poisson**2)

      print*,'********************************************'
      print*,'in set_wim_prams()...'
      print*,'Brine volume fraction:      ',vbf
      print*,'Youngs modulus (Pa):        ',young
      print*,'BRK_OPT:                    ',BRK_OPT
      print*,'Flexural strength (Pa):     ',sigma_c
      print*,'Breaking stress (Pa):       ',stress_c
      print*,'Breaking strain:            ',epsc
      print*,'Damping (Pa.s/m):           ',visc_rp
      print*,'********************************************'
      print*,' '

      end subroutine set_wim_prams

      end module mod_wim_prams
