      module mod_wim_prams

      implicit none

      !! Model Parameters
      real,parameter :: rhowtr   = 1025   ! Ice density      [kg/m^3]
      real,parameter :: rhoice   = 922.5  ! Ice density      [kg/m^3]
      real,parameter :: gravity  = 9.81   ! Gravity          [m/s^2]
      real,parameter :: poisson  = .3     ! Poisson's ratio
      real,parameter :: visc_rp  = 13     ! Robinson-Palmer viscosity coefficient [Pa/(m/s)]

      ! Parameters for the floe size distribution
      real,parameter :: Dmin        = 20   ! min floe size [m]
      real,parameter :: xi          = 2    ! [-]
      real,parameter :: fragility   = 0.9  ! [-]

      !! parameters that still need to be set
      !! in set_wim_prams
      real  :: young,epsc!!young's modulus,breaking strain

      contains

      subroutine set_wim_prams()

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      implicit none

      integer,parameter :: BRK_OPT     = 1
      integer,parameter :: YOUNG_OPT   = 1

      real  :: sigma_c,stress_c,flex_rig_coeff
         !! flex strength
         !! critical stress (for Marcheko criterion)
      real  :: vbf,vb,salt,temp
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

      !!Young's modulus
      if (YOUNG_OPT.eq.1) then
         !!just set Young's modulus
         young = 2.0d9     !!lower ~ Marchenko
         !young = 5.49d9    !!higher ~ Vernon
      else
         young = 10.0d9*(1-3.51*vbf)-1.0d9   ! Young's modulus [Pa]
      end if

      !!Flexural strength (Timco and O'Brien 1994)
      sigma_c  = 1.76d6*exp(-5.88*sqrt(vbf)) ! [Pa]

      !!Breaking criterion option
      if (BRK_OPT.eq.1) then !!beam test
         epsc  = sigma_c/young
      else!!Marchenko's stress criterion
          !! - convert to strain criterion
         stress_c = 2.6*sigma_c!!=E/(1-nu^2)*strain_c = thin plate (plane stress: \sigma_33=0)
         epsc     = (1-poisson**2)*stress_c/young
      end if

      !! flex rigidity = flex_rig_coeff*h^3
      flex_rig_coeff = young/12.0/(1-poisson**2)

      print*,'********************************************'
      print*,'in set_wim_prams()...'
      print*,'Brine volume fraction:      ',vbf
      print*,'Youngs modulus (Pa):        ',young
      print*,'Flexural strength (Pa):     ',sigma_c
      print*,'Breaking strain:            ',epsc
      print*,'Damping (Pa.s/m):           ',visc_rp
      print*,'********************************************'
      print*,' '

      end subroutine set_wim_prams

      end module mod_wim_prams
