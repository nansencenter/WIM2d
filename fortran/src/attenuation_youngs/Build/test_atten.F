      program test_atten

      use mod_atten_young, only: RTparam_outer
      use mod_common_wim

      implicit none

      integer, parameter   :: nw = 25
      real, parameter      :: iceh  = 2.0
      real                 :: freq,inputs(2)
      real,parameter       :: Tmin  = 2.5
      real,parameter       :: Tmax  = 25.0
      real                 :: fmin,fmax,df,T
      real                 :: om,PI
      integer              :: w
      !!
      !real  :: kice,kwtr,int_adm,modT,argR,argT
      type(atten_outputs)  :: outputs
      real  :: kw_inf,wlng,guess
     
      PI    = 4.0*atan(1.0)
      !!
      fmin  = 1./Tmax
      fmax  = 1./Tmin
      df    = (fmax-fmin)/(nw-1.0)
      do w=1,nw
         freq  = fmin+(w-1)*df
         om    = 2*PI*freq
         T     = 1.0/freq
         !!
         inputs(1)   = iceh
         inputs(2)   = om
         !inputs(3)   = wim_info%young
         !inputs(4)   = wim_info%drag_rp

         !! initial guess for ice wavenumber
         kw_inf   = om**2/wim_info%gravity
         wlng     = 2*PI/kw_inf
         if (w.eq.1) guess = kw_inf

         !! get atten coeff
         call RTparam_outer(outputs,inputs,guess)
!        call RTparam_outer(damping,kice,kwtr,int_adm,                  &
!    &                      alp_nd,modT,argR,argT,                      &
!    &                      inputs,guess)

         !! update guess for ice wavenumber
         guess = outputs%kice

         !!check attenuation coeff's, wavenumbers
         if (w.eq.1) then
            print*,'***********************************************'
            print*,'check outputs from RTparam_outer:'
            print*,'T,h',T,iceh
            print*,'atten',outputs%atten_nond,outputs%damping
            print*,'ki,kw,2pi/wlng_wtr/guess'
     &         ,outputs%kice,outputs%kwtr,kw_inf,guess
            print*,'|T|,argRT_s',outputs%modT
     &         ,outputs%argR,outputs%argT,outputs%int_adm
            print*,'***********************************************'
            print*,' '
         end if

         print*,'T (s), atten (per floe), damping (/m):'                &
     &         ,T,outputs%atten_nond,outputs%damping
      end do

      end program test_atten
