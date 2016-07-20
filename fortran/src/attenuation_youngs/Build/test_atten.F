      program test_atten

      use mod_wim_prams

      implicit none

      integer, parameter   :: nw = 25
      real, parameter      :: iceh  = 2.0
      real                 :: freq,alp_nd,damping,inputs(4)
      real,parameter       :: Tmin  = 2.5
      real,parameter       :: Tmax  = 25.0
      real                 :: fmin,fmax,df,T
      real                 :: om,PI
      integer              :: w
      !!
      real  :: kice,kwtr,int_adm,modT,argR,argT,guess
      real  :: kw_inf,wlng
     
      PI    = 4.0*atan(1.0)
      call set_wim_prams()
      !!
      fmin  = 1/Tmax
      fmax  = 1/Tmin
      df    = (fmax-fmin)/(nw-1.0)
      do w=1,nw
         freq  = fmin+(w-1)*df
         om    = 2*PI*freq
         T     = 1.0/freq
         !!
         inputs(1)   = iceh
         inputs(2)   = om
         inputs(3)   = young
         inputs(4)   = visc_rp

         !! initial guess for ice wavenumber
         kw_inf   = om**2/gravity
         wlng     = 2*PI/kw_inf
         if (w.eq.1) guess = kw_inf

         !! get atten coeff
         call RTparam_outer(damping,kice,kwtr,int_adm,                  &
     &                      alp_nd,modT,argR,argT,                      &
     &                      inputs,guess)

         !! update guess for ice wavenumber
         guess = kice

         !!check attenuation coeff's, wavenumbers
         if (w.eq.1) then
            print*,'***********************************************'
            print*,'check outputs from RTparam_outer:'
            print*,'T,h',T,iceh
            print*,'atten',alp_nd,damping
            print*,'ki,kw,2pi/wlng_wtr/guess',kice,kwtr,kw_inf,guess
            print*,'|T|,argRT_s',modT,argR,argT,int_adm
            print*,'***********************************************'
            print*,' '
         end if

         print*,'T (s), atten (per floe), damping (/m):'                &
     &         ,T,alp_nd,damping
      end do

      end program test_atten
