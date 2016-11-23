      module mod_wim_init

      use mod_gridinfo
      use mod_common_wim

      implicit none

      contains

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine set_wim_prams(BRK_OPT,YOUNG_OPT)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      implicit none

      integer,intent(in)   :: BRK_OPT
      !real,optional,intent(in)   :: wim_prams_in(2)

      integer,optional,intent(in)   :: YOUNG_OPT

      real  :: vb,salt,temp
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!Brine vol fraction -> Young's modulus and flexural strength
      if (.true.) then
         !!just set
         wim_info%vbf   = .1        !! brine volume fraction 
         vb    = wim_info%vbf*1.0   !! brine volume [ppt]
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
         wim_info%vbf   = vb/1.0d3                         ! brine volume fraction 
      end if

      if (present(YOUNG_OPT)) then
         !!Robinson-Palmer viscosity parameter
         wim_info%drag_rp  = 13.

         !!Young's modulus
         if (YOUNG_OPT.eq.1) then
            !!just set Young's modulus
            wim_info%young = 2.0e9     !!lower ~ Marchenko
            !wim_info%young = 5.49e9    !!higher ~ Vernon
         else
            wim_info%young = 10.0e9*(1-3.51*wim_info%vbf)-1.0e9   ! Young's modulus [Pa]
         end if
      end if

      !!Flexural strength (Timco and O'Brien 1994)
      wim_info%sigma_c  = 1.76d6*exp(-5.88*sqrt(wim_info%vbf)) ! [Pa]

      !!Breaking criterion option
      if ((BRK_OPT.eq.0).or.(BRK_OPT.eq.1)) then !!beam test
         !! not used if BRK_OPT==0 (no breaking)
         wim_info%epsc     = wim_info%sigma_c/wim_info%young
         call fn_cohesion(1)

      elseif (BRK_OPT.eq.2) then!!Marchenko's stress criterion (plate)
         !! - convert to strain criterion
         wim_info%stress_c = 2.6*wim_info%sigma_c!!=E/(1-nu^2)*strain_c = thin plate (plane stress: \sigma_33=0)
         call fn_cohesion(BRK_OPT)

      elseif (BRK_OPT.eq.3) then!!Mohr-Coulomb criterion (plate)
         !! - convert to strain criterion
         call fn_cohesion(BRK_OPT)

      else
         print*,'Unknown value for BRK_OPT: ',BRK_OPT
         stop
      end if

      !! flex rigidity = flex_rig_coeff*h^3
      wim_info%flex_rig_coeff = wim_info%young
     &                           /12.0/(1-wim_info%poisson**2)

      print*,'********************************************'
      print*,'in set_wim_prams()...'
      print*,'Brine volume fraction:      ',wim_info%vbf
      print*,'Youngs modulus (Pa):        ',wim_info%young
      print*,'BRK_OPT:                    ',BRK_OPT
      print*,'Flexural strength (Pa):     ',wim_info%sigma_c
      print*,'Breaking stress (Pa):       ',wim_info%stress_c
      print*,'Breaking strain:            ',wim_info%epsc
      print*,'RP drag (Pa.s/m):           ',wim_info%drag_rp
      print*,'WS viscosity (m^2/s):       ',wim_info%visc_ws
      print*,'********************************************'
      print*,' '

      end subroutine set_wim_prams
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine fn_cohesion(BRK_OPT)

      !!make epsc,stress_c and cohesion consistent
      implicit none
      integer, intent(in)  :: BRK_OPT
      real                 :: alpha,sig_N

      alpha    = (1+wim_info%poisson)/(1-wim_info%poisson)
      if (BRK_OPT.eq.1) then
         wim_info%stress_c = wim_info%young
     &                        *wim_info%epsc/(1-wim_info%poisson**2)
         sig_N             = (1-wim_info%poisson)/2.*wim_info%stress_c
         wim_info%cohesion = (alpha+wim_info%friction)*sig_N
         
      elseif (BRK_OPT.eq.2) then
         wim_info%epsc     = (1-wim_info%poisson**2)
     &                           *wim_info%stress_c/wim_info%young
         sig_N             = (1-wim_info%poisson)/2.*wim_info%stress_c
         wim_info%cohesion = (alpha+wim_info%friction)*sig_N

      elseif (BRK_OPT.eq.3) then
         sig_N    = wim_info%cohesion/(alpha+wim_info%friction)
            !compressive stress  = .5*(sig_11+sig_22) = .5*(1-poisson)*sig11 at breaking point
         wim_info%stress_c = 2*sig_N/(1-wim_info%poisson)!!=E/(1-nu^2)*strain_c = thin plate (plane stress: \sigma_33=0)
         wim_info%epsc     = (1-wim_info%poisson**2)
     &               *wim_info%stress_c/wim_info%young
      end if

      end subroutine fn_cohesion


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine wavecons_ideal(init_cons)

      !use mod_common_wim

      implicit none
      real,intent(in)   :: init_cons(3)

      integer  :: i,j
      real     :: xm,x0,x_edge
      logical  :: critter

      !! wave parameters for ideal cases
      real :: Hs_inc,Tp_inc,mwd_inc

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Set WAVEMASK & wave info
      !! TODO: enable reading this in from file?

      !! wave parameters for ideal cases
      Hs_inc   = init_cons(1)    !! significant wave heig
      Tp_inc   = init_cons(2)    !! mean period
      mwd_inc  = init_cons(3)   !! waves-from direction )
      print*,'Initial Hs_inc,Tp_inc,mwd_inc'
      print*,Hs_inc,Tp_inc,mwd_inc
      print*,' '

      !! non-IO version of code
      !! - setting wave mask, Hs, Tp, mwd here
      print*,'***************************************'
      print*,'In waves_init:'
      print*,'setting WAVE_MASK, Hs, Tp, mwd'
      print*,'***************************************'
      print*,' '

      xm       = maxval(X)
      x0       = minval(X)
      x_edge  = .5*(x0+xm)-.8*(.5*(xm-x0))
      !!
      WAVE_MASK   = 0.0
      Hs          = 0.0
      Tp          = 0.0
      mwd         = 0.0

      !!strip where waves are initially defined on left
      where (X.lt.x_edge)
         WAVE_MASK   = 1.0
         Hs          = Hs_inc
         Tp          = Tp_inc
         mwd         = mwd_inc
      end where

      end subroutine wavecons_ideal
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine waves_init_freq()

      implicit none

      integer  :: w

      integer,parameter :: nw    = n_wave_freq

      real  :: fmin,fmax,df
      real  :: Tmean,Tvar,Tstdev

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!set frequency grid
      print*,'***************************************'
      print*,'In waves_init_freq,n_wave_freq=',nw

      if (nw.gt.1) then

         !! multiple frequencies
         !! - fixed grid
         fmin  = 1/Tmax
         fmax  = 1/Tmin
         df    = (fmax-fmin)/(n_wave_freq-1.0)
         print*,'T range',Tmin,Tmax
         print*,'freq range',fmin,fmax
         do w=1,nw
            freq_vec(w) = fmin+(w-1)*df
            !print*,freq_vec(w)
         end do

      else

         !! single frequency
         !! - initialise from Tp array
         !! - set it to mean over the wave mask

         !print*,sum(WAVE_MASK)
         Tmean       = sum(WAVE_MASK*Tp)/sum(WAVE_MASK)!mean period
         freq_vec(1) = 1.0/Tmean
         !!
         Tvar     = sum(WAVE_MASK*(Tp-Tmean)**2)/sum(WAVE_MASK)!variance of the periods
         Tstdev   = sqrt(Tvar)
         print*,'mean of peak periods : ',Tmean
         print*,'std dev of peak periods : ',Tstdev
         !!
         if (abs(Tstdev).gt.1.0e-2) then
            print*,'***********************************************'
            print*,'warning: trying to run multiple frequencies'
            print*,'  with n_wave_freq = 1'
            print*,'***********************************************'
            print*,' '
            STOP
         end if

      end if

      print*,'***************************************'
      print*,' '

      end subroutine waves_init_freq
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine waves_init_dir()

      implicit none

      integer  :: wth

      integer,parameter :: ndir  = n_wavdir

      real  :: theta_min,theta_max,dtheta
      real  :: PI
      real  :: Tmean,Tvar,Tstdev

      PI    = 4.0*atan(1.0)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!set directional grid
      print*,'***************************************'
      print*,'In waves_init_dir, n_wavdir=',ndir


      !!set directional grid
      if (ndir.gt.1) then

         !! fixed directional grid

         theta_max   = 90.0
         theta_min   = -270.0
         dtheta      = (theta_min-theta_max)/(1.0*ndir)
         do wth=1,ndir
            wavdir(wth) = theta_max+(wth-1)*dtheta
         end do

      else

         !! single direction:
         !! - set it to mean over the wave mask
         wavdir(1)   = sum(WAVE_MASK*mwd)/sum(WAVE_MASK)!mean direction over wave mask
         !!
         Tvar     = sum(WAVE_MASK*(mwd-wavdir(1))**2)/sum(WAVE_MASK)!variance of the directions
         Tstdev   = sqrt(Tvar)
         print*,'std dev of wave directions : ',Tstdev
         !!
         if (abs(Tstdev).gt.1.0e-2) then
            print*,'***********************************************'
            print*,'warning: trying to run multiple directions'
            print*,'  with n_wavdir = 1'
            print*,'***********************************************'
            print*,' '
            STOP
         end if

      end if

      print*,'***************************************'
      print*,' '

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end subroutine waves_init_dir
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine set_inc_waves()

      implicit none

      integer,parameter :: nw    = n_wave_freq
      integer,parameter :: ndir  = n_wavdir
      integer           :: i,j,w,wth
      !!
      real,dimension(nw)   :: Sfreq(nw)
      real,dimension(ndir) :: theta_fac
      real  :: PI,om,chi,adv_dir,dtheta

      PI       = 4.0*atan(1.0)
      sdf_inc  = 0.0

      do j=1,jj
      do i=1,ii
      if (WAVE_MASK(i,j).eq.1.0) then
         !!set the incident waves

         !!freq
         if (nw.gt.1) then
            !! Bretschneider spectrum
            !! - freq spectrum
            do w=1,nw
               om = 2*PI*freq_vec(w)
               call SDF_Bretschneider(                                  &
     &                  Sfreq(w),om,Tp(i,j),Hs(i,j))
!              if (i*j.eq.1) then
!                 print*,'T,S_B',2*PI/om,Sfreq(w)
!              end if
            end do
         else
            Sfreq(1) = (Hs(i,j)/4.0)**2!!Hs=4*sqrt(Sfreq*wt_om), wt_om=1
!           if (i*j.eq.1) then
!              print*,'T,S_inc',1.0/freq_vec(1),Sfreq(1)
!           end if
         end if

         !!direcional spreading
         if (ndir.eq.1) then
            theta_fac   = 1.0
         else
            !theta_fac   = 0.0
            dtheta   = abs(wavdir(2)-wavdir(1)) !degrees
            !print*,dtheta
            do wth=1,ndir
#if 1
               theta_fac(wth) = theta_dirfrac(                          &
     &            wavdir(wth)-dtheta/2.,dtheta,mwd(i,j))                &
     &               *180./(pi*dtheta)
#else
               chi   = PI/180.0*(wavdir(wth)-mwd(i,j))
               if (cos(chi).gt.0.0) then
                  theta_fac(wth) = 2.0/PI*cos(chi)**2
               else
                  theta_fac(wth) = 0.0
               end if

#endif
               !print*,wavdir(wth),mwd(i,j),theta_fac(wth)

            end do
            !print*,'test dir spreading factor'
            !print*,2*PI/(1.0*ndir)*sum(theta_fac),1
         end if

         !!combine freq & dir
         do w=1,nw
         do wth=1,ndir
            sdf_inc(i,j,wth,w)   = Sfreq(w)*theta_fac(wth)
            if (i*j.eq.1) then
               adv_dir  = -PI/180.0*(wavdir(wth)+90.0)
!              print*,adv_dir,sdf_inc(i,j,wth,w),Sfreq(w),theta_fac(wth)
            end if
         end do
         end do

      end if!WAVE_MASK==1
      end do!j - spatial loop (x)
      end do!j - spatial loop (y)

      end subroutine set_inc_waves
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine icecons_ideal(init_cons)

      implicit none

      real,intent(in)   :: init_cons(3)
      real           :: x0,xm,x_edge
      real  :: unif_c,unif_h,Dmax_init

      unif_c      = init_cons(1) !initial const conc
      unif_h      = init_cons(2) !initial const thickness
      Dmax_init   = init_cons(3) !initial floe size
      print*,'Initial c,h,Dmax'
      print*,unif_c,unif_h,Dmax_init
      print*,' '

      
      WTR_MASK = 0.0
      !! strip of water on left
      !! - NB waves not everywhere yet
      !! - see WAVE_MASK definition
      xm       = maxval(X)
      x0       = minval(X)
      x_edge  = .5*(x0+xm)-.7*(.5*(xm-x0))
      where (X.lt.x_edge)
         WTR_MASK = 1.0
      end where

      ICE_MASK = (1-WTR_MASK)*(1-LANDMASK)!!0 on land & water
      icec     = unif_c*ICE_MASK
      iceh     = unif_h*ICE_MASK
      dfloe    = Dmax_init*ICE_MASK

      end subroutine icecons_ideal
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine SDF_Bretschneider(S,omega,T_peak,H_sig)

      implicit none

      real,intent(out)  :: S
      real,intent(in)   :: omega,T_peak,H_sig
      real              :: PI,T,om_m,f1,f2,f3

      PI    = 4.0*atan(1.0)
      T     = 2*PI/omega
      om_m  = 2*PI/T_peak
      !!
      f1 = 5.0/16.0*H_sig**2*om_m**4
      f2 = 1.0/omega**5
      f3 = exp(-1.25*(T/T_peak)**4)
      !!
      S  = f1*f2*f3

      end subroutine SDF_Bretschneider
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%!!!!!!!!!!!!!!


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real function theta_dirfrac(th1_,dtheta_,mwd_)
      !! chi=pi/180*(theta-mwd) 
      !!  (convert to radians and centre about mwd, the mean wave direction)
      !! chi1=pi/180*(th1-mwd)
      !! chi2=pi/180*(th2-mwd)
      !! D(chi)=2/pi*cos^2(chi) if chi\in[-pi/2,pi/2] 
      !! D(chi)=0 if chi\in[-pi/2,pi/2] (no backwards waves)
      !! theta_dirfrac  = \int_{chi1}^{chi2}D(chi)dchi
      !!                = 1 if -chi1=chi2=pi/2
      !!                = the fraction of the energy going in the directions in [chi1,chi2]

      implicit none
      
      real, intent(in) :: th1_,dtheta_,mwd_
      real             :: th2_,integral
      real             :: chi1,chi2,phi1,phi2,L1,L2
      real, parameter  :: pi  = 3.1415927

      !!get mwd\pm90 inside [th1_,th1_+360)
      phi1  = theta_in_range(mwd_-90,th1_)!>=th1_
      phi2  = theta_in_range(mwd_+90,th1_)!>=th1_
      th2_  = th1_+dtheta_

      integral     = 0.
      if (phi2.gt.phi1) then
         !! th1_,phi1,phi2, and th2_
         L1       = max(th1_,phi1)
         L2       = min(th2_,phi2)
         L2       = max(L1,L2)         !make L2>=L1
         chi1     = pi*(L1-mwd_)/180.
         chi2     = pi*(L2-mwd_)/180.
         integral = integral
     &      + 2*(chi2-chi1)+sin(2*chi2)-sin(2*chi1)
      else
         !! th1_,phi2,phi1, and th2_
         !! 1st consider (th1_,phi2) interval
         L1       = th1_
         L2       = min(th2_,phi2)
         chi1     = pi*(L1-mwd_)/180.
         chi2     = pi*(L2-mwd_)/180.
         integral = integral
     &      + 2*(chi2-chi1)+sin(2*chi2)-sin(2*chi1)
      
         !! 2nd consider (phi1,th2_) interval
         L1       = phi1
         L2       = max(L1,th2_)       !make L2>=L1
         chi1     = pi*(L1-mwd_)/180.
         chi2     = pi*(L2-mwd_)/180.
         integral = integral
     &      + 2*(chi2-chi1)+sin(2*chi2)-sin(2*chi1)
      end if

      theta_dirfrac  = integral/(2.*pi)
      end function theta_dirfrac
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real function theta_in_range(th_,th1)
      !! th = theta_in_range(th_)
      !! where th1<th<=th1+360.

      implicit none
      
      real, intent(in)  :: th_
      real, parameter   :: pi  = 3.1415927
      real              :: th,th1,th2,dth
      integer           :: njump

      th2   = th1+360.
      if (th_.lt.th1) then
         dth   = th1-th_
         njump = ceiling(dth/360.)
         th    = th_+njump*360.
      elseif (th_.gt.th2) then
         dth   = th_-th2
         njump = ceiling(dth/360.)
         th    = th_-njump*360.
      elseif (th_.eq.th2) then
         th = th1
      else
         th = th_
      end if

      theta_in_range = th
      end function
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module mod_wim_init
