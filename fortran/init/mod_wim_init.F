      module mod_wim_init

      use mod_gridinfo
      use mod_common_wim

      implicit none

      contains

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine wavecons_ideal()

      !use mod_common_wim

      implicit none

      integer  :: i,j
      real     :: xm,x0,x_edge
      logical  :: critter

      !! wave parameters for ideal cases
      real,parameter :: Hs_inc   = 2.0    !! significant wave height
      real,parameter :: Tp_inc   = 12.0   !! mean period
      real           :: mwd_inc  = -90.   !! waves-from direction

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Set WAVEMASK & wave info
      !! TODO: enable reading this in from file?


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

         Tmean       = sum(WAVE_MASK*Tp)/sum(WAVE_MASK)!mean period
         freq_vec(1) = 1.0/Tmean
         !!
         Tvar     = sum(WAVE_MASK*(Tp-Tmean)**2)/sum(WAVE_MASK)!variance of the periods
         Tstdev   = sqrt(Tvar)
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
      real  :: PI,om,chi,adv_dir

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
            do wth=1,ndir
               chi   = PI/180.0*(wavdir(wth)-mwd(i,j))
               if (cos(chi).gt.0.0) then
                  theta_fac(wth) = 2.0/PI*cos(chi)**2
               else
                  theta_fac(wth) = 0.0
               end if
            end do
!           print*,'test dir spreading factor'
!           print*,2*PI/(1.0*ndir)*sum(theta_fac),1
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
      subroutine icecons_ideal()

      implicit none

      real           :: x0,xm,x_edge
      real,parameter :: unif_c   = .7
      real,parameter :: unif_h   = 2.0
      
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
      dfloe    = dfloe_pack_init*ICE_MASK

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
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module mod_wim_init
