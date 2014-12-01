      module mod_wim_init

      use mod_gridinfo
      use mod_common_wim

      implicit none

      real,parameter :: Hs_inc   = 2.0
      real,parameter :: Tp_inc   = 12.0
      real           :: mwd_inc                       !mwd - set in waves_init

      contains

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine waves_init()

      use mod_common_wim

      implicit none

      integer  :: i,j,w,wth
      real     :: xm,ym
      logical  :: critter

      integer,parameter :: nw    = n_wave_freq
      integer,parameter :: ndir  = n_wavdir

      real,parameter :: Tmin  = 2.5
      real,parameter :: Tmax  = 25.0
      real           :: fmin,fmax,df
      real           :: theta_min,theta_max,dtheta
      real           :: PI

      PI    = 4.0*atan(1.0)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!set frequency grid
      print*,'***************************************'
      print*,'In waves_init,n_wave_freq=',nw
      if (nw.eq.1) then
         freq_vec(1) = 1.0/Tp_inc
         print*,'freq,T',freq_vec(1),1.0/freq_vec(1)
      else
         fmin  = 1/Tmax
         fmax  = 1/Tmin
         df    = (fmax-fmin)/(n_wave_freq-1.0)
         print*,'T range',Tmin,Tmax
         print*,'freq range',fmin,fmax
         do w=1,nw
            freq_vec(w) = fmin+(w-1)*df
            !print*,freq_vec(w)
         end do
      end if
      print*,'***************************************'
      print*,' '

      !!set directional grid
      mwd_inc  = -90.0!!waves-from direction

      if (ndir.eq.1) then
         wavdir(1)   = mwd_inc
      else
         theta_max   = 90.0
         theta_min   = -270.0
         dtheta      = (theta_min-theta_max)/(1.0*ndir)
         do wth=1,ndir
            wavdir(wth) = theta_max+(wth-1)*dtheta
         end do
      end if
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Set WAVEMASK & wave info
      !! TODO: enable reading this in from file?

      xm = maxval(X)
      ym = maxval(Y)
      !!
      WAVE_MASK   = 0.0
      Hs          = 0.0
      Tp          = 0.0
      mwd         = 0.0

      !!column of land on right
      where (X.lt.-.8*xm)
         WAVE_MASK   = 1.0
         Hs          = Hs_inc
         Tp          = Tp_inc
         mwd         = mwd_inc
      end where
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      end subroutine waves_init
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine set_inc_waves()

      use mod_common_wim

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
      subroutine ice_init()

      use mod_common_wim
      use mod_gridinfo

      implicit none

      real,dimension(ii,jj)   :: Rsq
      real     :: x0,y0,NL,R0,xm
      
      WTR_MASK = 0.0
      !!column of land on right
      xm       = maxval(X)
      where (X.lt.-.7*xm)
         WTR_MASK = 1.0
      end where

      ICE_MASK = (1-WTR_MASK)*(1-LANDMASK)!!0 on land & water
      icec     = unif_c*ICE_MASK
      iceh     = unif_h*ICE_MASK
      dfloe    = dfloe_pack_init*ICE_MASK

      end subroutine ice_init
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
