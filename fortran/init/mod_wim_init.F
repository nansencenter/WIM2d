      module mod_wim_init

      use mod_gridinfo
      use mod_common_wim

      implicit none

      real,parameter :: Hs_inc   = 2.0
      real,parameter :: Tp_inc   = 12.0
      real           :: mwd_inc                       !mwd - set in waves_init

      contains

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine waves_init(GRID_OPT)

      use mod_common_wim

      implicit none

      integer,intent(in)   :: GRID_OPT
      integer  :: i,j,w,wth
      real     :: xm,ym
      logical  :: critter

      integer,parameter :: nw    = n_wavdir
      integer,parameter :: ndir  = n_wave_freq

      real,parameter :: Tmin  = 2.5
      real,parameter :: Tmax  = 25.0
      real           :: fmin,fmax,df
      real           :: theta_min,theta_max,dtheta
      real           :: PI

      PI    = 4.0*atan(1.0)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!set frequency grid
      if (nw.eq.1) then
         freq_vec(1) = 1.0/Tp_inc
      else
         fmin  = 1/Tmax
         fmax  = 1/Tmin
         df    = (fmax-fmin)/(n_wave_freq-1.0)
         do w=1,nw
            freq_vec(w) = fmin+(w-1)*df
         end do
      end if

      !!set directional grid
      if (GRID_OPT.eq.0) then
         mwd_inc  = 135.0!!waves-from direction
      elseif (GRID_OPT.eq.1) then
         mwd_inc  = -90.0!!waves-from direction
      elseif (GRID_OPT.eq.2) then
         mwd_inc  = 135.0!!waves-from direction
      end if

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
      !! Set WAVEMASK

      xm = maxval(X)
      ym = maxval(Y)
      !!
      WAVE_MASK   = 0.0
      Hs          = 0.0
      Tp          = 0.0
      mwd         = 0.0

      if (GRID_OPT.eq.0) then

         !!island in corner, with some ice around it
         where ((X.lt.-0.4*xm).or.(Y.gt.0.4*ym))
            WAVE_MASK   = 1.0
         end where

      elseif (GRID_OPT.eq.1) then

         !!column of land on right
         where (X.lt.-.8*xm)
            WAVE_MASK   = 1.0
         end where

      elseif (GRID_OPT.eq.2) then

         !!land in southeast corner
         where ((X.lt.-0.25*xm).or.(Y.gt.0.25*ym))
            WAVE_MASK   = 1.0
         end where

      end if

      where (WAVE_MASK.eq.1.0)
         Hs    = Hs_inc
         Tp    = Tp_inc
         mwd   = mwd_inc
      end where
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      end subroutine waves_init
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine set_inc_waves()

      use mod_common_wim

      implicit none

      integer,parameter :: nw    = n_wavdir
      integer,parameter :: ndir  = n_wave_freq
      integer           :: i,j,w,wth
      !!
      real,dimension(nw)   :: Sfreq(nw)
      real,dimension(ndir) :: theta_fac
      real  :: PI,om,chi

      PI       = 4.0*atan(1.0)
      sdf_inc  = 0.0
      sdf_dir  = 0.0

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
            end do
         else
            Sfreq(1) = (Hs(i,j)/4.0)**2!!Hs=4*sqrt(Sfreq*wt_om), wt_om=1
         end if

         !!direcional spreading
         if (ndir.eq.1) then
            theta_fac   = 1.0
         else
            theta_fac   = 0.0
            do wth=1,ndir
               chi   = PI/180.0*(wavdir(wth)-mwd(i,j))
               if (cos(chi).gt.0.0) then
                  theta_fac(wth) = 2.0/PI*cos(chi)**2
               end if
            end do
         end if

         !!combine freq & dir
         do w=1,nw
         do wth=1,ndir
            sdf_inc(i,j,wth,w)   = Sfreq(w)*theta_fac(wth)
         end do
         end do

      end if
      end do
      end do

      end subroutine set_inc_waves
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine ice_init(GRID_OPT)

      use mod_common_wim
      use mod_gridinfo

      implicit none

      integer, intent(in)     :: GRID_OPT
      real,dimension(ii,jj)   :: Rsq
      real     :: x0,y0,NL,R0,xm
      
      WTR_MASK = 0.0
      if (GRID_OPT.eq.0) then

         !!island in corner, with some ice around it
         NL = sum(LANDMASK)        !!number of land cells
         x0 = sum(LANDMASK*X)/NL   !!mean x value of land cells
         y0 = sum(LANDMASK*Y)/NL   !!mean y value of land cells
         !!
         R0    = 70e3       !!radius
         Rsq   = (X-x0)**2+(Y-y0)**2
         where (Rsq.ge.R0**2)
            WTR_MASK = 1.0
         end where

      elseif (GRID_OPT.eq.1) then

         !!column of land on right
         xm       = maxval(X)
         where (X.lt.-.7*xm)
            WTR_MASK = 1.0
         end where

      elseif (GRID_OPT.eq.2) then

         !!land in southeast corner
         where ((X.lt.0.0).or.(Y.gt.0.0))
            WTR_MASK = 1.0
         end where

      end if

      ICE_MASK = (1-WTR_MASK)*(1-LANDMASK)!!0 on land & water
      icec     = unif_c*ICE_MASK
      iceh     = unif_h*ICE_MASK
      dfloe    = dfloe_pack_init*ICE_MASK

      end subroutine ice_init
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine SDF_Bretschneider(S,T,T_peak,H_sig)

      implicit none

      real,intent(out)  :: S
      real,intent(in)   :: T,T_peak,H_sig
      real              :: PI,omega,om_m,f1,f2,f3

      PI    = 4.0*atan(1.0)
      omega = 2*PI/T
      om_m  = 2*PI/T_peak
      !!
      f1 = 5.0/16.0*H_sig**2*om_m**4
      f2 = 1.0/omega**4
      f3 = exp(-1.25*(T/T_peak)**4)
      !!
      S  = f1*f2*f3

      end subroutine SDF_Bretschneider
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module mod_wim_init
