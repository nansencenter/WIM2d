      !mod_WIM2d_run.F
      !main code
      !Author: Timothy Williams
      !Date:   20141118, 17:20:34 CET

      module mod_wim2d_run

      implicit none
      
      contains
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine wim2d_run

      use mod_gridinfo     !! grid info
      use mod_wim_prams    !! ice parameters
      use mod_wim_init     !! initial wave fields
      use mod_common_wim   !! initial ice fields
      use mod_adv_atten    !! advection routine

      !! main code for WIM
      use mod_wavesice, only:
     &    amax,amin,CFL,dt    !need these to determine number of time steps
     &   ,wim_init,wim_step   !main subroutines
     &   ,check_grid,check_init,check_prog   !subroutines for dumping fields

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! import main code for WIM
      use mod_wavesice, only:
     &    amax,amin,CFL,dt                                              &
     &   ,wim_init,wim_step                                             &
     &   ,check_grid,check_init,check_prog
            !! need these to:
            !! - determine number of time steps
            !! - get main subroutines
            !! - get subroutines for dumping fields
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!variables
      integer,parameter :: SOLVER   = 1
      integer,parameter :: ADV_DIM  = 2
      !!
      real     :: PI
      real     :: x_ext,y_ext,u_ref
      integer  :: nt,n,w,wth,reps
      logical  :: critter

      integer,parameter :: DO_CHECK_PROG  = 1
      integer,parameter :: DO_CHECK_INIT  = 1
      integer,parameter :: ndir  = n_wavdir
      integer,parameter :: nw    = n_wave_freq

      !!test outputs
      real,dimension(idm,jdm) :: MIZ_MASK,tmp1
      real                    :: W_MIZ,taux_max,tauy_max                &
     &   ,taux_min,tauy_min,Dmax_min,Dmax_max
      integer                 :: nMIZ

      !!
      character(len=80)          :: outdir
      character(len=80)          :: log_file
      character(len=*),parameter :: this_subr   = 'wim2d_run'

      real :: et,etime,tictoc(2) !!print time
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      PI       = 4.0*atan(1.0)
      et       = etime(tictoc)!!startup etime
      outdir   = 'out'
      log_file = trim(outdir)//"/log/wim2d.log"

      print*,'****************************************'
      print*,'Initialising...'
      print*,'****************************************'
      print*,' '

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!initialisation of grid, ice con's
      call set_wim_prams() ! - ice properties like Young's modulus
      call get_grid_file() ! - grid: X,Y,LAND_MASK
      call ice_init()      ! - initial icec,iceh,dfloe

      !!incident waves
      call waves_init()    ! - initial Tp,Hs,mwd
      call set_inc_waves() ! - fill in the incident wave spectrum
                           !   with a Bretschneider spectrum

      !! initialise WIM (mod_wavesice)
      !! set the attenuation coefficients;
      call wim_init()
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (DO_CHECK_INIT.eq.1) then
         call check_grid(outdir)
         call check_init(outdir)
      end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! fill wave spectra arrays using (Hs,Tp,mwd)
      !!
      !! could just set sdf_dir = sdf_inc,
      !! but this way lets us retain any waves
      !! that are outside the incident wave mask
      sdf_dir  = 0.0
      do w=1,nw
      do wth=1,ndir
         where (WAVE_MASK.eq.1.0)
            sdf_dir(:,:,wth,w)   = sdf_inc(:,:,wth,w)
         end where
      end do
      end do
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!set number of time steps
      u_ref = amin+.7*(amax-amin);
      dt    = CFL*dx/amax
      x_ext = idm*dx/1.0e3!km
      y_ext = jdm*dy/1.0e3!km
      nt    = nint( 1.0e3*x_ext/u_ref/dt )
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!print info
      include "print_screen_info.h"

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      print*,'Beginning main integration...'

      !! text that is common to both
      !!    wim2d_run 
      !!    wim2d_run_io
      !! subroutines:
      include "WIM2d_run.h"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!print info to screen again
      include "print_screen_info.h"

      !!print  some outputs to log file
      include "print_log_info.h"

      end subroutine wim2d_run
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine wim2d_run_io(in_arrays,out_arrays,int_prams,nx,ny)

      use mod_gridinfo     !! grid info
      use mod_wim_prams    !! ice parameters
      use mod_wim_init     !! initial wave fields
      use mod_common_wim   !! initial ice fields
      use mod_adv_atten    !! advection routine

      !! main code for WIM
      use mod_wavesice, only:
     &    amax,amin,CFL,dt    !need these to determine number of time steps
     &   ,wim_init,wim_step   !main subroutines
     &   ,check_grid,check_init,check_prog   !subroutines for dumping fields

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! import main code for WIM
      use mod_wavesice, only:
     &    amax,amin,CFL,dt                                              &
     &   ,wim_init,wim_step                                             &
     &   ,check_grid,check_init,check_prog
            !! need these to:
            !! - determine number of time steps
            !! - get main subroutines
            !! - get subroutines for dumping fields
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!variables
      integer,intent(in)                  :: nx,ny,int_prams(2)
      real,dimension(nx,ny,6),intent(in)  :: in_arrays!!c,h,Dmax,Hs,Tp,Dmax
      real,dimension(nx,ny,5),intent(out) :: out_arrays!!Dmax,Hs,Tp,taux,tauy
!f2py intent(hide)   :: nx,ny

      integer   :: SOLVER,ADV_DIM
      !!
      real     :: PI
      real     :: x_ext,y_ext,u_ref
      integer  :: nt,n,w,wth,reps
      logical  :: critter

      integer,parameter :: DO_CHECK_PROG  = 1
      integer,parameter :: DO_CHECK_INIT  = 1
      integer,parameter :: ndir  = n_wavdir
      integer,parameter :: nw    = n_wave_freq

      real,parameter    :: min_conc = 0.05

      !!test outputs
      real,dimension(idm,jdm) :: MIZ_MASK,tmp1
      real                    :: W_MIZ,taux_max,tauy_max                &
     &   ,taux_min,tauy_min,Dmax_min,Dmax_max                           &
     &   ,Tmean,Tvar,Tstdev
      integer                 :: nMIZ

      !!
      character(len=80)          :: outdir
      character(len=80)          :: log_file
      character(len=*),parameter :: this_subr   = 'wim2d_run_io'

      real :: et,etime,tictoc(2) !!print time
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SOLVER   = int_prams(1) ! solver to use for advection/attenuation
      ADV_DIM  = int_prams(2) ! 1d advection option can be useful for checking
                              ! and also testing many parameters at once
                              ! ie each 'y' row can be a different experiment
                              !    if we use this option
                              !    - in practice we set speed in 'y' dirn to 0.

      PI       = 4.0*atan(1.0)
      et       = etime(tictoc)!!startup etime
      outdir   = "out_io/"
      log_file = trim(outdir)//'/log/wim2d.log'

      print*,'****************************************'
      print*,'Initialising...'
      print*,'****************************************'
      print*,' '

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!initialisation of grid, ice con's
      call set_wim_prams() ! - ice properties like Young's modulus
      call get_grid_file() ! - grid: X,Y,LAND_MASK

      !!ice conditions from inputs
      icec     = in_arrays(:,:,1)
      iceh     = in_arrays(:,:,2)
      dfloe    = in_arrays(:,:,3)
      ICE_MASK = 0.0
      where (icec.gt.min_conc)
         ICE_MASK = 1.0
      end where

      call waves_init()

      !! override incident waves from inputs
      !! - still need waves_init to set freq_vec
      Hs          = in_arrays(:,:,4)
      Tp          = in_arrays(:,:,5)
      mwd         = in_arrays(:,:,6)
      WAVE_MASK   = 0.0
      where ((ICE_MASK.lt.1.0).and.(Hs.gt.0.0))
         !!kill waves inside ice and propagate them in ourselves
         WAVE_MASK   = 1.0
      end where

      if (nw.eq.1) then
         !!single frequency:
         !!-set it to mean over the wave mask
         Tmean       = sum(WAVE_MASK*Tp)/sum(WAVE_MASK)!mean period
         freq_vec(1) = 1.0/Tmean
         !!
         Tvar  = sum(WAVE_MASK*(Tp-Tmean)**2)/sum(WAVE_MASK)!variance of the periods
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

      call set_inc_waves()       ! - fill in the incident wave spectrum
                                 !   with a Bretschneider spectrum

      !! initialise WIM (mod_wavesice)
      !! set the attenuation coefficients;
      call wim_init()
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (DO_CHECK_INIT.eq.1) then
         call check_grid(outdir)
         call check_init(outdir)
      end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! fill wave spectra arrays using (Hs,Tp,mwd)
      !!
      !! could just set sdf_dir = sdf_inc,
      !! but this way lets us retain any waves
      !! that are outside the incident wave mask
      sdf_dir  = 0.0
      do w=1,nw
      do wth=1,ndir
         where (WAVE_MASK.eq.1.0)
            sdf_dir(:,:,wth,w)   = sdf_inc(:,:,wth,w)
         end where
      end do
      end do
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!set number of time steps
      u_ref = amin+.7*(amax-amin);
      dt    = CFL*dx/amax
      x_ext = idm*dx/1.0e3!km
      y_ext = jdm*dy/1.0e3!km
      nt    = nint( 1.0e3*x_ext/u_ref/dt )
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!print info
      include "print_screen_info.h"

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      print*,'Beginning main integration...'

      !! text that is common to both
      !!    wim2d_run 
      !!    wim2d_run_io
      !! subroutines:
      include "WIM2d_run.h"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!print info to screen again
      include "print_screen_info.h"

      !!print  some outputs to log file
      include "print_log_info.h"

      !! put in same order as binary outputs:
      out_arrays(:,:,1) = dfloe
      out_arrays(:,:,2) = tau_x
      out_arrays(:,:,3) = tau_y
      out_arrays(:,:,4) = Hs
      out_arrays(:,:,5) = Tp

      end subroutine wim2d_run_io
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module mod_wim2d_run
