      !! test_waveadv_weno.F
      !! Author: Timothy Williams
      !! Date:   20140821, 12:22:17 CEST
      program main

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

      implicit none

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!variables
      integer,parameter :: SOLVER   = 1
      integer,parameter :: GRID_OPT = 1!option for land/ice/water config
      !!
      real     :: PI
      real     :: x_ext,y_ext,u_ref
      integer  :: nt,n,w,wth

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
      character(len=*),parameter :: log_file = "log/wim2d.log"
      character(len=*),parameter :: outdir   = "out/"

      real :: et,etime,tictoc(2) !!print time
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      PI       = 4.0*atan(1.0)
      et       = etime(tictoc)!!startup etime

      print*,'****************************************'
      print*,'Initialising...'
      print*,'****************************************'
      print*,' '

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!initialisation of grid, ice con's
      call set_wim_prams()
      call get_grid(GRID_OPT)
      call ice_init(GRID_OPT)

      !!incident waves
      call waves_init(GRID_OPT)
      call set_inc_waves()

      !! initialise WIM (mod_wavesice)
      !! set the attenuation coefficients;
      call wim_init(GRID_OPT)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (DO_CHECK_INIT.eq.1) then
         call check_grid(GRID_OPT)
         call check_init(GRID_OPT)
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

      print*,'Beginning main integration...'

      !nt = 13!stop at a certain point for testing
      do n=2,nt
         print*,'n/nt',n,nt
         call wim_step(SOLVER)

         if (DO_CHECK_PROG.eq.1) then
            !!dump fields to binary
            call check_prog(SOLVER,GRID_OPT,n)
         end if

      end do!! n - finish time stepping

      !!print info again
      include "print_screen_info.h"

      MIZ_MASK = 0.0
      where ((dfloe.gt.0.0).and.(dfloe.lt.dfloe_pack_init))
         MIZ_MASK = 1.0
      end where

      !!Dmax range in MIZ
      Dmax_max = maxval(MIZ_MASK*dfloe)
      tmp1     = dfloe+1.0e4*(1.0-MIZ_MASK)
      Dmax_min = minval(tmp1)

      W_MIZ = 0.0
      if (GRID_OPT.eq.1) then
         nMIZ  = sum(MIZ_MASK(:,1))
         W_MIZ = nMIZ*dx/1.0e3
         !!
         print*,' '
         print*,'MIZ width (km)',W_MIZ
      end if
      !!
      print*,'Dmax range in MIZ (m): ',Dmax_min,Dmax_max

      !!range of stresses
      taux_max = maxval(tau_x)
      taux_min = minval(tau_x)
      tauy_max = maxval(tau_y)
      tauy_min = minval(tau_y)
      print*,'tau_x range (Pa)',taux_min,taux_max
      print*,'tau_y range (Pa)',tauy_min,tauy_max
      print*,' '

      et = etime(tictoc)!!finish etime
      print*,' '
      print*,'******************************************************'
      print*,'p_WIM2d.F finished'
      !print *,'Elapsed time (s):',et                                    !&
      print *,'Elapsed time (min):',et/60.0                             !&
!    &       ,', user:', tictoc(1),                                     &
!    &       ,', sys:' ,tictoc(2)
      print*,'******************************************************'
      print*,' '

      !!print  some outputs to log file
      include "print_log_info.h"

      end program main
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
