      !mod_WIM2d_run.F
      !main code
      !Author: Timothy Williams
      !Date:   20141118, 17:20:34 CET

      module mod_wim2d_run

      implicit none

      private

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! variables common to subroutines here
      integer  :: SCATMOD,ADV_DIM,ADV_OPT,STEADY
      integer  :: DO_CHECK_PROG,DO_CHECK_INIT,DO_CHECK_FINAL
      integer  :: DO_BREAKING !can turn off breaking for testing of waves

      real     :: x_ext,y_ext,u_ref,duration
      integer  :: nt

      character(len=80) :: log_file,this_subr,outdir,indir
      character(len=*),parameter :: infile_dirs = 'infile_dirs.txt'
      character(len=*),parameter :: infile      = 'infile.txt'

      real  :: W_MIZ,taux_max,tauy_max                                  &
     &   ,taux_min,tauy_min,Dmax_min,Dmax_max                           &

      real :: et,tictoc(2) !!print time
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!only subroutines are public
      public   :: wim2d_run,wim2d_run_io,wim2d_run_vSdir                &
     &            ,mex_io_gate,mex_io_gate_v2,mex_io_gate_vSdir
      
      contains
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine wim2d_run

      use mod_gridinfo     !! grid info
      use mod_wim_prams    !! ice parameters
      use mod_wim_init     !! initial wave fields
      use mod_common_wim   !! initial ice fields
      use mod_adv_atten    !! advection routine

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
            !! - NB amax,amin,CFL,dt set in wim_init
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!variables
      !!
      real     :: PI

      logical,parameter :: SET_WAVES   = .true.
      integer,parameter :: ndir        = n_wavdir
      integer,parameter :: nw          = n_wave_freq
      integer           :: w,wth

      !! function to print time
      real     :: etime
      real     :: wim_prams_in(2)

      !! does infile* exist?
      logical  :: fex,fex2
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      PI       = 4.0*atan(1.0)
      et       = etime(tictoc)!!startup etime
      !!
      wim_prams_in   = 0.0

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!locations of inputs/outputs
      inquire(file=infile_dirs,exist=fex)
      if (fex) then
         call read_infile_dirs(indir,outdir)
      else
         indir    = 'inputs'
         outdir   = 'out'
      end if
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      log_file    = trim(outdir)//"/log/wim2d.log"
      this_subr   = 'wim2d_run'

      print*,'****************************************'
      print*,'Initialising...'
      print*,'****************************************'
      print*,' '


      inquire(file=infile,exist=fex2)
      if (fex2) then
         call read_infile(wim_prams_in,duration)

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!initialisation of grid, ice con's
         call set_wim_prams(wim_prams_in) ! - ice properties like Young's modulus
      else
         SCATMOD     = 1
         ADV_DIM     = 2
         ADV_OPT     = 2
            !!0, not periodic; 1, periodic in both x,y; 2, periodic in y only;
         DO_BREAKING = 1
         STEADY      = 1

         DO_CHECK_PROG  = 1
         DO_CHECK_INIT  = 1
         DO_CHECK_FINAL = 1

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!initialisation of grid, ice con's
         call set_wim_prams() ! - ice properties like Young's modulus

      end if
      call get_grid_file(indir) ! - grid: X,Y,LAND_MASK
      call ice_init()      ! - initial icec,iceh,dfloe

      !! set wave grid
      !! SET_WAVES = .true. so also set initial Tp,Hs,mwd,WAVE_MASK
      call waves_init(SET_WAVES)

      !! create full incident wave spectrum
      !! - use preset wave parameters (Tp,Hs,mwd,WAVE_MASK)
      !!   to fill in the incident wave spectrum
      !!   with a Bretschneider spectrum
      call set_inc_waves()

      !! initialise WIM (mod_wavesice)
      !! set the attenuation coefficients;
      call wim_init(STEADY)

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
      dt       = CFL*dx/amax
      x_ext    = idm*dx/1.0e3!km
      y_ext    = jdm*dy/1.0e3!km
      if (.not.fex2) then
         !!duration of simulation
         u_ref    = amin+.7*(amax-amin)
         duration = 1.0e3*x_ext/u_ref
      end if
      nt       = nint( duration/dt )
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!print info
      call print_screen_info()

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      print*,'Beginning main integration...'
      call wim2d_run_common()

      et = etime(tictoc)!!finish etime
      print*,' '
      print*,'******************************************************'
      print*,'run of 2d WIM2d finished'
      !print *,'Elapsed time (s):',et
      print *,'Elapsed time (min):',et/60.0                             !&
!    &       ,', user:', tictoc(1),                                     &
!    &       ,', sys:' ,tictoc(2)
      print*,'******************************************************'
      print*,' '
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!print info to screen again
      call print_screen_info()

      !!print  some outputs to log file
      call print_log_info()

      end subroutine wim2d_run
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine wim2d_run_io(in_arrays,out_arrays                      &
     &                       ,int_prams,real_prams,nx,ny)

      use mod_gridinfo     !! grid info
      use mod_wim_prams    !! ice parameters
      use mod_wim_init     !! initial wave fields
      use mod_common_wim   !! initial ice fields
      use mod_adv_atten    !! advection routine

!     !! main code for WIM
!     use mod_wavesice, only:
!    &    amax,amin,CFL,dt    !need these to determine number of time steps
!    &   ,wim_init,wim_step   !main subroutines
!    &   ,check_grid,check_init,check_prog   !subroutines for dumping fields

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
      integer,intent(in)                  :: nx,ny
      integer,dimension(8),intent(in)     :: int_prams
      ! SCATMOD,ADV_DIM,ADV_OPT,
      ! DO_CHECK_FINAL,DO_CHECK_PROG,DO_CHECK_INIT,
      ! DO_BREAKING,STEADY

      real,dimension(3),intent(in)        :: real_prams  !! young,visc_rp,duration
      real,dimension(nx,ny,6),intent(in)  :: in_arrays   !! c,h,Dmax,Hs,Tp,mwd
      real,dimension(nx,ny,5),intent(out) :: out_arrays  !! Dmax,Hs,Tp,taux,tauy

      real              :: PI
      real,dimension(2) :: wim_prams_in  !! (young,visc_rp) -> inputs to set_wim_prams

      logical,parameter :: SET_WAVES   = .false.
      integer,parameter :: ndir  = n_wavdir
      integer,parameter :: nw    = n_wave_freq
      integer           :: w,wth

      real,parameter    :: min_conc = 0.05

      real  :: Tmean,Tvar,Tstdev

      !! function to print time
      real  :: etime

      !! does infile_dirs exist?
      logical :: fex
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SCATMOD   = int_prams(1) ! solver to use for advection/attenuation
                               ! 0: scattering treated like attenuation
                               ! 1: directional scattering
      ADV_DIM  = int_prams(2)  ! 1d advection option can be useful for checking
                               ! and also testing many parameters at once
                               ! ie each 'y' row can be a different experiment
                               !    if we use this option
                               !    - in practice we set speed in 'y' dirn to 0.
      ADV_OPT  = int_prams(3)  ! 0, not periodic;
                               ! 1, periodic in both x,y;
                               ! 2, periodic in y only;

      !! Inputs to determine if results
      !! should be dumped to binaries or not
      DO_CHECK_FINAL = int_prams(4)
      DO_CHECK_PROG  = int_prams(5)
      DO_CHECK_INIT  = int_prams(6)
      DO_BREAKING    = int_prams(7)
      STEADY         = int_prams(8)

      !!real inputs
      wim_prams_in   = real_prams(1:2)    !!young,visc_rp
      duration       = real_prams(3)      !!length of sim [s]

      PI       = 4.0*atan(1.0)
      et       = etime(tictoc)!!startup etime

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!locations of inputs/outputs
      inquire(file=infile_dirs,exist=fex)
      if (fex) then
         call read_infile_dirs(indir,outdir)
      else
         indir    = 'inputs'
         outdir   = 'out_io'
      end if
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      log_file    = trim(outdir)//"/log/wim2d.log"
      this_subr   = 'wim2d_run_io'

      print*,'****************************************'
      print*,'Initialising...'
      print*,'****************************************'
      print*,' '

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!initialisation of grid, ice con's
      call set_wim_prams(wim_prams_in) ! - ice properties like Young's modulus
      call get_grid_file(indir) ! - grid: X,Y,LAND_MASK

      !!ice conditions from inputs
      icec     = in_arrays(:,:,1)
      iceh     = in_arrays(:,:,2)
      dfloe    = in_arrays(:,:,3)
      ICE_MASK = 0.0
      where (icec.gt.min_conc)
         ICE_MASK = 1.0
      end where

      !! set wave grid
      !! SET_WAVES = .false. so don't set initial Tp,Hs,mwd,WAVE_MASK as well
      call waves_init(SET_WAVES)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! set incident waves from inputs
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

      if (ndir.eq.1) then
         !!single direction:
         !!-set it to mean over the wave mask
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
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !! create full incident wave spectrum
      !! - use preset wave parameters (Tp,Hs,mwd,WAVE_MASK)
      !!   to fill in the incident wave spectrum
      !!   with a Bretschneider spectrum
      call set_inc_waves()

      !! initialise WIM (mod_wavesice)
      !! set the attenuation coefficients;
      call wim_init(STEADY)
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
      dt    = CFL*dx/amax
      x_ext = idm*dx/1.0e3!km
      y_ext = jdm*dy/1.0e3!km
      nt    = nint(duration/dt)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!print info
      call print_screen_info()

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      print*,'Beginning main integration...'
      call wim2d_run_common()

      et = etime(tictoc)!!finish etime
      print*,' '
      print*,'******************************************************'
      print*,'run of 2d WIM2d finished'
      !print *,'Elapsed time (s):',et
      print *,'Elapsed time (min):',et/60.0                             !&
!    &       ,', user:', tictoc(1),                                     &
!    &       ,', sys:' ,tictoc(2)
      print*,'******************************************************'
      print*,' '
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!print info to screen again
      call print_screen_info()

      !!print  some outputs to log file
      call print_log_info()

      !! put in same order as binary outputs:
      out_arrays(:,:,1) = dfloe     ! maximum floe size
      out_arrays(:,:,2) = tau_x     ! wave stress (x dirn)
      out_arrays(:,:,3) = tau_y     ! wave stress (y dirn)
      out_arrays(:,:,4) = Hs        ! significant wave height
      out_arrays(:,:,5) = Tp        ! mean wave period

      end subroutine wim2d_run_io
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine wim2d_run_vSdir(out_arrays,int_prams,real_prams)

      use mod_gridinfo     !! grid info
      use mod_wim_prams    !! ice parameters
      use mod_wim_init     !! initial wave fields
      use mod_common_wim   !! initial ice fields
      use mod_adv_atten    !! advection routine

!     !! main code for WIM
!     use mod_wavesice, only:
!    &    amax,amin,CFL,dt    !need these to determine number of time steps
!    &   ,wim_init,wim_step   !main subroutines
!    &   ,check_grid,check_init,check_prog   !subroutines for dumping fields

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! import main code for WIM
      use mod_wavesice, only:
     &    amax,amin,CFL,dt                                              &
     &   ,wim_init,wim_step                                             &
     &   ,check_grid,check_init,check_prog                              &
     &   ,do_wim_integrals
            !! need these to:
            !! - determine number of time steps
            !! - get main subroutines
            !! - get subroutines for dumping fields
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!variables
      integer,parameter                :: ndir  = n_wavdir
      integer,parameter                :: nw    = n_wave_freq
      integer,dimension(8),intent(in)  :: int_prams
      ! SCATMOD,ADV_DIM,ADV_OPT,
      ! DO_CHECK_FINAL,DO_CHECK_PROG,DO_CHECK_INIT,
      ! DO_BREAKING,STEADY

      real,dimension(ii,jj,5),intent(out) :: out_arrays        !! Dmax,Hs,Tp,taux,tauy
      real,dimension(3),intent(in)        :: real_prams        !! young,visc_rp,duration

      real              :: PI
      real,dimension(2) :: wim_prams_in  !! (young,visc_rp) -> inputs to set_wim_prams

      logical,parameter :: SET_WAVES   = .false.
      integer           :: w,wth

      real,parameter    :: min_conc = 0.05

      real  :: Tmean,Tvar,Tstdev

      !! function to print time
      real  :: etime

      !! does infile_dirs exist?
      logical :: fex
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SCATMOD   = int_prams(1) ! solver to use for advection/attenuation
                               ! 0: scattering treated like attenuation
                               ! 1: directional scattering
      ADV_DIM  = int_prams(2)  ! 1d advection option can be useful for checking
                               ! and also testing many parameters at once
                               ! ie each 'y' row can be a different experiment
                               !    if we use this option
                               !    - in practice we set speed in 'y' dirn to 0.
      ADV_OPT  = int_prams(3)  ! 0, not periodic;
                               ! 1, periodic in both x,y;
                               ! 2, periodic in y only;

      !! Inputs to determine if results
      !! should be dumped to binaries or not
      DO_CHECK_FINAL = int_prams(4)
      DO_CHECK_PROG  = int_prams(5)
      DO_CHECK_INIT  = int_prams(6)
      DO_BREAKING    = int_prams(7)
      STEADY         = int_prams(8)

      !!real inputs
      wim_prams_in   = real_prams(1:2)    !!young,visc_rp
      duration       = real_prams(3)      !!length of simulation [s]

      PI       = 4.0*atan(1.0)
      et       = etime(tictoc)!!startup etime

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!locations of inputs/outputs
      inquire(file=infile_dirs,exist=fex)
      if (fex) then
         call read_infile_dirs(indir,outdir)
      else
         indir    = 'inputs'
         outdir   = 'out_2'
      end if
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      log_file    = trim(outdir)//"/log/wim2d.log"
      this_subr   = 'wim2d_run_io'

      print*,'****************************************'
      print*,'Initialising...'
      print*,'****************************************'
      print*,' '

      !print*,freq_vec
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!initialisation of grid, ice con's
      call set_wim_prams(wim_prams_in) ! - ice properties like Young's modulus
      call get_grid_file(indir) ! - grid: X,Y,LAND_MASK

      !! ice conditions set from inputs in WIM2d_f2py or in mex gateway function
      !! > icec,iceh,dfloe
      ICE_MASK = 0.0
      where (icec.gt.min_conc)
         ICE_MASK = 1.0
      end where

      !! set wave grid
      !! SET_WAVES = .false. so don't set initial Tp,Hs,mwd,WAVE_MASK as well
      call waves_init(SET_WAVES)
      !print*,freq_vec

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! initialise WIM (mod_wavesice)
      !! set the attenuation coefficients;
      call wim_init(STEADY)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! get dir spec directly from Sdir_io in WIM2d_f2py

      call do_wim_integrals(0) ! calc Hs,Tp,mwd,WAVE_MASK from sdf_dir
         !! The "0" means it is relative to the water (displacement is
         !!  not reduced in the ice)
      !print*,minval(Hs),maxval(Hs)

      if (ndir.eq.1) then
         !!single direction:
         !!-set it to mean over the wave mask
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
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (DO_CHECK_INIT.eq.1) then
         print*,'Saving initial conditions to '//                       &
     &            trim(outdir)//'/binaries'
         print*,' '
         call check_grid(outdir)
         call check_init(outdir)
      end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!set number of time steps
      dt    = CFL*dx/amax
      x_ext = idm*dx/1.0e3!km
      y_ext = jdm*dy/1.0e3!km
      nt    = nint(duration/dt)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!print info
      call print_screen_info()

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      print*,'Beginning main integration...'
      call wim2d_run_common()

      et = etime(tictoc)!!finish etime
      print*,' '
      print*,'******************************************************'
      print*,'run of 2d WIM2d finished'
      !print *,'Elapsed time (s):',et
      print *,'Elapsed time (min):',et/60.0                             !&
!    &       ,', user:', tictoc(1),                                     &
!    &       ,', sys:' ,tictoc(2)
      print*,'******************************************************'
      print*,' '
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!print info to screen again
      call print_screen_info()

      !!print  some outputs to log file
      call print_log_info()

      !! put in same order as binary outputs:
      out_arrays(:,:,1) = dfloe     ! maximum floe size
      out_arrays(:,:,2) = tau_x     ! wave stress (x dirn)
      out_arrays(:,:,3) = tau_y     ! wave stress (y dirn)
      out_arrays(:,:,4) = Hs        ! significant wave height
      out_arrays(:,:,5) = Tp        ! mean wave period
      !out_arrays(:,:,6) = mwd       ! mean wave direction

      end subroutine wim2d_run_vSdir
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine wim2d_run_common()

      use mod_gridinfo
      use mod_wim_prams
      use mod_common_wim
      use mod_wim_init, only: freq_vec,n_wavdir,n_wave_freq
      use mod_wavesice, only: check_prog,wim_step

      implicit none

      integer  :: nMIZ

      !!TODO remove from wim2d_run_io between **
      real,dimension(idm,jdm) :: MIZ_MASK,tmp1

      logical :: critter
      integer  :: n,reps

      !! MAIN INTEGRATION
      !! - Common to WIM2d_run & WIM2d_run_io subroutines

      !nt = 13!stop at a certain point for testing
      !nt = 10*nt!extend to see true steady state
      reps  = 10!frequency of dumping (every 10 time steps dump binary file)

      do n=1,nt

         critter  = (DO_CHECK_PROG.eq.1).and.(mod(n-1,reps).eq.0)

         if (critter) then
            !!dump fields to binary
            call check_prog(SCATMOD,n-1,0,outdir)
         end if

         print*,'n/nt',n,nt
         call wim_step(SCATMOD,ADV_DIM,ADV_OPT                          &
     &                 ,DO_BREAKING,STEADY)

      end do!! n - finish time stepping

      if (DO_CHECK_FINAL.eq.1) then
         !!dump final output
         call check_prog(SCATMOD,n-1,1,outdir)
      end if

      MIZ_MASK = 0.0
      where ((dfloe.gt.0.0).and.(dfloe.lt.dfloe_pack_init))
         MIZ_MASK = 1.0
      end where

      !!Dmax range in MIZ
      Dmax_max = maxval(MIZ_MASK*dfloe)
      tmp1     = dfloe+1.0e4*(1.0-MIZ_MASK)
      Dmax_min = minval(tmp1)

      W_MIZ = 0.0
      if (.true.) then
         !! NB this definition won't always work
         !! for all configurations
         nMIZ  = sum(MIZ_MASK(:,1))
         W_MIZ = nMIZ*dx/1.0e3
         !!
         print*,' '
         print*,'MIZ width (km)',W_MIZ
      endif
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

      end subroutine wim2d_run_common
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine mex_io_gate(y1,y2,y3,y4,y5                             &
     &                  ,x1,x2,x3,x4,x5,x6                              &
     &                  ,int_prams,real_prams                           &
     &                  ,nx,ny,ni,nr)
      
      !! gateway between run_wim2d_io_mex.m,
      !! mex function WIM2d_run_io_mex.F,
      !! and wim2d_run_io subroutine in this file
      !! - 3d arrays hard to do in mex/fortran

      implicit none

      integer,intent(in)                  :: nx,ny,ni,nr
      integer,dimension(ni)               :: int_prams
      real,dimension(nr)                  :: real_prams
      real,dimension(nx*ny),intent(in)    :: x1,x2,x3,x4,x5,x6
      real,dimension(nx*ny),intent(out)   :: y1,y2,y3,y4,y5

      integer                 :: i,j,s
      real,dimension(nx,ny,6) :: in_arrays
      real,dimension(nx,ny,5) :: out_arrays

      do j=1,ny
      do i=1,nx
         s                 = i+(j-1)*nx
         in_arrays(i,j,1)  = x1(s)
         in_arrays(i,j,2)  = x2(s)
         in_arrays(i,j,3)  = x3(s)
         in_arrays(i,j,4)  = x4(s)
         in_arrays(i,j,5)  = x5(s)
         in_arrays(i,j,6)  = x6(s)
      end do
      end do

      call wim2d_run_io(in_arrays,out_arrays                            &
     &                  ,int_prams,real_prams,nx,ny)

      do j=1,ny
      do i=1,nx
         s     = i+(j-1)*nx
         y1(s) = out_arrays(i,j,1)
         y2(s) = out_arrays(i,j,2)
         y3(s) = out_arrays(i,j,3)
         y4(s) = out_arrays(i,j,4)
         y5(s) = out_arrays(i,j,5)
      end do
      end do

      end subroutine mex_io_gate
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine mex_io_gate_v2(y3d,x3d,int_prams,real_prams            &
     &                  ,nx,ny,ni,nr,nx3d,ny3d)
      
      !! gateway between run_wim2d_io_mex.m,
      !! mex function WIM2d_run_io_mex.F,
      !! and wim2d_run_io subroutine in this file
      !! - 3d arrays hard to do in mex/fortran

      implicit none

      integer                    ,intent(in)    :: nx,ny,ni,nr,nx3d,ny3d
      integer,dimension(ni)      ,intent(in)    :: int_prams
      real,dimension(nr)         ,intent(in)    :: real_prams
      real,dimension(nx*ny*nx3d) ,intent(in)    :: x3d !all 2d input arrays merged into 1 3d array
      real,dimension(nx*ny*ny3d) ,intent(out)   :: y3d !all 2d output arrays merged into 1 3d array

      integer                    :: i,j,s
      real,dimension(nx,ny,nx3d) :: in_arrays
      real,dimension(nx,ny,ny3d) :: out_arrays

      in_arrays   = reshape(x3d,(/nx,ny,nx3d/))

      print*,'before wim2d_run_io:',shape(in_arrays),nx,ny,nx3d
      call wim2d_run_io(in_arrays,out_arrays                            &
     &                  ,int_prams,real_prams,nx,ny)

      print*,'after wim2d_run_io:',shape(out_arrays),nx,ny,ny3d
      y3d   = reshape(out_arrays,(/nx*ny*ny3d/))

      end subroutine mex_io_gate_v2
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine mex_io_gate_vSdir(xy4d,y3d,x3d,int_prams,real_prams    &
     &                  ,Tp_io,mwd_io,ni,nr,nx3d,ny3d,nxy4d)
      
      !! gateway between run_wim2d_io_mex.m,
      !! mex function WIM2d_run_io_mex.F,
      !! and wim2d_run_io subroutine in this file
      !! - 3d arrays hard to do in mex/fortran

      use mod_common_wim, only: freq_vec,wavdir,sdf_dir                 &
     &                          ,icec,iceh,dfloe

      implicit none

      include "grid_info.h"
      include "wave_info.h"

      integer,parameter :: nx    = ii
      integer,parameter :: ny    = jj
      integer,parameter :: nw    = n_wave_freq
      integer,parameter :: ndir  = n_wavdir
      integer                    ,intent(in)    :: ni,nr,nx3d,ny3d,nxy4d
      real                       ,intent(in)    :: Tp_io,mwd_io
      integer,dimension(ni)      ,intent(in)    :: int_prams
      real,dimension(nr)         ,intent(in)    :: real_prams
      real,dimension(nx*ny*nx3d) ,intent(in)    :: x3d !all 2d input arrays merged into 1 3d array
      real,dimension(nx*ny*ny3d) ,intent(out)   :: y3d !all 2d output arrays merged into 1 3d array
      real,dimension(nx*ny*nxy4d),intent(inout) :: xy4d !4d input array (Sdir)

      integer  :: i,j,s
      real,dimension(nx,ny,nx3d) :: in_arrays
      real,dimension(nx,ny,ny3d) :: out_arrays

      !! initial ice conditions
      in_arrays   = reshape(x3d,(/nx,ny,nx3d/))
      icec        = in_arrays(:,:,1)
      iceh        = in_arrays(:,:,2)
      dfloe       = in_arrays(:,:,3)

      !! get input directional wave spectrum
      ! sdf_dir  = Sdir_io
      sdf_dir  = reshape(xy4d,(/nx,ny,ndir,nw/))

      print*,nw,(nw.eq.1),Tp_io
      if (nw.eq.1) then
         freq_vec(1) = 1./Tp_io
         print*,Tp_io,freq_vec
      end if
      if (ndir.eq.1) then
         wavdir(1)   = mwd_io
         print*,mwd_io
      end if

      call wim2d_run_vSdir(out_arrays,int_prams,real_prams)

      !! get updated directional wave spectrum
      ! Sdir_io  = sdf_dir
      xy4d  = reshape(sdf_dir,(/nx*ny*ndir*nw/))
      y3d   = reshape(out_arrays,(/nx*ny*ny3d/))

      end subroutine mex_io_gate_vSdir

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine print_screen_info()

      use mod_gridinfo
      use mod_common_wim
      use mod_wim_init, only: freq_vec,n_wavdir,n_wave_freq
      use mod_wavesice, only: amax,amin,CFL,dt

      implicit none

      integer,parameter :: ndir  = n_wavdir
      integer,parameter :: nw    = n_wave_freq

      print*,' '
      print*,'******************************************************'
      print*,'CFL number:                       ' ,CFL
      print*,'Maximum wave group velocity (m/s):' ,amax
      print*,'Time step (s):                    ' ,dt
      print*,'Number of time steps:             ' ,nt
      print*,'Time interval (h):                ' ,duration/60.0/60.0
      print*,' '
      print*,'Grid dimensions:                  ' ,idm,jdm
      print*,'Spatial resolution (km):          ' ,dx/1.0e3,dy/1.0e3
      print*,'Extent of domain (km):            ' ,x_ext,y_ext
      print*,' '
      print*,'Minimum period (s):               ' ,1.0/freq_vec(nw)
      print*,'Maximum period (s):               ' ,1.0/freq_vec(1)
      print*,'Number of wave frequencies:       ' ,nw
      print*,'Number of wave directions:        ' ,ndir
      print*,'Directional resolution (degrees): ' ,360.0/(1.0*ndir)
      print*,'******************************************************'
      print*,' '

      end subroutine print_screen_info
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine print_log_info()
      !!print text to log file during run

      use mod_gridinfo
      use mod_wim_prams
      use mod_common_wim
      use mod_wim_init, only: freq_vec,n_wavdir,n_wave_freq
      use mod_wavesice, only: amax,amin,CFL,dt

      implicit none

      integer,parameter :: ndir  = n_wavdir
      integer,parameter :: nw    = n_wave_freq
      integer,parameter :: fid=99

      print*,'Printing info to log file: ',trim(log_file)
      open(unit=fid,file=trim(log_file),status = 'replace')

      write(fid,'(a)'),'***********************************************'
      write(fid,'(a)'),'Outer subroutine:'
      write(fid,'(a,a)'),'>> ',this_subr
      write(fid,'(a)'),'***********************************************'

      write(fid,'(a)'),' '
      write(fid,'(a)'),'***********************************************'
      write(fid,'(a)'),'Main parameters:'
      write(fid,'(a,i2.2)'),'SCATMOD:                          ',SCATMOD
      write(fid,'(a,i2.2)'),'ADV_DIM:                          ',ADV_DIM
      write(fid,'(a,i2.2)'),'ADV_OPT:                          ',ADV_OPT
      write(fid,'(a,i2.2)'),'DO_BREAKING:                      '        &
     &   ,DO_BREAKING
      write(fid,'(a,i2.2)'),'STEADY:                           ',STEADY
      write(fid,'(a)'),'***********************************************'

      write(fid,'(a)'),' '
      write(fid,'(a)'),'***********************************************'
      write(fid,'(a)'),'WIM parameters:'
      write(fid,'(a,f4.2)'),'Brine volume fraction:      '  ,vbf
      write(fid,'(a,e10.3)'),'Youngs modulus (Pa):        ' ,young
      write(fid,'(a,e10.3)'),'Flexural strength (Pa):     ' ,sigma_c
      write(fid,'(a,e10.3)'),'Breaking strain:            ' ,epsc
      write(fid,'(a,f5.2)'),'Damping (Pa.s/m):           '  ,visc_rp
      write(fid,'(a)'),'***********************************************'

      write(fid,'(a)'),' '
      write(fid,'(a)'),'***********************************************'
      write(fid,'(a)'),'Other Parameters:'
      write(fid,'(a,f6.1)'),'Time step (s):                    ' ,dt
      write(fid,'(a,f4.3)'),'CFL number:                       ' ,CFL
      write(fid,'(a,f5.2)'),'Maximum wave group velocity (m/s):' ,amax
      write(fid,'(a,f6.1)'),'Time step (s):                    ' ,dt
      write(fid,'(a,i4.4)'),'Number of time steps:             ' ,nt
      write(fid,'(a,f5.2)'),'Time interval (h):                '        &
     &    ,duration/60.0/60.0
      write(fid,'(a)'),'***********************************************'

      write(fid,'(a)'),' '
      write(fid,'(a)'),'***********************************************'
      write(fid,'(a,i4.4,a,i4.4)'),'Grid dimensions:                  ' &
     &                        ,idm,' ',jdm
      write(fid,'(a,f4.1,a,f4.1)'),'Spatial resolution (km):          ' &
     &                        ,dx/1.0e3,' ',dy/1.0e3
      write(fid,'(a,f6.1,a,f6.1)'),'Extent of domain (km):            ' &
     &                        ,x_ext,' ',y_ext

      write(fid,'(a)'),' '
      write(fid,'(a,f5.2)'),'Minimum period (s):               '        &
     &                   ,1.0/freq_vec(nw)
      write(fid,'(a,f5.2)'),'Maximum period (s):               '        &
     &                   ,1.0/freq_vec(1)
      write(fid,'(a,i2.2)'),'Number of wave frequencies:       ' , nw
      write(fid,'(a,i3.3)'),'Number of wave directions:        ' , ndir
      write(fid,'(a,f5.2)'),'Directional resolution (degrees): '
     &                   ,360.0/(1.0*ndir)
      write(fid,'(a)'),'***********************************************'

      write(fid,'(a)'),' '
      write(fid,'(a)'),'***********************************************'
      write(fid,'(a)'),'Diagnostics:'
      write(fid,'(a,f6.1)'),'MIZ width (km): ',W_MIZ
      write(fid,'(a,f6.1,a,f6.1)'),'Dmax range in MIZ (m): '            &
     &                          ,Dmax_min,' ',Dmax_max
      write(fid,'(a,e10.3,a,e10.3)'),'tau_x range (Pa): '               &
     &                          ,taux_min,' ',taux_max
      write(fid,'(a,e10.3,a,e10.3)'),'tau_y range (Pa): '               &
     &                          ,tauy_min,' ',tauy_max
      write(fid,'(a)'),'***********************************************'

      write(fid,'(a,f6.1)'),' '
      write(fid,'(a)'),'***********************************************'
      write(fid,'(a,f7.1)'),'Elapsed time (min):',et/60.0
      write(fid,'(a)'),'***********************************************'
      write(fid,'(a)'),' '

      !!close file
      close(fid)

      end subroutine print_log_info
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine read_infile_dirs(indir,outdir)
      !!determines where to read inputs/outputs from
      
      implicit none
      
      character(len=80),intent(out) :: indir,outdir
      !character(len=80),intent(in)  :: infile
      integer,parameter :: fid=99

      open(unit=fid,file=infile_dirs,action='read')
      read(fid,'(a)'), indir
      read(fid,'(a)'), outdir
      close(fid)

      end subroutine read_infile_dirs
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine read_infile(wim_prams_in,duration)
      !!determines where to read inputs/outputs from
      
      implicit none
      
      real,intent(out)  :: wim_prams_in(2),duration
      integer,parameter :: fid=99
      character(len=80) :: ctmp

      wim_prams_in   = 0.0
      print*,' '
      print*,'******************************************************'
      print*,'Opening infile: ',infile

      open(unit=fid,file=infile,action='read')

      !! int_prams
      read(fid,'(a)'), ctmp
      read(fid,'(i2.2,a)'), SCATMOD     ,ctmp
      print*,'SCATMOD: ',SCATMOD
      read(fid,'(i2.2,a)'), ADV_DIM     ,ctmp
      print*,'ADV_DIM: ',ADV_DIM
      read(fid,'(i2.2,a)'), ADV_OPT     ,ctmp
      print*,'ADV_OPT: ',ADV_OPT
      read(fid,'(i2.2,a)'), DO_BREAKING ,ctmp
      print*,'DO_BREAKING: ',DO_BREAKING
      read(fid,'(i2.2,a)'), STEADY      ,ctmp
      print*,'STEADY: ',STEADY
      !!
      read(fid,'(i2.2,a)'), DO_CHECK_FINAL,ctmp
      print*,'DO_CHECK_FINAL: ',DO_CHECK_FINAL
      read(fid,'(i2.2,a)'), DO_CHECK_PROG ,ctmp
      print*,'DO_CHECK_PROG: ',DO_CHECK_PROG
      read(fid,'(i2.2,a)'), DO_CHECK_INIT ,ctmp
      print*,'DO_CHECK_INIT: ',DO_CHECK_INIT

      !! real_prams
      read(fid,'(a)'), ctmp
      read(fid,'(a)'), ctmp
      read(fid,'(e10.3,a)'),  wim_prams_in(1) ,ctmp
      print*,"Young's modulus: ",wim_prams_in(1)
      read(fid,'(f5.2,a)'),   wim_prams_in(2) ,ctmp
      print*,"RP damping: ",wim_prams_in(2)
      read(fid,'(f7.2,a)'),   duration        ,ctmp
      print*,"Duration (s): ",duration
      print*,'******************************************************'
      print*,' '
      !!
      close(fid)

      end subroutine read_infile
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module mod_wim2d_run
