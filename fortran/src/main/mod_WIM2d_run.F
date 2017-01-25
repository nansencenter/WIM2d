      !mod_WIM2d_run.F
      !main code
      !Author: Timothy Williams
      !Date:   20141118, 17:20:34 CET

      module mod_wim2d_run

      use mod_year_info
      use mod_common_wim

      implicit none

      private

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! variables common to subroutines here
!     integer  :: SCATMOD,ADV_DIM,ADV_OPT,STEADY
!     integer  :: DO_CHECK_PROG,DO_CHECK_INIT,DO_CHECK_FINAL
!     integer  :: BRK_OPT  !breaking criterion to use
!                          ! 0: no breaking; 1: Vernon; 2: Marchenko; 3: Mohr-Coulomb
!     integer  :: DO_ATTEN !can turn off attenuation to test advection
!     !integer  :: FSD_OPT  ! which FSD parameterisation to use (0,1)
!                          ! 0: RG (Dumont et al, 2011);1: smooth FSD
!     integer  :: FSD_OPT=1 !TODO make this an input
      integer  :: YOUNG_OPT=1
      integer,parameter :: n_param_vec   = 30

      real     :: x_ext,y_ext,u_ref,duration!,CFL
      integer  :: nt
      integer  :: dumpfreq = 10!frequency of dumping (every 10 time steps dump binary file)

      character(len=80) :: log_file,this_subr,indir
      character(len=*),parameter :: infile_dirs = 'infile_dirs.txt'

      real  :: W_MIZ,taux_max,tauy_max                                  &
     &   ,taux_min,tauy_min,Dmax_min,Dmax_max                           &

      real*4 :: et,tictoc(2) !!print time
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!subroutines are public (& n_param_vec)
      !! -also:
      !!   n_param_vec (for mex functions)
      !!   STEADY      (for f2py)
      public   :: wim2d_run
     &            ,wim2d_run_io                                         &
     &            ,wim2d_run_vSdir                                      &
     &            ,py_wim2d_run_vSdir_gate                              &
     &            ,mex_io_gate_v2                                       &
     &            ,mex_io_gate_vSdir                                    &
     &            ,mex_io_gate_vSdir_mesh                               &
     &            ,read_param_vec                                       &
     &            ,n_param_vec,STEADY
      
      contains
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine wim2d_run

      use mod_gridinfo     !! grid info
      use mod_wim_init     !! initial wave fields
      use mod_adv_atten    !! advection routine

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! import main code for WIM
      use mod_wavesice, only:
     &    amax,amin,dt                                                  &
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

      integer,parameter :: ndir        = n_wavdir
      integer,parameter :: nw          = n_wave_freq
      integer           :: w,wth

      !! function to print time
      real     :: etime
      real     :: wim_prams_in(2)
      real*8   :: init_cons(6)

      !! does infile* exist?
      logical  :: fex,fex2
      character(len=*),parameter :: infile = 'infile_nonIO.txt'
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      PI       = 4.0*atan(1.0)
      et       = etime(tictoc)!!startup etime
      !!
      wim_prams_in   = 0.0

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!locations of inputs/outputs
      inquire(file=infile_dirs,exist=fex)
      if (fex) then
         call read_infile_dirs(indir)
      else
         indir    = 'grid'
         outdir   = 'out'
      end if

      !!check if dir's exist & make them if necessary
      call create_directory(trim(outdir))
      call create_directory(trim(outdir)//'/diagnostics')
      call create_directory(trim(outdir)//'/diagnostics/global')
      call create_directory(trim(outdir)//'/diagnostics/local')
      call create_directory(trim(outdir)//'/binaries/')
      call create_directory(trim(outdir)//'/binaries/init')
      call create_directory(trim(outdir)//'/binaries/prog')
      call create_directory(trim(outdir)//'/binaries/final')
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      !log_file    = trim(outdir)//"/log/wim2d.log"
      this_subr   = 'wim2d_run'

      print*,'****************************************'
      print*,'Initialising...'
      print*,'****************************************'
      print*,' '

      inquire(file=infile,exist=fex2)
      if (fex2) then
         call read_infile(init_cons,infile)

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!initialisation of grid, ice con's
         call set_wim_prams() ! - ice properties like breaking strain
      else
         stop(infile//' not present')
!        SCATMOD     = 0
!        ADV_DIM     = 2
!        ADV_OPT     = 2
!           !!0, not periodic; 1, periodic in both x,y; 2, periodic in y only;
!        BRK_OPT     = 1
!        STEADY      = 1
!        DO_ATTEN    = 1

!        DO_CHECK_PROG  = 1
!        DO_CHECK_INIT  = 1
!        DO_CHECK_FINAL = 1

!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        !!initialisation of grid, ice con's
!        call set_wim_prams(YOUNG_OPT) ! - ice properties like Young's modulus

      end if

      if ((ADV_DIM.eq.2).and.(jj.eq.1)) then
         print*,'Invalid value of ADV_DIM (2) for jj = 1'
         print*,'Stopping'
         print*,''
         STOP
      end if

      call get_grid_file(indir) ! - grid: X,Y,LAND_MASK

      !! set idealised incident ice parameters
      call icecons_ideal(real(init_cons(4:6)))      ! - initial icec,iceh,dfloe

      !! set idealised incident wave parameters
      call wavecons_ideal(real(init_cons(1:3)))

      !! set wave grid
      call waves_init_freq()
      call waves_init_dir()

      !! create full incident wave spectrum
      !! - use preset wave parameters (Tp,Hs,mwd,WAVE_MASK)
      !!   to fill in the incident wave spectrum
      !!   with a Bretschneider spectrum
      call set_inc_waves()

      !! initialise WIM (mod_wavesice)
      !! set the attenuation coefficients;
      !call wim_init(STEADY,DO_ATTEN,CFL)
      call wim_init()
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (DO_CHECK_INIT.eq.1) then
         call check_grid()
         call check_init()
      end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! fill wave spectra arrays using sdf_inc
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
      nt = ceiling(duration/dt)
      dt = duration/nt
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!print info
      call print_screen_info()

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      print*,'Beginning main integration...'
      call wim2d_run_common(1) !wim2d_run

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
      !call print_log_info()

      end subroutine wim2d_run
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine wim2d_run_io(in_arrays,out_arrays                      &
     &                       ,real_prams,nx,ny)

      use mod_gridinfo     !! grid info
      use mod_wim_init     !! initial wave fields
      use mod_adv_atten    !! advection routine

!     !! main code for WIM
!     use mod_wavesice, only:
!    &    amax,amin,CFL,dt    !need these to determine number of time steps
!    &   ,wim_init,wim_step   !main subroutines
!    &   ,check_grid,check_init,check_prog   !subroutines for dumping fields

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! import main code for WIM
      use mod_wavesice, only:
     &    amax,amin,dt                                                  &
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
      integer,intent(in)   :: nx,ny
      real*8,intent(in)    :: real_prams(n_param_vec)
      real*8,dimension(nx,ny,6),intent(in)  :: in_arrays   !! c,h,Dmax,Hs,Tp,mwd
      real*8,dimension(nx,ny,6),intent(out) :: out_arrays  !! Dmax,Hs,Tp,taux,tauy,mwd

      real  :: PI,Tdiff,mwd_diff
      real  :: wim_prams_in(2)
      real*8   :: init_cons(6)
      integer,parameter :: ndir  = n_wavdir
      integer,parameter :: nw    = n_wave_freq
      integer           :: w,wth

      real,parameter    :: min_conc = 0.05

      !real  :: Tmean,Tvar,Tstdev

      !! function to print time
      real  :: etime

      !! does infile_dirs exist?
      logical :: fex
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call read_param_vec(init_cons,real_prams)
      !stop

      PI       = 4.0*atan(1.0)
      et       = etime(tictoc)!!startup etime

      if ((ADV_DIM.eq.2).and.(jj.eq.1)) then
         print*,'Invalid value of ADV_DIM (2) for jj = 1'
         print*,'Stopping'
         print*,''
         STOP
      end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!locations of inputs/outputs
      inquire(file=infile_dirs,exist=fex)
      if (fex) then
         call read_infile_dirs(indir)
      else
         indir    = 'grid'
         outdir   = 'out_io_1'
      end if

      !!check if dir's exist & make them if necessary
      call create_directory(trim(outdir))
      call create_directory(trim(outdir)//'/diagnostics')
      call create_directory(trim(outdir)//'/diagnostics/global')
      call create_directory(trim(outdir)//'/diagnostics/local')
      call create_directory(trim(outdir)//'/binaries/')
      call create_directory(trim(outdir)//'/binaries/init')
      call create_directory(trim(outdir)//'/binaries/prog')
      call create_directory(trim(outdir)//'/binaries/final')
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      this_subr   = 'wim2d_run_io'
      !print*,this_subr
      !print*,outdir

      print*,'****************************************'
      print*,'Initialising...'
      print*,'****************************************'
      print*,' '

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!initialisation of grid, ice con's
      call set_wim_prams()   ! - ice properties like breaking strain
      call get_grid_file(indir)     ! - grid: X,Y,LAND_MASK

      !!ice conditions from inputs
      icec     = real(in_arrays(:,:,1))
      iceh     = real(in_arrays(:,:,2))
      dfloe    = real(in_arrays(:,:,3))
      ICE_MASK = 0.0
      where (icec.gt.min_conc)
         ICE_MASK = 1.0
      end where

      !! SET_WAVES = .false. so don't set initial Tp,Hs,mwd,WAVE_MASK as well

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! set incident waves from inputs
      Hs          = real(in_arrays(:,:,4))
      Tp          = real(in_arrays(:,:,5))
      mwd         = real(in_arrays(:,:,6))
      WAVE_MASK   = 0.0
      where ((ICE_MASK.lt.1.0).and.(Hs.gt.0.0))
         !!kill waves inside ice and propagate them in ourselves
         WAVE_MASK   = 1.0
      end where
      !print*,WAVE_MASK
      ! print*,'# of wave cells',sum(WAVE_MASK)
      ! stop

      ! set waves grid
      call waves_init_freq()
      Tdiff = abs(1./freq_vec(1)-init_cons(2))
      if ((nw.eq.1).and.(Tdiff.gt.init_cons(2)*1.e-5) ) then
         print*,"Inconsistent values of parameter T_init and array Tp"
         print*,init_cons(2)
         print*,1./freq_vec(1)
         stop
      end if

      call waves_init_dir()
      mwd_diff = abs(wavdir(1)-init_cons(3))
      if ((ndir.eq.1).and.(mwd_diff.gt.abs(init_cons(3))*1.e-5) ) then
         print*,"Inconsistent values of parameter mwd_init"//
     &            " and array mwd"
         print*,init_cons(3)
         print*,wavdir(1)
         stop
      end if
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !! create full incident wave spectrum
      !! - use preset wave parameters (Tp,Hs,mwd,WAVE_MASK)
      !!   to fill in the incident wave spectrum
      !!   with a Bretschneider spectrum
      call set_inc_waves()

      !! initialise WIM (mod_wavesice)
      !! set the attenuation coefficients;
      !call wim_init(STEADY,DO_ATTEN,CFL)
      call wim_init()
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if ((DO_CHECK_INIT.eq.1).or.
     &    (DO_CHECK_PROG.eq.1).or.
     &    (DO_CHECK_FINAL.eq.1)) then
         !!save the grid
         call check_grid()
      end if

      if (DO_CHECK_INIT.eq.1) then
         call check_init()
      end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! fill wave spectra arrays using sdf_inc
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
      nt    = ceiling(duration/dt)
      dt    = duration/nt
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!print info
      call print_screen_info()

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      print*,'Beginning main integration...'
      call wim2d_run_common(1) !wim2d_run_io

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
      !call print_log_info()

      !! put in same order as binary outputs:
      out_arrays(:,:,1) = dble(dfloe)  ! maximum floe size
      out_arrays(:,:,2) = dble(tau_x)  ! wave stress (x dirn)
      out_arrays(:,:,3) = dble(tau_y)  ! wave stress (y dirn)
      out_arrays(:,:,4) = dble(Hs   )  ! significant wave height
      out_arrays(:,:,5) = dble(Tp   )  ! mean wave period
      out_arrays(:,:,6) = dble(mwd  )  ! mean wave period

      !!print  some outputs to log file
      !call print_log_info()

      end subroutine wim2d_run_io
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine wim2d_run_vSdir(out_arrays,nmesh_e,mesh_e)

      use mod_gridinfo     !! grid info
      use mod_wim_init     !! initial wave fields
      use mod_adv_atten    !! advection routine

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! import main code for WIM
      use mod_wavesice, only:
     &    amax,amin,dt                                                  &
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
      integer,parameter :: ndir  = n_wavdir
      integer,parameter :: nw    = n_wave_freq

      real*8,dimension(ii,jj,6),intent(out) :: out_arrays  !!Dmax,Hs,Tp,taux,tauy,mwd

      integer,intent(in)            :: nmesh_e
      real*8,intent(inout),optional :: mesh_e(nmesh_e,nmesh_vars)
      real :: mesh_e_dummy(nmesh_e,nmesh_vars)

      logical  :: INTERP_MESH
      real     :: PI

      logical,parameter :: SET_WAVES   = .false.
      integer           :: w,wth

      real,parameter    :: min_conc = 0.05

      real  :: Tmean,Tvar,Tstdev

      !! function to print time
      real  :: etime

      !! does infile_dirs exist?
      logical :: fex
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!using mesh?
      INTERP_MESH = .false.
      if (present(mesh_e)) then
         INTERP_MESH = .true.
      end if

      PI = 4.0*atan(1.0)
      et = etime(tictoc)!!startup etime

      if ((ADV_DIM.eq.2).and.(jj.eq.1)) then
         print*,'Invalid value of ADV_DIM (2) for jj = 1'
         print*,'Stopping'
         print*,''
         STOP
      end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!locations of inputs/outputs
      inquire(file=infile_dirs,exist=fex)
      if (fex) then
         call read_infile_dirs(indir)
      else
         indir    = 'grid'
         if (.not.INTERP_MESH) then
            outdir   = 'out_io_2'
         else
            outdir   = 'out_io_3'
         end if
      end if

      !!check if dir's exist & make them if necessary
      call create_directory(trim(outdir))
      call create_directory(trim(outdir)//'/diagnostics')
      call create_directory(trim(outdir)//'/diagnostics/global')
      call create_directory(trim(outdir)//'/diagnostics/local')
      call create_directory(trim(outdir)//'/binaries/')
      call create_directory(trim(outdir)//'/binaries/init')
      call create_directory(trim(outdir)//'/binaries/prog')
      call create_directory(trim(outdir)//'/binaries/final')
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      this_subr   = 'wim2d_run_io'

      print*,'****************************************'
      print*,'Initialising...'
      print*,'****************************************'
      print*,' '

      !print*,freq_vec
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!initialisation of grid, ice con's
      call set_wim_prams() ! - ice properties like breaking strain
      call get_grid_file(indir) ! - grid: X,Y,LAND_MASK

      !! ice conditions set from inputs in WIM2d_f2py or in mex gateway function
      !! > icec,iceh,dfloe
      ICE_MASK = 0.0
      where (icec.gt.min_conc)
         ICE_MASK = 1.0
      end where

      !! wave frequency grid (have already set if nw==1)
      if (nw.gt.1) then
         call waves_init_freq()
      end if

      !! wave directional grid (have already set if ndir==1)
      if (ndir.gt.1) then
         call waves_init_dir()
      end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! initialise WIM (mod_wavesice)
      !! set the attenuation coefficients;
      call wim_init()
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! have already got sdf_dir,sdf_inc directly from
      !! python/matlab input
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call do_wim_integrals(0) ! calc Hs,Tp,mwd,WAVE_MASK from sdf_dir
         !! The "0" means it is relative to the water (displacement is
         !!  not reduced in the ice)
      !print*,minval(Hs),maxval(Hs)

      if (DO_CHECK_INIT.eq.1) then
         print*,'Saving initial conditions to '//                       &
     &            trim(outdir)//'/binaries'
         print*,' '
         call check_grid()
         call check_init()
      end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!set number of time steps
      dt    = CFL*dx/amax
      x_ext = idm*dx/1.0e3!km
      y_ext = jdm*dy/1.0e3!km
      nt    = ceiling(duration/dt)
      dt    = duration/nt
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!print info
      call print_screen_info()

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      print*,'Beginning main integration...'
      if (INTERP_MESH) then
         ! run while doing interpolation onto the mesh points
         mesh_e_dummy   = real(mesh_e)
         call wim2d_run_common(nmesh_e,mesh_e_dummy)   !wim2d_run_vSdir
         mesh_e         = dble(mesh_e_dummy)
      else
         call wim2d_run_common(1)   !wim2d_run_vSdir
      end if

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

      !! put in same order as binary outputs:
      out_arrays(:,:,1) = dble(dfloe) ! maximum floe size
      out_arrays(:,:,2) = dble(tau_x) ! wave stress (x dirn)
      out_arrays(:,:,3) = dble(tau_y) ! wave stress (y dirn)
      out_arrays(:,:,4) = dble(Hs   ) ! significant wave height
      out_arrays(:,:,5) = dble(Tp   ) ! mean wave period
      out_arrays(:,:,6) = dble(mwd  ) ! mean wave direction

      !!print  some outputs to log file
      !call print_log_info()

      end subroutine wim2d_run_vSdir
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine wim2d_run_common(nmesh_e,mesh_e)

      use mod_gridinfo
      use mod_wim_init, only: freq_vec,n_wavdir,n_wave_freq
      use mod_wavesice, only: check_prog,wim_step,dt,DUMP_DIAG

      implicit none

      integer,intent(in)            :: nmesh_e
      real,intent(inout),optional   :: mesh_e(nmesh_e,nmesh_vars)
      integer  :: nMIZ,ndays_fwd
      !integer,parameter                   :: n_wim_step_prams  = 7
      !integer,dimension(n_wim_step_prams) :: wim_step_prams

      !!TODO remove from wim2d_run_io between **
      real,dimension(idm,jdm) :: MIZ_MASK,tmp1

      logical :: critter,WANT_DIAG
      integer :: n

      !! MAIN INTEGRATION
      !! - Common to WIM2d_run & WIM2d_run_io subroutines

      !nt = 13!stop at a certain point for testing
      !nt = 10*nt!extend to see true steady state
      WANT_DIAG   = (itest.gt.0).and.(jtest.gt.0)
      DUMP_DIAG   = WANT_DIAG

      call model_time_to_year_info(model_year_info
     &         ,model_day,model_seconds)

      log_file = trim(outdir)
     &   //'/diagnostics/global/WIM_diagnostics'
     &   //model_year_info%cdate//'T'
     &   //model_year_info%chour
     &   //model_year_info%cminute
     &   //model_year_info%csecond//'Z.txt'

      do n=1,nt


         critter  = (DO_CHECK_PROG.eq.1).and.(mod(n-1,dumpfreq).eq.0)

         if (critter) then
            !!dump fields to binary
            call check_prog(0)
         end if

         !!test if want to dump diagnostic text file at (itest,jtest)
         DUMP_DIAG   = (WANT_DIAG).and.(mod(n-1,dumpfreq).eq.0)

         print*,'n/nt',n,nt
!        wim_step_prams = (/SCATMOD,ADV_DIM,ADV_OPT                     &
!    &                     ,BRK_OPT,STEADY,DO_ATTEN,FSD_OPT/)
         if (.not.present(mesh_e)) then
            call wim_step(nmesh_e)
         else
            call wim_step(nmesh_e,mesh_e)
         end if

         if (.true.) then
            taux_max = maxval(tau_x)
            taux_min = minval(tau_x)
            tauy_max = maxval(tau_y)
            tauy_min = minval(tau_y)
            print*,' '
            print*,'time step no:',n,nt
            print*,'tau_x range (Pa)',taux_min,taux_max
            print*,'tau_y range (Pa)',tauy_min,tauy_max
            print*,' '
         end if

         !update global variables model_time,model_time_year_info
         model_seconds  = model_seconds+dt
         ndays_fwd      = floor(model_seconds/(24.*3600.))
         model_day      = model_day+ndays_fwd
         model_seconds  = model_seconds-ndays_fwd*(24.*3600.)
         call model_time_to_year_info(model_year_info
     &            ,model_day,model_seconds)

      end do!! n - finish time stepping

      if (DO_CHECK_FINAL.eq.1) then
         !!dump final output
         call check_prog(1)
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
         W_MIZ = nMIZ*dx
         !!
         print*,' '
         print*,'MIZ width (km)',W_MIZ/1.0e3
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

      call print_log_info()

      end subroutine wim2d_run_common
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine py_wim2d_run_vSdir_gate(Sdir_out,out_arrays            &
     &                 ,Sdir_in,ice_arrays,param_vec)

      use mod_wavesice, only: WAVE_MASK2

      implicit none

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!inputs
      real*8,dimension(ii*jj*n_wavdir*n_wave_freq)
     &   ,intent(in) :: Sdir_in                                         !! pass in directional spec
      real*8,dimension(ii,jj,3),intent(in)  :: ice_arrays              !! pass in  icec,iceh,dfloe
      real*8,intent(in)                     :: param_vec(n_param_vec)  !! parameters

      !! outputs
      real*8,dimension(ii,jj,6),intent(out) :: out_arrays  !! pass out Dmax,Hs,Tp,taux,tauy,mwd
      real*8,dimension(ii*jj*n_wavdir*n_wave_freq)
     &   ,intent(out) :: Sdir_out    !! pass in/out directional spec
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer,parameter :: nw    = n_wave_freq
      integer,parameter :: ndir  = n _wavdir
      !!
      integer  :: i,j,w,wth
      real     :: PI,adv_dir
      real     :: Tp_io,mwd_io,init_cons(6)

      PI = 4.*atan(1.)

      !!get parameters
      call read_param_vec(init_cons,param_vec)
      Tp_io    = init_cons(2)
      mwd_io   = init_cons(3)

      !! initial ice conditions
      icec     = real(ice_arrays(:,:,1))
      iceh     = real(ice_arrays(:,:,2))
      dfloe    = real(ice_arrays(:,:,3))

      !! get input directional wave spectrum
      sdf_dir  = reshape(real(Sdir_in),(/ii,jj,ndir,nw/))
      sdf_inc  = sdf_dir

      print*,nw,(nw.eq.1),Tp_io
      if (nw.eq.1) then
         freq_vec(1) = 1./Tp_io
         print*,Tp_io,freq_vec
      end if
      if (ndir.eq.1) then
         wavdir(1)   = mwd_io
         print*,mwd_io
      end if

      ! mesh_e is optional, but need a place holder for nmesh_e
      call wim2d_run_vSdir(out_arrays,1)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! get updated directional wave spectrum

      !! correction for STEADY (copied from mod_waves_ice.F)
      if (STEADY.eq.1) then
         do w=1,nw
         do wth=1,ndir
            adv_dir  = -PI/180.0*(wavdir(wth)+90.0)
            if (cos(adv_dir).ge.0.0) then
               where (WAVE_MASK2.gt.0.0)
                   sdf_dir(:,:,wth,w)   = sdf_inc(:,:,wth,w)
               end where
            end if
         end do
         end do
      end if

      !!TODO update Hs,Tp,... also?
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      Sdir_out = reshape(dble(sdf_dir),(/ii*jj*ndir*nw/))
      end subroutine py_wim2d_run_vSdir_gate
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     subroutine mex_io_gate_v2(y3d,x3d,int_prams,real_prams            &
!    &                  ,nx,ny,ni,nr,nx3d,ny3d)
      subroutine mex_io_gate_v2(y3d,x3d,real_prams
     &                  ,nx,ny,nx3d,ny3d)
      
      !! gateway between run_wim2d_io_mex.m,
      !! mex function WIM2d_run_io_mex.F,
      !! and wim2d_run_io subroutine in this file
      !! - 3d arrays hard to do in mex/fortran

      implicit none

      integer              ,intent(in) :: nx,ny,nx3d,ny3d!,ni,nr
      !integer,dimension(ni),intent(in) :: int_prams
      !!
      real*8,dimension(n_param_vec),intent(in)    :: real_prams
      real*8,dimension(nx*ny*nx3d) ,intent(in)    :: x3d !all 2d input arrays merged into 1 3d array
      real*8,dimension(nx*ny*ny3d) ,intent(out)   :: y3d !all 2d output arrays merged into 1 3d array

      integer                      :: i,j,s
      real*8,dimension(nx,ny,nx3d) :: in_arrays
      real*8,dimension(nx,ny,ny3d) :: out_arrays

      in_arrays   = reshape(x3d,(/nx,ny,nx3d/)) !convert to real*4 if necessary

      print*,'before wim2d_run_io:',shape(in_arrays),nx,ny,nx3d
      call wim2d_run_io(in_arrays,out_arrays                            &
     &                  ,real_prams,nx,ny)

      print*,'after wim2d_run_io:',shape(out_arrays),nx,ny,ny3d
      y3d   = reshape(out_arrays,(/nx*ny*ny3d/)) !convert to real*8 if necessary

      end subroutine mex_io_gate_v2
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine mex_io_gate_vSdir(xy4d,y3d,x3d,real_prams
     &                  ,nx3d,ny3d,nxy4d)
      
      !! gateway between run_wim2d_io_mex.m,
      !! mex function WIM2d_run_io_mex.F,
      !! and wim2d_run_io subroutine in this file
      !! - 3d arrays hard to do in mex/fortran

      use mod_wavesice, only: WAVE_MASK2

      implicit none

      include "grid_info.h"
      include "wave_info.h"

      integer,parameter :: nx    = ii
      integer,parameter :: ny    = jj
      integer,parameter :: nw    = n_wave_freq
      integer,parameter :: ndir  = n_wavdir
      !!
      integer              ,intent(in)    :: nx3d,ny3d,nxy4d!,ni,nr
      !!
      real*8,dimension(n_param_vec),intent(in)    :: real_prams
      real*8,dimension(nx*ny*nx3d) ,intent(in)    :: x3d !all 2d input arrays merged into 1 3d array
      real*8,dimension(nx*ny*ny3d) ,intent(out)   :: y3d !all 2d output arrays merged into 1 3d array
      real*8,dimension(nx*ny*nxy4d),intent(inout) :: xy4d !4d input array (Sdir)
      !!
      real     :: Tp_io,mwd_io
      real*8   :: init_cons(6)
      integer  :: i,j,s,w,wth
      integer  :: STEADY
      real     :: PI,adv_dir
      !!
      real*8,dimension(nx,ny,nx3d) :: in_arrays
      real*8,dimension(nx,ny,ny3d) :: out_arrays

      PI = 4.*atan(1.)

      call read_param_vec(init_cons,real_prams)
      Tp_io    = real(init_cons(2))
      mwd_io   = real(init_cons(3))

      !! initial ice conditions
      in_arrays   = reshape(real(x3d),(/nx,ny,nx3d/))!convert to real*4 if necessary
      icec        = in_arrays(:,:,1)
      iceh        = in_arrays(:,:,2)
      dfloe       = in_arrays(:,:,3)

      !! get input directional wave spectrum
      sdf_dir  = reshape(real(xy4d),(/nx,ny,ndir,nw/))!convert to real*4 if necessary
      sdf_inc  = sdf_dir

      if (nw.eq.1) then
         freq_vec(1) = 1./Tp_io
      end if
      if (ndir.eq.1) then
         wavdir(1)   = mwd_io
      end if

      !call wim2d_run_vSdir(out_arrays,int_prams,real(real_prams),1)
      call wim2d_run_vSdir(out_arrays,1)
         !1=dummy arg for nmesh_e

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! get updated directional wave spectrum

      !! correction for STEADY (copied from mod_waves_ice.F)
      ! STEADY   = int_prams(8)
      if (STEADY.eq.1) then
         do w=1,nw
         do wth=1,n_wavdir
            adv_dir  = -PI/180.0*(wavdir(wth)+90.0)
            if (cos(adv_dir).ge.0.0) then
               where (WAVE_MASK2.gt.0.0)
                   sdf_dir(:,:,wth,w)   = sdf_inc(:,:,wth,w)
               end where
            end if
         end do
         end do
      end if

      !!TODO update Hs,Tp,... also?
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!outputs
      xy4d  = reshape(dble(sdf_dir),(/nx*ny*ndir*nw/))!convert to real*8 if necessary
      y3d   = reshape(out_arrays,(/nx*ny*ny3d/))!convert to real*8 if necessary

      end subroutine mex_io_gate_vSdir
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine mex_io_gate_vSdir_mesh(xy4d,y3d,ymesh                  &
     &                  ,x3d,real_prams,nx3d,ny3d,nxy4d,nmesh_e)
      
      !! gateway between run_wim2d_io_mex.m,
      !! mex function WIM2d_run_io_mex.F,
      !! and wim2d_run_io subroutine in this file
      !! - 3d arrays hard to do in mex/fortran

      use mod_wavesice, only: WAVE_MASK2

      implicit none

      include "grid_info.h"
      include "wave_info.h"

      integer,parameter  :: nx    = ii
      integer,parameter  :: ny    = jj
      integer,parameter  :: nw    = n_wave_freq
      integer,parameter  :: ndir  = n_wavdir
      integer,parameter  :: nmesh_vars  = 6
      integer,intent(in) :: nx3d,ny3d,nxy4d,nmesh_e!,nr,ni
      !!
      real*8,dimension(n_param_vec),intent(in)  :: real_prams
      real*8,dimension(nx*ny*nx3d) ,intent(in)  :: x3d !all 2d input arrays merged into 1 3d array
      real*8,dimension(nx*ny*ny3d) ,intent(out) :: y3d !all 2d output arrays merged into 1 3d array
      !!
      real*8,dimension(nx*ny*nxy4d),intent(inout)        :: xy4d !4d input array (Sdir)
      real*8,dimension(nmesh_e*nmesh_vars),intent(inout) :: ymesh !
      !!
      integer  :: i,j,s,w,wth
      real     :: PI,adv_dir
      real     :: Tp_io,mwd_io
      real*8   :: init_cons(6)
      !!
      real*8,dimension(nx,ny,nx3d) :: in_arrays
      real*8,dimension(nx,ny,ny3d) :: out_arrays
      real*8,dimension(nmesh_e,nmesh_vars) :: mesh_e

      PI = 4.*atan(1.)

      call read_param_vec(init_cons,real_prams)
      Tp_io    = real(init_cons(2))
      mwd_io   = real(init_cons(3))

      !! initial ice conditions
      in_arrays   = reshape(real(x3d),(/nx,ny,nx3d/))!convert to real*4 if necessary
      icec        = in_arrays(:,:,1)
      iceh        = in_arrays(:,:,2)
      dfloe       = in_arrays(:,:,3)

      !! get input directional wave spectrum
      sdf_dir  = reshape(real(xy4d),(/nx,ny,ndir,nw/))!convert to real*4 if necessary
      sdf_inc  = sdf_dir

      if (nw.eq.1) then
         freq_vec(1) = 1./Tp_io
      end if
      if (ndir.eq.1) then
         wavdir(1)   = mwd_io
      end if

      mesh_e   = reshape(ymesh,(/nmesh_e,nmesh_vars/))
!     call wim2d_run_vSdir(out_arrays,int_prams,real(real_prams)
!    &         ,nmesh_e,mesh_e)!convert to real*4 if necessary
      
      call wim2d_run_vSdir(out_arrays,nmesh_e,mesh_e)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! get updated directional wave spectrum

      !! correction for STEADY (copied from mod_waves_ice.F)
      if (STEADY.eq.1) then
         do w=1,nw
         do wth=1,n_wavdir
            adv_dir  = -PI/180.0*(wavdir(wth)+90.0)
            if (cos(adv_dir).ge.0.0) then
               where (WAVE_MASK2.gt.0.0)
                   sdf_dir(:,:,wth,w)   = sdf_inc(:,:,wth,w)
               end where
            end if
         end do
         end do
      end if

      !!TODO update Hs,Tp,... also?
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!outputs
      xy4d  = reshape(dble(sdf_dir),(/nx*ny*ndir*nw/))
      y3d   = reshape(out_arrays   ,(/nx*ny*ny3d/))
      ymesh = reshape(mesh_e       ,(/nmesh_e*nmesh_vars/))

      end subroutine mex_io_gate_vSdir_mesh
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine print_screen_info()

      use mod_gridinfo
      use mod_wim_init, only: freq_vec,n_wavdir,n_wave_freq
      use mod_wavesice, only: amax,amin,dt

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
      use mod_wim_init, only: freq_vec,n_wavdir,n_wave_freq
      use mod_wavesice, only: amax,amin,dt

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
      write(fid,'(a,i2.2)'),'BRK_OPT:                          ',BRK_OPT
      if (BRK_OPT.eq.0) then
         write(fid,'(a)'),'(No breaking)'
      elseif (BRK_OPT.eq.1) then
         write(fid,'(a)'),'(Williams et al, 2013, Oc Mod)'
      elseif (BRK_OPT.eq.2) then
         write(fid,'(a)'),'(Marchenko)'
      elseif (BRK_OPT.eq.3) then
         write(fid,'(a)'),'(Mohr-Coulomb)'
      end if
      write(fid,'(a,i2.2)'),'STEADY:                           ',STEADY
      write(fid,'(a,i2.2)'),'DO_ATTEN:                         '        &
     &   ,DO_ATTEN
      write(fid,'(a)'),'***********************************************'

      write(fid,'(a)'),' '
      write(fid,'(a)'),'***********************************************'
      write(fid,'(a)'),'Other integer parameters:'
      write(fid,'(a,i2.2)'),'FSD_OPT:                          '        &
     &   ,FSD_OPT
      write(fid,'(a)'),'***********************************************'

      write(fid,'(a)'),' '
      write(fid,'(a)'),'***********************************************'
      write(fid,'(a)'),'WIM parameters:'
      write(fid,'(a,f4.2)'), 'Brine volume fraction:      '
     &                           ,wim_info%vbf
      write(fid,'(a,e10.3)'),'Youngs modulus (Pa):        '
     &                           ,wim_info%young
      write(fid,'(a,e10.3)'),'Flexural strength (Pa):     '
     &                           ,wim_info%sigma_c
      write(fid,'(a,e10.3)'),'Breaking stress (Pa):       '
     &                           ,wim_info%stress_c
      write(fid,'(a,e10.3)'),'Breaking strain:            '
     &                           ,wim_info%epsc
      write(fid,'(a,e10.3)'),'Cohesion (Pa):              '
     &                           ,wim_info%cohesion
      write(fid,'(a,e10.3)'),'Friction coefficient:       '
     &                           ,wim_info%friction
      write(fid,'(a,f5.2)'), 'RP drag (Pa.s/m):           '
     &                           ,wim_info%drag_rp
      write(fid,'(a,f5.2)'), 'WS viscosity (m^2/s):       '
     &                           ,wim_info%visc_ws
      write(fid,'(a)'),'***********************************************'

      write(fid,'(a)'),' '
      write(fid,'(a)'),'***********************************************'
      write(fid,'(a)'),'Other Parameters:'
      write(fid,'(a,f9.4)'),'Time step (s):                    ' ,dt
      write(fid,'(a,f4.3)'),'CFL number:                       ' ,CFL
      write(fid,'(a,f9.4)'),'Maximum wave group velocity (m/s):' ,amax
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
      write(fid,'(a,f6.3)'),'Minimum period (s):               '        &
     &                   ,1.0/freq_vec(nw)
      write(fid,'(a,f6.3)'),'Maximum period (s):               '        &
     &                   ,1.0/freq_vec(1)
      write(fid,'(a,i2.2)'),'Number of wave frequencies:       ' , nw
      write(fid,'(a,i3.3)'),'Number of wave directions:        ' , ndir
      write(fid,'(a,f6.3)'),'Directional resolution (degrees): '
     &                   ,360.0/(1.0*ndir)
      write(fid,'(a)'),'***********************************************'

      write(fid,'(a)'),' '
      write(fid,'(a)'),'***********************************************'
      write(fid,'(a)'),'Diagnostics:'
      write(fid,'(a,f6.1)'),'MIZ width (km): ',W_MIZ/1.e3
      write(fid,'(a,f9.4,a,f9.4)'),'Dmax range in MIZ (m): '            &
     &                          ,Dmax_min,' ',Dmax_max
      write(fid,'(a,e13.6,a,e13.6)'),'tau_x range (Pa): '               &
     &                          ,taux_min,' ',taux_max
      write(fid,'(a,e13.6,a,e13.6)'),'tau_y range (Pa): '               &
     &                          ,tauy_min,' ',tauy_max
      write(fid,'(a)'),'***********************************************'

!     write(fid,'(a,f6.1)'),' '
!     write(fid,'(a)'),'***********************************************'
!     write(fid,'(a,f7.1)'),'Elapsed time (min):',et/60.0
!     write(fid,'(a)'),'***********************************************'
!     write(fid,'(a)'),' '

      !!close file
      close(fid)

      end subroutine print_log_info
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine read_infile_dirs(indir)
      !!determines where to read inputs/outputs from
      
      implicit none
      
      character(len=80),intent(out) :: indir
      !character(len=80),intent(in)  :: infile
      integer,parameter :: fid=99

      open(unit=fid,file=infile_dirs,action='read')
      read(fid,'(a)'), indir
      read(fid,'(a)'), outdir
      close(fid)

      end subroutine read_infile_dirs
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine read_infile(init_cons,infile)
      !!determines where to read inputs/outputs from
      use mod_wavesice,    only: REF_HS_ICE,USE_ICE_VEL
      
      implicit none
      
      real*8,intent(out)  :: init_cons(6)
      character(len=*),intent(in)   :: infile
      integer,parameter :: fid=99
      integer           :: inf_version
      character(len=80) :: ctmp
      integer  :: m_year
      integer  :: m_month
      integer  :: m_day
      integer  :: m_hour
      integer  :: m_minute
      integer  :: m_second
      real*8   :: d_second

      print*,' '
      print*,'******************************************************'
      print*,'Opening infile: ',infile

      open(unit=fid,file=infile,action='read')

      ! ===================================================
      ! check infile version
      read(fid,'(i2.2,a)'), inf_version     ,ctmp

      ! set correct infile version number here
      ! - check if its correct
      if (inf_version.ne.4) then
         print*,'Wrong infile version'
         print*,'- get template and edit'
         stop
      end if
      ! ===================================================


      ! ===================================================
      !! int_prams
      read(fid,'(a)'), ctmp !blank
      read(fid,'(a)'), ctmp !comment
      read(fid,'(i2.2,a)'), SCATMOD     ,ctmp
      print*,'SCATMOD: ',SCATMOD
      read(fid,'(i2.2,a)'), ADV_DIM     ,ctmp
      print*,'ADV_DIM: ',ADV_DIM
      read(fid,'(i2.2,a)'), ADV_OPT     ,ctmp
      print*,'ADV_OPT: ',ADV_OPT
      !!
      read(fid,'(i2.2,a)'), BRK_OPT ,ctmp
      print*,'BRK_OPT: ',BRK_OPT
      read(fid,'(i2.2,a)'), STEADY      ,ctmp
      print*,'STEADY: ',STEADY
      read(fid,'(i2.2,a)'), DO_ATTEN ,ctmp
      print*,'DO_ATTEN: ',DO_ATTEN
      !!
      read(fid,'(i2.2,a)'), DO_CHECK_FINAL,ctmp
      print*,'DO_CHECK_FINAL: ',DO_CHECK_FINAL
      read(fid,'(i2.2,a)'), DO_CHECK_PROG ,ctmp
      print*,'DO_CHECK_PROG: ',DO_CHECK_PROG
      read(fid,'(i2.2,a)'), DO_CHECK_INIT ,ctmp
      print*,'DO_CHECK_INIT: ',DO_CHECK_INIT
      ! ===================================================

      ! ===================================================
      !! real_prams
      read(fid,'(a)'), ctmp !blank
      read(fid,'(a)'), ctmp !comment
      read(fid,'(e10.3,a)')      ,wim_info%young  ,ctmp
      print*,"Young's modulus: " ,wim_info%young
      read(fid,'(f5.2,a)')       ,wim_info%drag_rp,ctmp
      print*,"RP drag: "         ,wim_info%drag_rp
      read(fid,'(f5.2,a)')       ,wim_info%visc_ws,ctmp
      print*,"WS viscosity: "    ,wim_info%visc_ws
      read(fid,'(e10.3,a)')      ,wim_info%cohesion,ctmp
      print*,"cohesion: "        ,wim_info%cohesion
      read(fid,'(f5.2,a)')       ,wim_info%friction,ctmp
      print*,"friction: "        ,wim_info%friction
      read(fid,'(f7.2,a)')       ,CFL             ,ctmp
      print*,"CFL: "             ,CFL
      print*,'******************************************************'
      print*,' '
      ! ===================================================

      ! ===================================================
      !! other integer parameters
      print*,' '
      print*,'******************************************************'
      read(fid,'(a)'), ctmp !blank
      read(fid,'(a)'), ctmp !comment
      read(fid,'(i2.2,a)'), FSD_OPT ,ctmp
      print*,'FSD_OPT: ',FSD_OPT
      read(fid,'(i2.2,a)'), REF_HS_ICE ,ctmp
      print*,'REF_HS_ICE: ',REF_HS_ICE
      read(fid,'(i2.2,a)'), USE_ICE_VEL ,ctmp
      print*,'USE_ICE_VEL: ',USE_ICE_VEL
      print*,'******************************************************'
      ! ===================================================

      ! ===================================================
      ! initial conditions
      print*,' '
      print*,'******************************************************'
      read(fid,'(a)'), ctmp !blank
      read(fid,'(a)'), ctmp !comment
      read(fid,'(f5.2,a)'), init_cons(1),ctmp !Hs_init   
      print*,"Hs_init (m): ",init_cons(1)
      read(fid,'(f5.2,a)'), init_cons(2),ctmp !T_init    
      print*,"T_init (s): ",init_cons(2)
      read(fid,'(f5.2,a)'), init_cons(3),ctmp !dir_init  
      print*,"mwd_init (deg): ",init_cons(3)
      read(fid,'(f5.2,a)'), init_cons(4),ctmp !conc_init 
      print*,"conc_init: ",init_cons(4)
      read(fid,'(f5.2,a)'), init_cons(5),ctmp !h_init    
      print*,"h_init: ",init_cons(5)
      read(fid,'(f5.2,a)'), init_cons(6),ctmp !Dmax_init 
      print*,"Dmax_init: ",init_cons(6)
      print*,'******************************************************'
      ! ===================================================

      ! ===================================================
      ! time info
      read(fid,'(a)'), ctmp !blank
      read(fid,'(f7.2,a)')       ,duration        ,ctmp
      print*,"Duration (s): "    ,duration
      !!
      print*,' '
      print*,'******************************************************'
      read(fid,'(a)'), ctmp !blank
      read(fid,'(a)'), ctmp !comment
      read(fid,'(a)'), ctmp !start-date YYYYMMDDTHHMMSSZ + comments
      
      !!parse date
      read(ctmp(1:4)  ,'(i4.4)'), m_year  
      read(ctmp(5:6)  ,'(i2.2)'), m_month 
      read(ctmp(7:8)  ,'(i2.2)'), m_day   
      read(ctmp(10:11),'(i2.2)'), m_hour  
      read(ctmp(12:13),'(i2.2)'), m_minute
      read(ctmp(14:15),'(i2.2)'), m_second

      !!calc model start time etc
      d_second = m_second+60.*m_minute+3600.*m_hour
      call make_year_info(model_year_info,                              &
     &      m_year,m_month,m_day,d_second)
      model_day               = model_year_info%model_day
      model_seconds           = model_year_info%dsecond
      model_year_info_start   = model_year_info

      print*,'Start date & time: '//                                    &
     &   model_year_info%cyear //'/'//                                  &
     &   model_year_info%cmonth//'/'//                                  &
     &   model_year_info%cday//' '//
     &   model_year_info%chour  //':'//                                 &
     &   model_year_info%cminute//':'//                                 &
     &   model_year_info%csecond
      print*,'Starting model day: '                                     &
     &   ,model_year_info%model_day
      print*,'Seconds in starting model day: '                          &
     &   ,model_year_info%dsecond
      print*,'******************************************************'
      ! ===================================================


      !! diagnostic grid cell
      print*,' '
      print*,'******************************************************'
      read(fid,'(a)'), ctmp !blank
      read(fid,'(a)'), ctmp !comment
      read(fid,'(i6.6,a)'),itest,ctmp
      read(fid,'(i6.6,a)'),jtest,ctmp
      print*,'Diagnostics requested at (i,j):',itest,jtest
      read(fid,'(i6.6,a)'),dumpfreq,ctmp
      print*,'Diagnostic output freq:',dumpfreq
      print*,'******************************************************'

      close(fid)

      print*,' '
      !stop

      end subroutine read_infile
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine read_param_vec(init_cons,param_vec)
      !!determines where to read inputs/outputs from
      use mod_wavesice,    only: REF_HS_ICE,USE_ICE_VEL
      
      implicit none
      
      real*8,intent(out)   :: init_cons(6)
      real*8,intent(in)    :: param_vec(n_param_vec)
      integer  :: i,j,fid
      logical,parameter :: test_in=.false.

      if (test_in) then
         open(unit=fid,file='param_vec.txt',action='write')
         do i=1,n_param_vec
            write(fid,*) param_vec(i)
         end do
         close(unit=fid)
      end if
      i  = 0


      ! ===================================================
      !! int_prams
      i        = i+1
      SCATMOD  = int(param_vec(i))
      print*,'SCATMOD: ',SCATMOD
      i        = i+1
      ADV_DIM  = int(param_vec(i))
      print*,'ADV_DIM: ',ADV_DIM
      i        = i+1
      ADV_OPT  = int(param_vec(i))
      print*,'ADV_OPT: ',ADV_OPT
      !!
      i        = i+1
      BRK_OPT  = int(param_vec(i))
      print*,'BRK_OPT: ',BRK_OPT
      i        = i+1
      STEADY   = int(param_vec(i))
      print*,'STEADY: ',STEADY
      i        = i+1
      DO_ATTEN = int(param_vec(i))
      print*,'DO_ATTEN: ',DO_ATTEN
      !!
      i              = i+1
      DO_CHECK_FINAL = int(param_vec(i))
      print*,'DO_CHECK_FINAL: ',DO_CHECK_FINAL
      i              = i+1
      DO_CHECK_PROG  = int(param_vec(i))
      print*,'DO_CHECK_PROG: ',DO_CHECK_PROG
      i              = i+1
      DO_CHECK_INIT  = int(param_vec(i))
      print*,'DO_CHECK_INIT: ',DO_CHECK_INIT
      ! ===================================================

      ! ===================================================
      !! real_prams
      i              = i+1
      wim_info%young = real(param_vec(i))
      print*,"Young's modulus: ",wim_info%young
      !!
      i                 = i+1
      wim_info%drag_rp  = real(param_vec(i))
      print*,"RP drag: ",wim_info%drag_rp
      !!
      i                 = i+1
      wim_info%visc_ws  = real(param_vec(i))
      print*,"WS viscosity: ",wim_info%visc_ws
      !!
      i                 = i+1
      wim_info%cohesion = real(param_vec(i))
      print*,"cohesion: ",wim_info%cohesion
      !!
      i                 = i+1
      wim_info%friction = real(param_vec(i))
      print*,"friction: ",wim_info%friction
      !!
      i     = i+1
      CFL   = real(param_vec(i))
      print*,"CFL: ",CFL
      print*,'******************************************************'
      print*,' '
      ! ===================================================

      ! ===================================================
      !! other integer parameters
      print*,' '
      print*,'******************************************************'
      i        = i+1
      FSD_OPT  = int(param_vec(i))
      print*,'FSD_OPT: ',FSD_OPT
      i           = i+1
      REF_HS_ICE  = int(param_vec(i))
      print*,'REF_HS_ICE: ',REF_HS_ICE
      i           = i+1
      USE_ICE_VEL = int(param_vec(i))
      print*,'USE_ICE_VEL: ',USE_ICE_VEL
      print*,'******************************************************'
      ! ===================================================

      ! ===================================================
      ! initial conditions
      print*,' '
      print*,'******************************************************'
      do j=1,6
         i     = i+1
         init_cons(j)   = param_vec(i)
      end do
      print*,"Hs_init (m): "     ,init_cons(1)
      print*,"T_init (s): "      ,init_cons(2)
      print*,"mwd_init (deg): "  ,init_cons(3)
      print*,"conc_init: "       ,init_cons(4)
      print*,"h_init: "          ,init_cons(5)
      print*,"Dmax_init: "       ,init_cons(6)
      print*,'******************************************************'
      ! ===================================================

      ! ===================================================
      ! time info
      i        = i+1
      duration = real(param_vec(i))
      print*,"Duration (s): ",duration
      !!
      i              = i+1
      model_day      = int(param_vec(i))
      !!
      i              = i+1
      model_seconds  = dble(param_vec(i))
      call model_time_to_year_info(model_year_info,
     &      model_day,model_seconds)
      model_year_info_start   = model_year_info

      print*,' '
      print*,'******************************************************'
      print*,'Start date & time: '//                                    &
     &   model_year_info%cyear //'/'//                                  &
     &   model_year_info%cmonth//'/'//                                  &
     &   model_year_info%cday//' '//
     &   model_year_info%chour  //':'//                                 &
     &   model_year_info%cminute//':'//                                 &
     &   model_year_info%csecond
      print*,'Starting model day: '                                     &
     &   ,model_year_info%model_day
      print*,'Seconds in starting model day: '                          &
     &   ,model_year_info%dsecond
      print*,'******************************************************'
      ! ===================================================

      !! diagnostic grid cell
      print*,' '
      print*,'******************************************************'
      i     = i+1
      itest = int(param_vec(i))
      i     = i+1
      jtest = int(param_vec(i))
      print*,'Diagnostics requested at (i,j):',itest,jtest
      i        = i+1
      dumpfreq = int(param_vec(i))
      print*,'Diagnostic output freq:',dumpfreq
      print*,'******************************************************'
      print*,' '
      !stop


      end subroutine read_param_vec
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine create_directory( newDirPath )
      ! Author:  Jess Vriesema
      ! Date:    Spring 2011
      ! Purpose: Creates a directory at ./newDirPath
      
      implicit none
      
      character(len=*), intent(in) :: newDirPath
      character(len=256)           :: mkdirCmd
      logical                      :: dirExists
      
      ! Check if the directory exists first
      inquire( file=trim(newDirPath)//'/.', exist=dirExists )  ! Works with gfortran, but not ifort
      !inquire( directory=newDirPath, exist=dirExists )         ! Works with ifort, but not gfortran
      
      
      if (dirExists) then
         write (*,*) "Directory already exists: '"//                    &
     &      trim(newDirPath)//"'"
      else
         mkdirCmd = 'mkdir -p '//trim(newDirPath)
         write(*,'(a)') "Creating new directory: '"//                   &
     &      trim(mkdirCmd)//"'"
         call system( mkdirCmd )
      endif
      end subroutine create_directory
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module mod_wim2d_run
