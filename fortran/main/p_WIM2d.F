      !! test_waveadv_weno.F
      !! Author: Timothy Williams
      !! Date:   20140821, 12:22:17 CEST
      program main

      use mod_gridinfo     !!grid info
      use mod_wim_prams    !!ice parameters
      use mod_wim_init   !!initial wave fields
      use mod_common_wim   !!initial ice fields
      !use mod_file_utils   !!file reading/writing
      !use mod_waveadv_weno !!advection routine

      implicit none

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!variables
      integer,parameter :: SOLVER      = 1
      integer,parameter :: STEADY      = 1
      integer,parameter :: DO_ATTEN    = 1
      integer,parameter :: DO_BREAKING = 1
      integer,parameter :: USE_ICE_VEL = 0
      !!
      integer,parameter :: OPT      = 1
      integer,parameter :: bc_opt   = 2
         !1: periodic in x,y; 2: periodic in y only; 0: not periodic
      integer,parameter :: GRID_OPT = 1!option for land/ice/water config
      real,parameter    :: CFL      = 0.4
      !!
      integer,parameter :: nw    = n_wave_freq
      integer,parameter :: ndir  = n_wavdir
      !!
      real,dimension(idm,jdm)    :: mom0,mom2,var_strain,tau_x,tau_y
      real,dimension(idm,jdm,nw) ::                                     &
     &   ag_eff,ap_eff,wlng_ice,atten_nond,damping,disp_ratio

      real     :: PI
      real     :: adv_dir,dt,dtheta,om,dom,guess
      real     :: kwtr,kice,modT,int_adm,argT,argR
      real,dimension(nw)   :: wt_om,wt_simp,T,wlng,ag,ap
      integer  :: nt,n,i,j,w,wth
      logical  :: critter

      character(len=*),parameter :: outfile = "init"
      character(len=*),parameter :: outdir = "out/"
      character(len=3)           :: cts
      !character(len=80)          :: afile,bfile,stat1,stat2

      real :: et,etime,tictoc(2) !!print time
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      PI = 4.0*atan(1.0)
      et = etime(tictoc)!!startup etime

      print*,'****************************************'
      print*,'Initialising...'
      print*,'****************************************'
      print*,' '

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!initialisation
      call set_wim_prams()
      call get_grid(GRID_OPT)
      call ice_init(GRID_OPT)

      !!incident waves
      call waves_init(GRID_OPT)
      call set_inc_waves()

      do w=1,nw
      do wth=1,ndir
         where (WAVE_MASK.eq.1.0)
            sdf_dir(:,:,wth,w)   = sdf_inc(:,:,wth,w)
         end where
      end do
      end do

      if (.true.) then
         !!check initialisation
         !!(dump binaries to look at with test_WIM2d_F.m)
         call check_grid(GRID_OPT)
         call check_init(GRID_OPT)
      end if
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! define weights for numerical quadrature;
      if (nw.gt.1) then!! weights for integral over frequency
         !! (Simpson's rule);
         wt_simp     = 2.0
         wt_simp(1)  = 1.0
         wt_simp(nw) = 1.0

         w  = 0
         do while ((w+2).lt.nw)
            w  = w+2
            wt_simp(w)  = 4.0
         end do
      
         !!NB om needs to be equally spaced;
         dom   = 2*PI*(freq_vec(nw)-freq_vec(1))/(nw-1.0)
         wt_om = dom/3.0*wt_simp
      else
         wt_om = 1
      end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!get wave speeds and attenuation coefficients
      print*,'****************************************'
      print*,'Get attenuation coefficients...'
      print*,'****************************************'
      print*,' '

      !! Water wavelength and wave speed
      !! is a function only of wave period
      wlng     = gravity*T**2/(2*PI)
      ap       = sqrt(gravity*wlng/(2*PI))   ! Phase speed
      ag       = ap/2.0                      ! Group speed

      ag_eff      = 0.0
      ap_eff      = 0.0
      wlng_ice    = 0.0
      disp_ratio  = 1.0
      atten_nond  = 0.0
      damping     = 0.0

      do w=1,nw
         
         om    = 2*PI*freq_vec(1)
         guess = om**2/gravity!guess for root

         do j=1,jj
         do i=1,ii
         if (ICE_MASK(i,j).eq.1.0) then

            om = 2*PI*freq_vec(w)
            call RTparam_outer(damping(i,j,w),kice,kwtr,int_adm,        &
     &                  atten_nond(i,j,w),modT,argR,argT,               &
     &                  iceh(i,j),om,guess)

            guess             = kice
            disp_ratio(i,j,w) = kice/kwtr*modT
            wlng_ice(i,j,w)   = 2*PI/kice

            if (USE_ICE_VEL.eq.0) then
               !!use water group and phase velocity
               ag_eff(i,j,w)  = ag(w)
               ap_eff(i,j,w)  = ap(w)
            end if

         else

            !!no ice
            wlng_ice(i,j,w)   = wlng(w)
            ag_eff(i,j,w)     = ag(w)
            ap_eff(i,j,w)     = ap(w)

         end if
         end do
         end do
      end do

      if (DO_ATTEN.eq.0) then
         atten_nond  = 0.0
         damping     = 0.0
      end if

      nt = 100
      print*,' '
      print*,'******************************************************'
      print*,'Number of time steps:       ',nt
      print*,'Number of wave frequencies: ',nw
      print*,'Number of wave directions:  ',ndir
      print*,'******************************************************'
      print*,' '
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      

      print*,'Beginning main integration...'

      do n=2,nt
         print*,'n/nt',n,nt

         !! spectral moments;
         mom0  = 0.0
         mom2  = 0.0

         !! wave stresses;
         tau_x = 0.0
         tau_y = 0.0

         !! variance of strain;
         var_strain     = 0.0
         do w=1,nw

            !! CALC DIMENSIONAL ATTEN COEFF;
            atten_dim   = 0.0
            damp_dim    = 0.0
            do j = 1,jj
            do i = 1,ii

               !!top-up waves in wave mask if STEADY==1
               !!(steady-state solution);
               if ((WAVE_MASK(i,j).ge.0).and.(STEADY.eq.1))
                  !sdf_dir(i,j,J_STEADY,:) = sdf_inc(i,j,J_STEADY,:)
                  sdf_dir(i,j,:,:)  = sdf_inc(i,j,:,:)
               end if
               
               if ((ICE_MASK(i,j).eq.1.0).and.(DO_ATTEN.eq.1)) then
                  call floe_scaling(Dave,Dmax(i,j))

                  !! get expected no of floes met per unit
                  !!  distance if travelling in a line;
                  if (Dmax(i,j).lt.200.0) then
                     c1d = cice(i,j)/Dave
                     !! floes per unit length;
                  else
                     c1d = cice(i)/Dmax(i,j)
                     !! uniform lengths
                  end if

                  !! ENERGY attenuation coeff;
                  atten_dim(i,j) = atten_nond(i,j,jw)*c1d      !!scattering
                  damp_dim(i,j)  = 2*damping(i,j,jw)*cice(i,j) !!damping
               end if
            end do!! j
            end do!! i

            if (ndir==1) then
               if (SOLVER.ne.0) then
                  print*,'warning: changing SOLVER option'
     &                  ,'as not enough directions'
                  print*,'(ndir = ',ndir,')'
               end if
               SOLVER   = 0
            end if

            if (SOLVER.eq.0) then
               !! Simple attenuation scheme - doesn't conserve scattered energy
!              call adv_atten_simple(S_freq,sdf_dir(:,:,:,w)            &
!    &            ,tau_x_om,tau_y_om,                                   &
!    &            ,atten_dim,damp_dim,ag_eff(:,:,w),wavdir,dt)
            elseif (SOLVER.eq.1) then
               !! same as SOLVER==0, but scattered energy
               !! is distributed isotropically
!              call adv_atten_isotropic(S_freq,sdf_dir(:,:,:,w)         &
!    &            ,tau_x_om,tau_y_om,                                   &
!    &            ,atten_dim,damp_dim,ag_eff(:,:,w),wavdir,dt)
            end if

            !! integrate stress densities over frequency
            !! TODO: check if this is correct for ice-covered water
!           tmp1  = rhowtr*gravity*tau_x_om/ap_eff(:,:,w)   !![Pa*s]
!           tau_x = tau_x+wt_om(w)*tmp1                     !![Pa]
!           tmp1  = rhowtr*gravity*tau_y_om./ap_eff(:,:,w)  !![Pa*s]
!           tau_y = tau_y+wt_om(w)*tmp1                     !![Pa]

            !! INTEGRALS FOR BREAKING PROB:
            do j=1,jj
            do i=1,ii

               !! convert from water amp's to ice amp's;
               F     = disp_ratio(i,j,w)
               k_ice = 2*PI/wlng_ice(i,j,w)

               !! SPECTRAL MOMENTS;
               !!take abs as small errors can make S_freq negative
               om          = 2*PI*freq_vec(w)
               tmp         = wt_om(w)*S_freq(i,j)*F**2
               mom0(i,j)   = mom0(i,j)+abs( tmp )!!variance of displacement
               tmp         = wt_om(w)*om**2*S_freq(i,j)*F**2
               mom2(i,j)   = mom2(i,j)+abs( tmp )!!variance of speed

               if (ICE_MASK(i,j).eq.1.0) then
                  !! VARIANCE OF STRAIN;
                  tmp   = F*k_ice**2*hice(i,j)/2.0    !!strain conversion factor
                  tmp   = wt_om(w)*S_freq(i,j)*tmp**2 !!strain density
                  var_strain(i,j)   = var_strain(i,j)+abs( tmp )
               end if
            end do!! i, end spatial loop x;
            end do!! j, end spatial loop y;

         end do!! end spectral loop;

      end do

      et = etime(tictoc)!!finish etime
      print*,' '
      print*,'******************************************************'
      print*,'p_WIM2d.F finished'
      print *,'Elapsed time (s):',et                                    &
!    &       ,', user:', tictoc(1),                                     &
!    &       ,', sys:' ,tictoc(2)
      print*,'******************************************************'
      print*,' '

!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!      do j=1,jj
!!      do i=1,ii
!!         R(i,j)      = sqrt(X(i,j)**2+Y(i,j)**2)
!!         Theta(i,j)  = atan2(Y(i,j),X(i,j))
!!      end do
!!      end do
!!
!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!      !!Test situations
!!      if (OPT.eq.1) then
!!         !!constant speed
!!         print*,'******************************************'
!!         print*,'using constant speed: OPT = ',OPT
!!         print*,'******************************************'
!!         print*,''
!!         uc       = 30!const speed m/s
!!         xc       = 2*xm/3.0
!!         !adv_dir = 180.0!deg straight across
!!         adv_dir  = 135.0!deg
!!         u        = uc*cos(pi/180.0*adv_dir)
!!         v        = uc*sin(pi/180.0*adv_dir)
!!         dt       = CFL*dx/uc
!!         if (bc_opt.eq.1) then
!!            nt = 2*aint(2*xm/(uc*dt))  !test periodic boundary condition
!!         elseif (bc_opt.eq.2) then
!!            nt = aint(3*xm/(uc*dt))  !test periodic boundary condition
!!         else
!!            nt = aint(2*xm/(uc*dt))
!!         end if
!!         !!
!!         h  = 0.0
!!         do j=1,jj
!!         do i=1,ii
!!            if (X(i,j).gt.xc) then
!!               h(i,j)   = 1.0
!!            end if
!!         end do
!!         end do
!!
!!      elseif (OPT.eq.2) then
!!         !!variable speed
!!         uc = 30!max speed m/s
!!         xc = 2*xm/3.0
!!         u  = -uc*X/xm!speed variable with x
!!         v  = 0.0
!!         dt = CFL*dx/uc
!!         nt = 2*xm/(uc*dt)
!!         !!
!!         h  = 0.0
!!         do j=1,jj
!!         do i=1,ii
!!            if (X(i,j).gt.xc) then
!!               h(i,j)   = 1.0
!!            end if
!!         end do
!!         end do
!!
!!      elseif (OPT.eq.3) then
!!         !!pacman
!!         Rc       = xm/3.0
!!         Yc       = ym/12.0
!!         angrot   = (1/20.0)*pi/180.0!!radian/s
!!         u        = 0.0
!!         v        = 0.0
!!         !!
!!         h  = 0.0
!!         do j=1,jj
!!         do i=1,ii
!!            critter  = ((R(i,j).lt.Rc).and.(X(i,j).le.0))
!!            if (critter) then
!!               h(i,j)   = 1.0
!!            end if
!!            critter  = ((R(i,j).lt.Rc).and.(X(i,j).ge.0)                &
!!     &                  .and.(abs(Y(i,j)).ge.Yc))
!!            if (critter) then
!!               h(i,j)   = 1.0
!!            end if
!!            critter  = (R(i,j).lt.(1.45*Rc))
!!            if (critter) then
!!               u(i,j)   = -Y(i,j)*angrot
!!               v(i,j)   =  X(i,j)*angrot
!!            end if
!!         end do
!!         end do
!!         !!
!!         uc       = Rc*angrot!max speed
!!         dt       = CFL*dx/uc
!!         nt       = aint(2*pi/(angrot*dt))
!!         print*,angrot,dt,CFL
!!         dtheta   = dt*angrot
!!      end if
!!
!!      !if (.false.) then
!!      if (.true.) then
!!         stat1 = 'replace'
!!         stat2 = 'old'
!!         afile = 'test_out/all0.a'
!!
!!         print*,'******************************************'
!!         print*,'testing initial conditions:'
!!         print*,'check: ',afile
!!         print*,'******************************************'
!!         print*,''
!!
!!         call wrt_afile_2d(afile,u,    idm,jdm,stat1,1)
!!         call wrt_afile_2d(afile,v,    idm,jdm,stat2,2)
!!         call wrt_afile_2d(afile,scp2, idm,jdm,stat2,3)
!!         call wrt_afile_2d(afile,scp2i,idm,jdm,stat2,4)
!!         call wrt_afile_2d(afile,scuy, idm,jdm,stat2,5)
!!         call wrt_afile_2d(afile,scvx, idm,jdm,stat2,6)
!!      end if
!!
!!      !!output initial conditions for reading by matlab/python
!!      n  = 0
!!      stat1 = 'replace'
!!      stat2 = 'old'
!!      !!
!!      write(cts,'(i3.3)') n
!!      afile = outdir//outfile//cts//'.a'
!!      bfile = outdir//outfile//cts//'.b'
!!      call wrt_bfile(bfile,ii,jj,OPT)
!!      call wrt_afile_2d(afile,x,ii,jj,stat1,1)
!!      call wrt_afile_2d(afile,y,ii,jj,stat2,2)
!!      call wrt_afile_2d(afile,h,ii,jj,stat2,3)
!!      !!
!!      print*,'******************************************'
!!      print*,'finished writing initial conditions to:'
!!      print*,afile
!!      print*,bfile
!!      print*,'******************************************'
!!      print*,''
!!
!!      do n = 1,nt
!!      !do n = 1,1
!!         !print*,'timestep = ',n
!!         !!advect for one time step
!!         call waveadv_weno(h,u,v,dt,bc_opt)
!!
!!         !!output initial conditions for reading by matlab/python
!!         write(cts,'(i3.3)') n
!!         afile = outdir//outfile//cts//'.a'
!!         bfile = outdir//outfile//cts//'.b'
!!         call wrt_bfile(bfile,ii,jj,OPT)
!!         call wrt_afile_2d(afile,x,ii,jj,stat1,1)
!!         call wrt_afile_2d(afile,y,ii,jj,stat2,2)
!!         call wrt_afile_2d(afile,h,ii,jj,stat2,3)
!!      end do
!!
!!      print*,'******************************************'
!!      print*,'finished writing results to:'
!!      print*,afile
!!      print*,bfile
!!      print*,'******************************************'
!!      print*,''

      end program main
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine check_grid(GRID_OPT)

      use mod_gridinfo
      !use mod_common_wim
      use mod_file_utils

      implicit none
      integer,intent(in)   :: GRID_OPT
      character(len=80)    :: afile,bfile
      character(len=10)    :: stat1,stat2

      afile = 'out/wim_grid.a'
      bfile = 'out/wim_grid.b'
      stat1 = 'replace'
      stat2 = 'old'

 111  format(i2.2,a)  
 112  format(i3.3,a)  
 113  format(a)

      !! write .b file (description of stuff in .a file)
      open(unit=2,file=trim(bfile),status = 'replace')

      write(2,111)                                                      &
     &  5,          "       Number of records"
      write(2,112)                                                      &
     &  ii,          "      Record length in x direction (elements)"
      write(2,112)                                                      &
     &  jj,          "      Record length in y direction (elements)"
      write(2,111)                                                      &
     &  GRID_OPT,   "       Option number for initial conditions"
      write(2,113) ""

      write(2,113)    "Record number and name:"
      write(2,111)                                                      &
     &  1,    "       X"
      write(2,111)                                                      &
     &  2,    "       Y"
      write(2,111)                                                      &
     &  3,    "       scuy"
      write(2,111)                                                      &
     &  4,    "       scvx"
      write(2,111)                                                      &
     &  5,    "       scp2"
      write(2,111)                                                      &
     &  6,    "       scp2i"
      write(2,111)                                                      &
     &  7,    "       LANDMASK"

      close(2)

      !! dump arrays to .a file
      !! - X,Y,scuy,scvx,scp2,scp2i,LANDMASK
      call wrt_afile_2d(afile ,X          ,ii,jj  ,stat1,1)
      call wrt_afile_2d(afile ,Y          ,ii,jj  ,stat2,2)
      call wrt_afile_2d(afile ,scuy       ,ii,jj  ,stat2,3)
      call wrt_afile_2d(afile ,scvx       ,ii,jj  ,stat2,4)
      call wrt_afile_2d(afile ,scp2       ,ii,jj  ,stat2,5)
      call wrt_afile_2d(afile ,scp2i      ,ii,jj  ,stat2,6)
      call wrt_afile_2d(afile ,LANDMASK   ,ii,jj  ,stat2,7)

      end subroutine check_grid
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine check_init(GRID_OPT)

      use mod_gridinfo
      use mod_common_wim
      use mod_file_utils

      implicit none
      integer,intent(in)   :: GRID_OPT
      character(len=80)    :: afile,bfile
      character(len=10)    :: stat1,stat2

      afile = 'out/wim_init.a'
      bfile = 'out/wim_init.b'
      stat1 = 'replace'
      stat2 = 'old'

 111  format(i2.2,a)  
 112  format(i3.3,a)  
 113  format(a)

      !! write .b file (description of stuff in .a file)
      open(unit=2,file=trim(bfile),status = 'replace')

      write(2,111)                                                      &
     &  5,          "       Number of records"
      write(2,112)                                                      &
     &  ii,          "      Record length in x direction (elements)"
      write(2,112)                                                      &
     &  jj,          "      Record length in y direction (elements)"
      write(2,111)                                                      &
     &  GRID_OPT,   "       Option number for initial conditions"
      write(2,113) ""

      write(2,113)    "Record number and name:"
      write(2,111)                                                      &
     &  1,    "       icec"
      write(2,111)                                                      &
     &  2,    "       iceh"
      write(2,111)                                                      &
     &  3,    "       Dmax"
      write(2,111)                                                      &
     &  4,    "       Hs"
      write(2,111)                                                      &
     &  5,    "       Tp"
      write(2,111)                                                      &
     &  6,    "       mwd"

      close(2)

      !! dump arrays to .a file
      !! - X,Y,scuy,scvx,scp2,scp2i,LANDMASK
      call wrt_afile_2d(afile ,icec    ,ii,jj  ,stat1,1)
      call wrt_afile_2d(afile ,iceh    ,ii,jj  ,stat2,2)
      call wrt_afile_2d(afile ,dfloe   ,ii,jj  ,stat2,3)
      call wrt_afile_2d(afile ,Hs      ,ii,jj  ,stat2,4)
      call wrt_afile_2d(afile ,Tp      ,ii,jj  ,stat2,5)
      call wrt_afile_2d(afile ,mwd     ,ii,jj  ,stat2,6)

      end subroutine check_init
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
