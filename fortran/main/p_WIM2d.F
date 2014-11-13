      !! test_waveadv_weno.F
      !! Author: Timothy Williams
      !! Date:   20140821, 12:22:17 CEST
      program main

      use mod_gridinfo     !! grid info
      use mod_wim_prams    !! ice parameters
      use mod_wim_init     !! initial wave fields
      use mod_common_wim   !! initial ice fields
      use mod_adv_atten    !! advection routine

      implicit none

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!variables
      integer           :: SOLVER
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
      real,dimension(idm,jdm)    :: mom0,mom2,var_strain                &
     &   ,S_freq,atten_dim,damp_dim                                     &
     &   ,taux_om,tauy_om,tmp1,WAVE_MASK2
      real,dimension(idm,jdm,nw) ::                                     &
     &   ag_eff,ap_eff,wlng_ice,atten_nond,damping,disp_ratio

      real     :: PI,tmp
      real     :: dom_width,amin,amax,u_ref
      real     :: adv_dir,dt,dtheta,om,dom,guess
      real     :: kwtr,kice,modT,int_adm,argT,argR,F,Dave,c1d
      real,dimension(nw)   :: wt_om,wt_simp,T,wlng,ag,ap
      real     :: wt,Hs_test1,Hs_test2
      integer  :: nt,n,i,j,w,wth

      !!for breaking
      real     :: sig_strain,Pstrain,P_crit,wlng_crest,Dc
      real     :: ommin,ommax,om1,lam1,lam2
      integer  :: jcrest
      logical  :: BREAK_CRIT


      character(len=*),parameter :: outfile = "init"
      character(len=*),parameter :: outdir = "out/"
      !character(len=80)          :: afile,bfile,stat1,stat2

      real :: et,etime,tictoc(2) !!print time
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SOLVER   = 1
      PI       = 4.0*atan(1.0)
      et       = etime(tictoc)!!startup etime

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

      sdf_dir  = sdf_inc
      !sdf_dir  = 0.0
      !do w=1,nw
      !do wth=1,ndir
      !   where (WAVE_MASK.eq.1.0)
      !      sdf_dir(:,:,wth,w)   = sdf_inc(:,:,wth,w)
      !   end where
      !end do
      !end do

      WAVE_MASK2  = 0.0
!     if (STEADY.eq.1) then
!        if (GRID_OPT.eq.1) then
!           WAVE_MASK2(1:3,:) = 1.0
!        end if
!        print*,'steady test',maxval(WAVE_MASK2),sum(WAVE_MASK2)
!     end if

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
      T        = 1.0/freq_vec
      wlng     = gravity/(2*PI*freq_vec**2)
      ap       = sqrt(gravity*wlng/(2*PI))   ! Phase speed
      ag       = ap/2.0                      ! Group speed
      !print*,'T,wlng,ag,ap',T,wlng,ag,ap

      ag_eff      = 0.0
      ap_eff      = 0.0
      wlng_ice    = 0.0
      disp_ratio  = 1.0
      atten_nond  = 0.0
      damping     = 0.0

      om    = 2*PI*freq_vec(1)

      do w=1,nw
         
         do j=1,jj
         do i=1,ii
         if (ICE_MASK(i,j).eq.1.0) then

            om = 2*PI*freq_vec(w)
            if (w.eq.1) then
               guess = om**2/gravity!guess for root
            else
               guess = 2*PI/wlng_ice(i,j,w-1)
!              if ((i.eq.24).and.(j.eq.1)) then
!                 print*,'previous wavelengths/nos'
!                 do wth=1,w-1
!                    print*,wlng_ice(i,j,wth),2*PI/wlng_ice(i,j,wth)
!                 end do
!              end if
            end if

            call RTparam_outer(damping(i,j,w),kice,kwtr,int_adm,        &
     &                  atten_nond(i,j,w),modT,argR,argT,               &
     &                  iceh(i,j),om,guess)

            disp_ratio(i,j,w) = kice/kwtr*modT
            wlng_ice(i,j,w)   = 2*PI/kice

            !!check attenuation coeff's, wavenumbers
!           if ((i.eq.24).and.(j.eq.1)) then
!              print*,'***********************************************'
!              print*,'check outputs from RTparam_outer, (i,j,w)=',i,j,w
!              print*,'atten',atten_nond(i,j,w),damping(i,j,w)
!              print*,'ki,kw,2pi/wlng_wtr',kice,kwtr,2*PI/wlng(w)
!              print*,'lam,|T|,disp_rat',                               &
!    &            wlng_ice(i,j,w),modT,disp_ratio(i,j,w)
!              print*,'argRT_s',argR,argT,int_adm
!              print*,'***********************************************'
!              print*,' '
!           end if

            if (USE_ICE_VEL.eq.0) then
               !!use water group and phase velocity
               ag_eff(i,j,w)  = ag(w)
               ap_eff(i,j,w)  = ap(w)
            end if

         else!end check for ice

            !!no ice
            wlng_ice(i,j,w)   = wlng(w)
            ag_eff(i,j,w)     = ag(w)
            ap_eff(i,j,w)     = ap(w)

         end if!check for no ice

         end do!i - spatial (x) loop
         end do!j - spatial (y) loop
      end do!w - spectral loop

      if (DO_ATTEN.eq.0) then
         atten_nond  = 0.0
         damping     = 0.0
      end if

      amax        = maxval(ag_eff(:,:,1 ))!lowest freq is fastest
      amin        = minval(ag_eff(:,:,nw))!highest freq is slowest
      u_ref       = amin+.7*(amax-amin);
      dt          = CFL*dx/amax
      dom_width   = maxval(X)-minval(X)
      nt          = aint( dom_width/u_ref/dt )
      
      print*,' '
      print*,'******************************************************'
      print*,'CFL number:                       '  ,CFL
      print*,'Spatial resolution (km):          '  ,dx/1.0e3,dy/1.0e3
      print*,'Maximum wave group velocity (m/s):'  ,amax
      print*,'Width of domain (km):             '  ,dom_width/1.0e3
      print*,'Number of time steps:             '  ,nt
      print*,'Time step (s):                    '  ,dt
      print*,'Time interval (h):                '  ,(nt*dt)/60.0/60.0
      print*,'Number of wave frequencies:       '  ,nw
      print*,'Minimum period (s):               '  ,1.0/freq_vec(nw)
      print*,'Maximum period (s):               '  ,1.0/freq_vec(1)
      print*,'Number of wave directions:        '  ,ndir
      print*,'Directional resolution (degrees): '  ,360.0/(1.0*ndir)
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

         !! top-up waves (with directions that correspond to incident
         !! waves) in appropriate area if STEADY==1
         !! (steady-state solution);
         do w=1,nw
         do wth=1,n_wavdir
            adv_dir  = -PI/180.0*(wavdir(wth)+90.0)
            if (cos(adv_dir).ge.0.0) then
               where (WAVE_MASK2.gt.0.0)
                  !sdf_dir(i,j,J_STEADY,:) = sdf_inc(i,j,J_STEADY,:)
                   sdf_dir(:,:,wth,w)   = sdf_inc(:,:,wth,w)
               end where
            end if
         end do
         end do

         i        = 2
         j        = 25
         Hs_test1 = 0.0
         Hs_test2 = 0.0
         do w=1,nw
            if (n_wavdir.eq.1) then
               Hs_test1 = Hs_test1+wt_om(w)*sdf_dir(i,j,1,w)
               Hs_test2 = Hs_test2+wt_om(w)*sdf_inc(i,j,1,w)
            else
               do wth=1,n_wavdir
                  wt       = (2*PI/(1.0*n_wavdir))*wt_om(w)
                  Hs_test1 = Hs_test1+wt*sdf_dir(i,j,wth,w)
                  Hs_test2 = Hs_test2+wt*sdf_inc(i,j,wth,w)
                  !!
                  adv_dir  = -PI/180.0*(wavdir(wth)+90.0)
!                 print*,sdf_dir(i,j,wth,w),sdf_inc(i,j,wth,w)          &
!    &                  ,adv_dir,(cos(adv_dir).ge.0.0)
               end do
            end if
         end do
         Hs_test1 = 4*sqrt(Hs_test1)
         Hs_test2 = 4*sqrt(Hs_test2)
         print*,'Hs_test inc',Hs_test1,Hs_test2

         do w=1,nw

            !! CALC DIMENSIONAL ATTEN COEFF;
            atten_dim   = 0.0
            damp_dim    = 0.0
            do j = 1,jj
            do i = 1,ii

               if ((ICE_MASK(i,j).eq.1.0).and.(DO_ATTEN.eq.1)) then
                  call floe_scaling(Dave,dfloe(i,j))

                  !! get expected no of floes met per unit
                  !!  distance if travelling in a line;
                  if (dfloe(i,j).lt.200.0) then
                     c1d = icec(i,j)/Dave
                     !! floes per unit length;
                  else
                     c1d = icec(i,j)/dfloe(i,j)
                     !! uniform lengths
                  end if

                  !! ENERGY attenuation coeff;
                  atten_dim(i,j) = atten_nond(i,j,w)*c1d      !!scattering
                  damp_dim(i,j)  = 2*damping(i,j,w)*icec(i,j) !!damping
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
               call adv_atten_simple(sdf_dir(:,:,:,w)                   &
     &            ,S_freq,taux_om,tauy_om                               &
     &            ,atten_dim,damp_dim,ag_eff(:,:,w),dt)
               !print*,'Spectrum at',itest,jtest
               !print*,S_freq(itest,jtest)
               !print*,sdf_dir(itest,jtest,:,w)
            elseif (SOLVER.eq.1) then
               !! same as SOLVER==0, but scattered energy
               !! is distributed isotropically
               call adv_atten_isotropic(sdf_dir(:,:,:,w)                &
     &            ,S_freq,taux_om,tauy_om                               &
     &            ,atten_dim,damp_dim,ag_eff(:,:,w),dt)
            end if

            !! integrate stress densities over frequency
            !! TODO: check if this is correct for ice-covered water
            tmp1  = rhowtr*gravity*taux_om/ap_eff(:,:,w)   !![Pa*s]
            tau_x = tau_x+wt_om(w)*tmp1                     !![Pa]
            tmp1  = rhowtr*gravity*tauy_om/ap_eff(:,:,w)  !![Pa*s]
            tau_y = tau_y+wt_om(w)*tmp1                     !![Pa]

            !! INTEGRALS FOR BREAKING PROB:
            do j=1,jj
            do i=1,ii

               !! convert from water amp's to ice amp's;
               F     = disp_ratio(i,j,w)
               kice  = 2*PI/wlng_ice(i,j,w)

               !! SPECTRAL MOMENTS;
               !!take abs as small errors can make S_freq negative
               om          = 2*PI*freq_vec(w)
               tmp         = wt_om(w)*S_freq(i,j)*F**2
               mom0(i,j)   = mom0(i,j)+abs( tmp )!!variance of displacement
               tmp         = wt_om(w)*om**2*S_freq(i,j)*F**2
               mom2(i,j)   = mom2(i,j)+abs( tmp )!!variance of speed

               if (ICE_MASK(i,j).eq.1.0) then
                  !! VARIANCE OF STRAIN;
                  tmp   = F*kice**2*iceh(i,j)/2.0    !!strain conversion factor
                  tmp   = wt_om(w)*S_freq(i,j)*tmp**2 !!strain density
                  var_strain(i,j)   = var_strain(i,j)+abs( tmp )
               end if
            end do!! i, end spatial loop x;
            end do!! j, end spatial loop y;

         end do!! end spectral loop;

         Hs = 4*sqrt(mom0)
         Tp = 0.0
         where (mom2.gt.0.0)
            Tp = 2*PI*sqrt(mom0/mom2)
         end where

         !! FINALLY DO FLOE BREAKING;
         do j=1,jj
         do i=1,ii
            if ((ICE_MASK(i,j).eq.1.0).and.(mom0(i,j).ge.0.0)) then
               !! only try breaking if ice is present
               !!  & some waves have arrived;
      
               !! significant strain amp
               sig_strain  = 2*sqrt(var_strain(i,j))
      
               !! probability of critical strain
               !! being exceeded from Rayleigh distribution;
               Pstrain  = exp( -epsc**2/(2*var_strain(i,j)) )
               P_crit   = (1-DO_BREAKING)+exp(-1.0)!!this is critical prob if monochromatic wave
      
               !! FLOE BREAKING:
               BREAK_CRIT  = ( Pstrain.ge.P_crit )!!breaks if larger than this
      
               if (BREAK_CRIT) then
                  !! use crest period to work out wavelength
                  !! - interpolate using pre-calc'd wavelengths
                  !! - half this is max poss floe length;
                  om    = 2*pi/Tp(i,j)
                  ommin = 2*PI*freq_vec(1)
                  ommax = 2*PI*freq_vec(nw)
                  if (om.le.ommin) then
                     wlng_crest  = wlng_ice(i,j,1)
                  elseif (om.ge.ommax) then
                     wlng_crest  = wlng_ice(i,j,nw)
                  else
                     jcrest      = floor((om-ommin+dom)/dom)
                     om1         = 2*PI*freq_vec(jcrest)
                     lam1        = wlng_ice(i,j,jcrest)
                     lam2        = wlng_ice(i,j,jcrest+1)
                     wlng_crest  = lam1+(om-om1)*(lam2-lam1)/dom
                  end if
      
                  !!change dfloe if wavelength/2 is shorter than it
                  !!(keeping it bigger than Dmin)
                  Dc          = max(Dmin,wlng_crest/2.0)
                  dfloe(i,j)  = min(Dc,dfloe(i,j))

               end if!! end breaking action;
      
            elseif (WTR_MASK(i,j).eq.1.0) then!! only water present
               dfloe(i,j)  = 0.0
            end if
            
         end do!! end spatial loop j in y;
         end do!!end spatial loop i in x;
         
         !!dump fields to binary
         !print*,'dfloe  lims',maxval(dfloe),minval(dfloe)
         !print*,'Hs     lims',maxval(Hs   ),minval(Hs   )
         call check_prog(GRID_OPT,n)

      end do!! n - finish time stepping

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

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine check_prog(GRID_OPT,n)

      use mod_gridinfo
      use mod_common_wim
      use mod_file_utils

      implicit none
      integer,intent(in)   :: GRID_OPT,n
      character(len=80)    :: afile,bfile
      character(len=10)    :: stat1,stat2
      character(len=3)     :: cts

      write(cts,'(i3.3)') n
      afile = 'out/wim_prog'//cts//'.a'
      bfile = 'out/wim_prog'//cts//'.b'
      stat1 = 'replace'
      stat2 = 'old'

 111  format(i2.2,a)  
 112  format(i3.3,a)  
 113  format(a)

      !! write .b file (description of stuff in .a file)
      open(unit=2,file=trim(bfile),status = 'replace')

      write(2,111)                                                      &
     &  4,          "       Number of records"
      write(2,112)                                                      &
     &  ii,          "      Record length in x direction (elements)"
      write(2,112)                                                      &
     &  jj,          "      Record length in y direction (elements)"
      write(2,111)                                                      &
     &  GRID_OPT,   "       Option number for initial conditions"
      write(2,113) ""

      write(2,113)    "Record number and name:"
      write(2,111)                                                      &
     &  1,    "       Dmax"
      write(2,111)                                                      &
     &  2,    "       tau_x"
      write(2,111)                                                      &
     &  3,    "       tau_y"
      write(2,111)                                                      &
     &  4,    "       Hs"
!     write(2,111)                                                      &
!    &  5,    "       Tp"
!     write(2,111)                                                      &
!    &  6,    "       mwd"

      close(2)

      !! dump arrays to .a file
      !! - X,Y,scuy,scvx,scp2,scp2i,LANDMASK
      call wrt_afile_2d(afile ,dfloe   ,ii,jj  ,stat1,1)
      call wrt_afile_2d(afile ,tau_x   ,ii,jj  ,stat2,2)
      call wrt_afile_2d(afile ,tau_y   ,ii,jj  ,stat2,3)
      call wrt_afile_2d(afile ,Hs      ,ii,jj  ,stat2,4)
!     call wrt_afile_2d(afile ,Tp      ,ii,jj  ,stat2,5)
!     call wrt_afile_2d(afile ,mwd     ,ii,jj  ,stat2,6)

      end subroutine check_prog
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine floe_scaling(Dave,dmax)

      use mod_wim_prams, only: fragility,xi,Dmin
      implicit none

      real,    intent(out) :: Dave
      real,    intent(in)  :: dmax
      integer              :: m,mm 
      real                 :: n,nsum,nd,ndsum,r,dfac

      mm    = 0  
      r     = dmax/dmin
      do while ( r .gt. xi ) 
         r  = r/xi 
         mm = mm+1 
      end do

      if ( mm .gt. 0 ) then 
      n     = 0.0
      nsum  = 0.0
      nd    = 0.0
      ndsum = 0.0

         do m = 0,1,mm
            n     = (1.0-fragility)*(fragility*xi**2)**m
            nd    = n/(xi**m)
            nsum  = nsum +n
            ndsum = ndsum+nd
            dfac  = ndsum/nsum
         end do
         Dave  = dfac*dmax

      else
         Dave = dmin
      end if

      end subroutine floe_scaling
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
