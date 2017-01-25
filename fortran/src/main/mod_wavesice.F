      !mod_wavesice.F
      !Author: Timothy Williams
      !Date: 20141119, 05:22:32 CET
      module mod_wavesice

      use mod_gridinfo     !! grid info
      use mod_common_wim   !! initial ice fields
      use mod_year_info

      implicit none
      private

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!variables

      !!TODO: make these inputs
      !integer,parameter :: STEADY      = 0
      !integer,parameter :: ADV_OPT     = 0
         !!0, not periodic; 2, periodic in y only; 1, periodic in both x,y;

      integer :: REF_HS_ICE  = 0 !1 = in ice, give Hs for displacement of the ice
      integer :: USE_ICE_VEL = 0 !0 = approx ice group vel by water group vel
      !!
      integer,parameter :: nw    = n_wave_freq
      integer,parameter :: ndir  = n_wavdir
      integer,parameter :: ne    = 2

      !!common to all subroutines
      real,dimension(idm,jdm)    :: WAVE_MASK2
      real,dimension(idm,jdm,nw) ::                                     &
     &   ag_eff,ap_eff,wlng_ice,atten_nond,damping,disp_ratio
      real,dimension(nw)   :: wt_om

      !!needed by outer routine
      real     :: amin,amax,dt
      
      logical  :: DUMP_DIAG

      type break_info
         real  :: conc,thick
     &            ,mom0,mom2,var_strain
     &            ,dmax,broken
      end type
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      public   :: wim_init,wim_step,do_wim_integrals
     &            ,check_grid,check_init,check_prog
     &            ,amin,amax,dt,WAVE_MASK2
     &            ,USE_ICE_VEL,REF_HS_ICE
     &            ,DUMP_DIAG
      !,calc_mwd

      contains

      !subroutine wim_init(STEADY,DO_ATTEN,CFL)
      subroutine wim_init()
      !! set initial waves
      !! set attenuation coefficients
      !! - need to do this each time ice is updated

      use mod_gridinfo     !! grid info
      use mod_common_wim,  only: ICE_MASK,iceh           !! initial ice fields

      implicit none

!     integer, intent(in) :: STEADY,DO_ATTEN
!     real,    intent(in) :: CFL
      !!
      real     :: PI,tmp
      real     :: kwtr,kice,modT,int_adm,argT,argR
      real,dimension(nw)   :: wt_simp,wlng,ag,ap
      integer  :: i,j,w
      real     :: adv_dir,dtheta,om,dom,guess,inputs(2)
!     integer  :: wth
!     real     :: wt,Hs_test1,Hs_test2

      PI       = 4.0*atan(1.0)

      !! TODO: pass in WAVE_MASK@ as an input??
      WAVE_MASK2  = 0.0
      if (STEADY.eq.1) then
         WAVE_MASK2(1:3,:) = 1.0
         print*,'steady test',maxval(WAVE_MASK2),sum(WAVE_MASK2)
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
            w           = w+2
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
      wlng     = wim_info%gravity/(2*PI*freq_vec**2)
      ap       = sqrt(wim_info%gravity*wlng/(2*PI))   ! Phase speed
      ag       = ap/2.0                      ! Group speed
      !print*,'period,wlng,ag,ap',1./freq_vec,wlng,ag,ap

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
               guess = om**2/wim_info%gravity!guess for root
            else
               guess = 2*PI/wlng_ice(i,j,w-1)
!              if ((i.eq.itest).and.(j.eq.jtest) then
!                 print*,'previous wavelengths/nos'
!                 do wth=1,w-1
!                    print*,wlng_ice(i,j,wth),2*PI/wlng_ice(i,j,wth)
!                 end do
!              end if
            end if

            inputs(1)   = iceh(i,j)
            inputs(2)   = om

            call RTparam_outer(damping(i,j,w),kice,kwtr,int_adm,        &
     &                  atten_nond(i,j,w),modT,argR,argT,               &
     &                  inputs,guess)

            disp_ratio(i,j,w) = kice/kwtr*modT
            wlng_ice(i,j,w)   = 2*PI/kice

!           !!check attenuation coeff's, wavenumbers
!           if ((i.eq.itest).and.(j.eq.jtest)) then
!              print*,'***********************************************'
!              print*,'check outputs from RTparam_outer, (i,j,w)=',i,j,w
!              print*,'T,h,c',1.0/freq_vec(w),iceh(i,j),icec(i,j)
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

      amax  = maxval(ag_eff(:,:,1 ))!lowest freq is fastest
      amin  = minval(ag_eff(:,:,nw))!highest freq is slowest
      dt    = CFL*dx/amax

      end subroutine wim_init
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine wim_step(nmesh_e,mesh_e)
      !! do 1 time step

      use mod_gridinfo     !! grid info
      use mod_common_wim   !! initial ice fields
      use mod_adv_atten    !! advection routine
      use mod_interp, only: pwl_interp_2d_real

      implicit none

      type(break_info)  :: brk_info
      !integer,parameter :: n_wim_step_prams  = 7
      !integer,dimension(n_wim_step_prams),intent(in) :: wim_step_prams  &
         !SCATMOD,ADV_DIM,ADV_OPT,BRK_OPT,STEADY,DO_ATTEN,FSD_OPT
      integer,intent(in)   :: nmesh_e

      real,dimension(nmesh_e,nmesh_vars),intent(inout),optional         &
     &   :: mesh_e
         ! this is [x_e,y_e,conc_e,thick_e,Dmax_var_e,broken_e]
         ! -1st 4 are unchanged; last 2 get updated

!     integer  :: SCATMOD,ADV_DIM,ADV_OPT                               &
!    &            ,BRK_OPT,STEADY,DO_ATTEN,FSD_OPT

      real,dimension(idm,jdm)    :: mom0,mom2,var_strain                &
     &   ,S_freq,atten_dim,damp_dim                                     &
     &   ,taux_om,tauy_om,tmp1                                          &
     &   ,mom0w,mom2w
      integer  :: i,j,w,ibrk

      !!to top up waves
      integer  :: wth
      real  :: adv_dir,PI

!     !!for testing
!     real     :: wt,Hs_test1,Hs_test2
      real  :: E_tot

      !!for breaking
      real,dimension(ii,jj)   :: Dave
      real     :: sig_strain,Pstrain,P_crit0,P_crit,wlng_crest,Dc
      real     :: om,ommin,ommax,om1,lam1,lam2,dom
      real     :: F,c1d,kice,tmp
      integer  :: jcrest,fid_diag
      logical  :: BREAK_CRIT

      logical                 :: INTERP_MESH,crit_e,test_ij
      real                    :: wim_break_damage,Tp_e,lam_e
      real,dimension(nmesh_e) :: x_e,y_e                    !centres of elements
      real,dimension(nmesh_e) :: mom0_e,mom2_e,var_strain_e !interpolated values onto elements
      real,dimension(nmesh_e) :: conc_e,thick_e,Dmax_var_e
     &                              ,broken_e,dmax_e
      real     :: kwtr,modT,int_adm,argT,argR,guess,inputs(2)
      character(len=800) :: diag_filename
      !
      INTERP_MESH = .false.
      if (present(mesh_e)) then
         ! check if there is any ice
         INTERP_MESH = maxval(mesh_e).gt.0.
      end if

      
      if (DUMP_DIAG) then
         diag_filename  = trim(outdir)//'/diagnostics/local/'//
     &      'WIMdiagnostics_local'//model_year_info%cdate//
     &      'T'//model_year_info%chour//model_year_info%cminute//
     &      model_year_info%csecond//'Z.txt'
         fid_diag = 91
         open(unit=fid_diag,file=trim(diag_filename),status = 'replace')
         write(fid_diag,'(a,a)'),model_year_info%cdate,' # date'
         write(fid_diag,'(a,a)'),model_year_info%ctime,' # time'
         write(fid_diag,'(i6.6,a)'),model_day,' # model day'
         write(fid_diag,'(f13.5,a)'),model_seconds,' # model second'
         write(fid_diag,'(i6.6,a)'),itest,' # itest'
         write(fid_diag,'(i6.6,a)'),jtest,' # jtest'
         write(fid_diag,'(f3.1,a)'),ICE_MASK(itest,jtest),' # ICE_MASK'
         write(fid_diag,'(a)'),' '
      end if


      !wim_step_prams   = (/
      !SCATMOD,ADV_DIM,ADV_OPT,BRK_OPT,STEADY,DO_ATTEN,FSD_OPT
      !/)
!     SCATMOD  = wim_step_prams(1)
!     ADV_DIM  = wim_step_prams(2)
!     ADV_OPT  = wim_step_prams(3)
!     BRK_OPT  = wim_step_prams(4)
!     STEADY   = wim_step_prams(5)
!     DO_ATTEN = wim_step_prams(6)
!     FSD_OPT  = wim_step_prams(7)

      PI    = 4*atan(1.0)
      dom   = 2*PI*(freq_vec(nw)-freq_vec(1))/(nw-1.0)
      E_tot = 0.

      !! spectral moments;
      mom0  = 0.0
      mom2  = 0.0
      mom0w = 0.0
      mom2w = 0.0

      !! wave stresses;
      tau_x = 0.0
      tau_y = 0.0

      !! variance of strain;
      var_strain     = 0.0

      !! top-up waves (with directions that correspond to incident
      !! waves) in appropriate area if STEADY==1
      !! (steady-state solution);
      if (STEADY.eq.1) then
         do w=1,nw
         do wth=1,n_wavdir
            adv_dir  = -PI/180.0*(wavdir(wth)+90.0)
            if (cos(adv_dir).ge.0.0) then
               where (WAVE_MASK2.gt.0.0)
                   sdf_dir(:,:,wth,w)   = sdf_inc(:,:,wth,w)
               end where
               !slower way to do above test:
               !do j=1,jj
               !do i=1,ii
               !if (WAVE_MASK2(i,j)==1.0) then
               !    sdf_dir(i,j,wth,w)   = sdf_inc(i,j,wth,w)
               !end if
               !end do
               !end do
            end if
         end do
         end do
      end if
!     !!test incident waves
!     i        = 2
!     j        = 25
!     Hs_test1 = 0.0
!     Hs_test2 = 0.0
!     do w=1,nw
!        if (n_wavdir.eq.1) then
!           Hs_test1 = Hs_test1+wt_om(w)*sdf_dir(i,j,1,w)
!           Hs_test2 = Hs_test2+wt_om(w)*sdf_inc(i,j,1,w)
!        else
!           do wth=1,n_wavdir
!              wt       = (2*PI/(1.0*n_wavdir))*wt_om(w)
!              Hs_test1 = Hs_test1+wt*sdf_dir(i,j,wth,w)
!              Hs_test2 = Hs_test2+wt*sdf_inc(i,j,wth,w)
!              !!
!              adv_dir  = -PI/180.0*(wavdir(wth)+90.0)
!              !print*,sdf_dir(i,j,wth,w),sdf_inc(i,j,wth,w)          &
!    &         !     ,adv_dir,(cos(adv_dir).ge.0.0)
!           end do
!        end if
!     end do
!     Hs_test1 = 4*sqrt(Hs_test1)
!     Hs_test2 = 4*sqrt(Hs_test2)
!     print*,'Hs_test inc',Hs_test1,Hs_test2

      !!get <D>
      Dave  = 0.
      do j = 1,jj
      do i = 1,ii

         if ((ICE_MASK(i,j).eq.1.0).and.(DO_ATTEN.eq.1)) then

            !! get expected no of floes met per unit
            !!  distance if travelling in a line;
            if (dfloe(i,j).le.200.0) then
               if (FSD_OPT.eq.0) then
                  !! renormalisation group method (Dumont et al, 2011)
                  call floe_scaling(Dave(i,j),dfloe(i,j),1)
               elseif (FSD_OPT.eq.1) then
                  !! smoothed power law method (Montiel et al, 2015)
                  call floe_scaling_smooth(Dave(i,j),dfloe(i,j),1)
               end if
            else
               Dave(i,j)   = dfloe(i,j)
               !! uniform lengths
            end if

            test_ij  = (i.eq.itest).and.(j.eq.jtest) 
            if (DUMP_DIAG.and.test_ij) then
               write(fid_diag,'(a)'),'Ice info: pre-breaking'
               write(fid_diag,'(f8.4,a)'),icec(i,j)
     &         ,' # conc'
               write(fid_diag,'(f8.4,a)'),iceh(i,j)
     &         ,' # h, m'
               write(fid_diag,'(f8.4,a)'),Dave(i,j)
     &         ,' # D_av, m'
               write(fid_diag,'(f9.5,a)'),dfloe(i,j)
     &         ,' # D_max, m'
            end if
         end if
      end do
      end do

      if (DUMP_DIAG.and.ICE_MASK(itest,jtest).eq.1.) then
         write(fid_diag,'(a)'),' '
         write(fid_diag,'(a)'),'# period, s | atten_dim, m^{-1}'//
     &                         '| damp_dim, m^{-1}'
      end if
      do w=1,nw

         !! CALC DIMENSIONAL ATTEN COEFF;
         atten_dim   = 0.0
         damp_dim    = 0.0
         do j = 1,jj
         do i = 1,ii

            if ((ICE_MASK(i,j).eq.1.0).and.(DO_ATTEN.eq.1)) then

               !! get expected no of floes met per unit
               !!  distance if travelling in a line;
               c1d = icec(i,j)/Dave(i,j)!! floes per unit length;

               !! ENERGY attenuation coeff;
               atten_dim(i,j) = atten_nond(i,j,w)*c1d      !!scattering
               damp_dim(i,j)  = 2*damping(i,j,w)*icec(i,j) !!damping

               test_ij  = (i.eq.itest).and.(j.eq.jtest) 
               if (DUMP_DIAG.and.test_ij) then
                  write(fid_diag,'(f8.4,a,e13.6,a,e13.6)')
     &                  ,1./freq_vec(w),' | '
     &                  ,atten_dim(i,j),' | '
     &                  ,damp_dim(i,j)
               end if

!                if ((i.eq.itest).and.(j.eq.jtest)) then
!                   print*,' '
!                   print*,'Hs (pre) (m)  = ',Hs(i,j)
!                   print*,'Tp (pre) (m)  = ',Tp(i,j)
!                   print*,'Dmax (m)      = ',dfloe(i,j)
!                   print*,'Dave (m)      = ',Dave
!                   print*,'c1d           = ',c1d
!                   print*,'q_scat (1/m)  = ',atten_dim(i,j)
!                   print*,'q_abs  (1/m)  = ',damp_dim(i,j)
!                   print*,' '
!                end if
               
            end if
         end do!! j
         end do!! i

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!stop if ndir/SCATMOD combo not OK
         if (ndir==1) then
            if (SCATMOD.ne.0) then
               print*,'warning: not enough directions'
               print*,'  for SCATMOD option:'
               print*,'  ndir   = ',ndir
               print*,'  SCATMOD = ',SCATMOD
               STOP
            end if
         end if
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !! advect 1 frequency
         call advect_1freq(sdf_dir(:,:,:,w),S_freq                      &
     &      ,ag_eff(:,:,w),dt,ADV_DIM,ADV_OPT)

         if (DO_ATTEN.eq.1) then

            !! *do attenuation/scattering
            !! *also calc momentum flux due to this process

            if (SCATMOD.eq.0) then
               !! Simple attenuation scheme
               !! - doesn't conserve scattered energy
               !print*,'Simple atten'
               call atten_simple(sdf_dir(:,:,:,w)                       &
     &            ,S_freq,taux_om,tauy_om                               &
     &            ,atten_dim,damp_dim,ag_eff(:,:,w),dt)

            elseif (SCATMOD.eq.1) then
               !! scattered energy is distributed isotropically
               !print*,'Isotropic atten'
               call atten_isotropic(sdf_dir(:,:,:,w)                    &
     &            ,S_freq,taux_om,tauy_om                               &
     &            ,atten_dim,damp_dim,ag_eff(:,:,w),dt)

            end if

            !! integrate stress densities over frequency
            !! TODO: check if this is correct for ice-covered water
            tmp1  = wim_info%rhowtr*wim_info%gravity
     &                  *taux_om/ap_eff(:,:,w)  !![Pa*s]
            tau_x = tau_x+wt_om(w)*tmp1                   !![Pa]
            tmp1  = wim_info%rhowtr*wim_info%gravity
     &                  *tauy_om/ap_eff(:,:,w)  !![Pa*s]
            tau_y = tau_y+wt_om(w)*tmp1                   !![Pa]

         end if!!DO_ATTEN


         !! INTEGRALS FOR BREAKING PROB:
         do j=1,jj
         do i=1,ii

            !! convert from water amp's to ice amp's;
            F     = disp_ratio(i,j,w)
            kice  = 2*PI/wlng_ice(i,j,w)

            !! SPECTRAL MOMENTS;
            !!take abs as small errors can make S_freq negative
            om          = 2*PI*freq_vec(w)
            tmp         = wt_om(w)*S_freq(i,j)
            mom0w(i,j)  = mom0w(i,j)+ abs( tmp )    !!variance of displacement (water)
            mom0(i,j)   = mom0(i,j) + abs( tmp*F**2 )!!variance of displacement (ice)
            tmp         = wt_om(w)*om**2*S_freq(i,j)
            mom2w(i,j)  = mom2w(i,j) + abs( tmp )    !!variance of speed (water)
            mom2(i,j)   = mom2(i,j)  + abs( tmp*F**2 )!!variance of speed (ice)

            if (ICE_MASK(i,j).eq.1.0) then
               !! VARIANCE OF STRAIN;
               tmp   = F*kice**2*iceh(i,j)/2.0     !!strain conversion factor
               tmp   = wt_om(w)*S_freq(i,j)*tmp**2 !!strain density
               var_strain(i,j)   = var_strain(i,j)+abs( tmp )
            end if
         end do!! i, end spatial loop x;
         end do!! j, end spatial loop y;

      end do!! w, end spectral loop;

      if (REF_HS_ICE.eq.1) then
         Hs = 4*sqrt(mom0)
         Tp = 0.0
         where (mom2.gt.0.0)
            Tp    = 2*PI*sqrt(mom0/mom2)
            !mwd   = -90. -180./PI*(mom_dir/mom0)!mwd (waves-from in degrees)
         end where
      else
         Hs = 4*sqrt(mom0w)
         Tp = 0.0
         where (mom2w.gt.0.0)
            Tp = 2*PI*sqrt(mom0w/mom2w)
            !mwd   = -90. -180./PI*(mom_dir_w/mom0w)!mwd (waves-from in degrees)
         end where
      end if
      call calc_mwd(REF_HS_ICE)

      if (DUMP_DIAG) then
         write(fid_diag,'(a)'),' '
         write(fid_diag,'(e13.6,a)')
     &   ,mom0w(itest,jtest),' # mom0w, m^2'
         write(fid_diag,'(e13.6,a)')
     &   ,mom2w(itest,jtest),' # mom2w, m^2/s^2'
         write(fid_diag,'(e13.6,a)')
     &   ,mom0(itest,jtest),' # mom0, m^2'
         write(fid_diag,'(e13.6,a)')
     &   ,mom2(itest,jtest),' # mom2, m^2/s^2'
         write(fid_diag,'(f8.4,a)')
     &   ,Hs(itest,jtest),' # Hs, m'
         write(fid_diag,'(f8.4,a)')
     &   ,Tp(itest,jtest),' # Tp, s'
         write(fid_diag,'(f10.4,a)')
     &   ,mwd(itest,jtest),' # mwd, deg'
         write(fid_diag,'(e13.6,a)')
     &   ,tau_x(itest,jtest),' # tau_x, Pa'
         write(fid_diag,'(e13.6,a)')
     &   ,tau_y(itest,jtest),' # tau_y, Pa'
         write(fid_diag,'(a)'),' '
      end if


      if ((STEADY.eq.0).and.(BRK_OPT.eq.0)) then
         !! STEADY==0: not steady-state so E should be conserved until
         !! it starts to exit the domain?
         !! TODO Try with periodicity?

         !! test conservation of total energy 
         E_tot = sum(Hs*Hs)
         print*,'********************************************'
         print*,'Total energy [m^2] = ',E_tot
         !print*,'alp nd/d',maxval(atten_nond),maxval(atten_dim)
         print*,'c1d,dave,dmax',c1d,Dave,maxval(dfloe)
         print*,'********************************************'
      end if


      !! FINALLY DO FLOE BREAKING (BRK_OPT>0);
      do j=1,jj
      do i=1,ii
         if ((BRK_OPT.gt.0)                                             &
     &       .and.(ICE_MASK(i,j).eq.1.0).and.(mom0(i,j).ge.0.0)) then
            !! only try breaking if ice is present
            !!  & some waves have arrived;
            brk_info%conc        = icec(i,j)
            brk_info%thick       = iceh(i,j)
            brk_info%mom0        = mom0(i,j)
            brk_info%mom2        = mom2(i,j)
            brk_info%var_strain  = var_strain(i,j)
            brk_info%dmax        = dfloe(i,j)
            brk_info%broken      = 0.

            call do_breaking(brk_info)
            dfloe(i,j)  = brk_info%dmax
   
         elseif (WTR_MASK(i,j).eq.1.0) then!! only water present
            dfloe(i,j)  = 0.0
         end if
         

         test_ij  = (i.eq.itest).and.(j.eq.jtest) 
         if (DUMP_DIAG.and.test_ij) then
            write(fid_diag,'(a)'),'Ice info: post-breaking'
            write(fid_diag,'(f9.5,a)'),dfloe(i,j)
     &      ,' # D_max, m'
         end if

      end do!! end spatial loop j in y;
      end do!!end spatial loop i in x;
      close(fid_diag)


      if ((INTERP_MESH).and.(BRK_OPT>0)) then
         x_e         = mesh_e(:,1)
         y_e         = mesh_e(:,2)
         conc_e      = mesh_e(:,3)
         thick_e     = mesh_e(:,4)
         Dmax_var_e  = mesh_e(:,5)
         broken_e    = mesh_e(:,6)

         call pwl_interp_2d_real(ii,jj,X(:,1),Y(1,:),mom0
     &      ,nmesh_e,x_e,y_e,mom0_e)
         call pwl_interp_2d_real(ii,jj,X(:,1),Y(1,:),mom2
     &         ,nmesh_e,x_e,y_e,mom2_e)
         call pwl_interp_2d_real(ii,jj,X(:,1),Y(1,:),var_strain
     &      ,nmesh_e,x_e,y_e,var_strain_e)

         dmax_e   = 0.
         if (Dmax_varname == 'Nfloes') then
            !Dmax stored as Nfloes = c/Dmax^2
            where (Dmax_var_e.gt.0)
               dmax_e   = sqrt(conc_e/Dmax_var_e)
            end where
         elseif (Dmax_varname == 'cDmax') then
            !Dmax stored as cDmax = c*Dmax
            where (conc_e.gt.0)
               dmax_e   = Dmax_var_e/conc_e
            end where
         else
            stop('Dmax_varname should be Nfloes or cDmax')
         end if

         ibrk  = 0
         do i=1,nmesh_e
            if (mom0_e(i).gt.0.) then
               brk_info%conc        = conc_e(i)
               brk_info%thick       = thick_e(i)
               brk_info%mom0        = mom0_e(i)
               brk_info%mom2        = mom2_e(i)
               brk_info%var_strain  = var_strain_e(i)
               brk_info%dmax        = dmax_e(i)
               brk_info%broken      = broken_e(i)

               call do_breaking(brk_info)
               
               broken_e(i) = brk_info%broken
               !if (broken_e(i).gt..5) then
               !   ibrk  = i
               !   print*,'TWbrk',i,dmax_e(i),brk_info%dmax
               !end if
               dmax_e(i)   = brk_info%dmax
            end if
         end do

         Dmax_var_e  = 0.
         if (Dmax_varname == 'Nfloes') then
            !Dmax stored as Nfloes = c/Dmax^2
            where (dmax_e.gt.0.)
               Dmax_var_e  = conc_e/dmax_e**2
            end where
            !if (ibrk.gt.0) then
            !   print*,'TWbrk2',dmax_e(ibrk),Dmax_var_e(ibrk)
            !end if
         elseif (Dmax_varname == 'cDmax') then
            !Dmax stored as cDmax = c*Dmax
            Dmax_var_e  = conc_e*dmax_e
         else
            stop('Dmax_varname should be Nfloes or cDmax')
         end if
         mesh_e(:,5) = Dmax_var_e
         mesh_e(:,6) = broken_e
         !print*,maxval(broken_e)

      end if
         
      end subroutine wim_step
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine do_wim_integrals(REF_Hs_ICE_)

      ! already used at top
!     use mod_common_wim,only: Hs,Tp,mwd,WAVE_MASK,                     &
!    &      sdf_dir,freq_vec,wavdir

      implicit none
      include "wave_info.h"

      integer,intent(in)      :: REF_Hs_ICE_
      real,dimension(idm,jdm) :: mom0,mom2,mom_dir
      real,dimension(idm,jdm) :: Sfreq,mom_dir0,F
      integer                 :: w,wth
      real                    :: adv_dir,PI,wt_theta,om

      PI = 4*atan(1.0)
      if (ndir.eq.1) then
         wt_theta = 1.
      else
         wt_theta = 2.*PI/real(ndir)
      end if

      !! spectral moments;
      mom0     = 0.0
      mom2     = 0.0
      mom_dir  = 0.0 !! for mwd

      do w=1,n_wave_freq

         om       = 2.*PI*freq_vec(w)
         Sfreq    = 0.  ! freq spec
         mom_dir0 = 0.  ! freq integral for mwd

         do wth=1,n_wavdir
            adv_dir  = -PI/180.0*(wavdir(wth)+90.0)
            Sfreq    = Sfreq     + wt_theta*sdf_dir(:,:,wth,w)
            mom_dir0 = mom_dir0  + wt_theta*sdf_dir(:,:,wth,w)*adv_dir
         end do! wth

         if (REF_Hs_ICE_.eq.1) then
            !! calc for ice
            F  = disp_ratio(:,:,w)
         else
            F  = 1.
         end if
         mom0     = mom0    + abs( wt_om(w)*F*F*Sfreq       )!should be >0
         mom2     = mom2    + abs( wt_om(w)*F*F*Sfreq*om**2 )
         mom_dir  = mom_dir + abs( wt_om(w)*F*F*mom_dir0    )

      end do! w

      !print*,minval(disp_ratio),minval(mom0),minval(mom2)
      !print*,maxval(disp_ratio),maxval(mom0),maxval(mom2)
      Hs          = 0.
      Tp          = 0.
      mwd         = 0.
      WAVE_MASK   = 0.

      where (mom2*mom2>0.)
         WAVE_MASK   = 1.
         Hs          = 4.*sqrt(mom0)
         Tp          = 2.*PI*sqrt(mom0/mom2)
         mwd         = -90. -180./PI*(mom_dir/mom0)!mwd (waves-from in degrees)
            !adv_dir  = -PI/180.0*(wavdir(wth)+90.0)
      end where
      !print*,minval(Hs),minval(Tp),minval(mwd),minval(WAVE_MASK)
      !print*,maxval(Hs),maxval(Tp),maxval(mwd),maxval(WAVE_MASK)
      !stop

      end subroutine do_wim_integrals
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calc_mwd(REF_Hs_ICE_)

      ! already used at top
!     use mod_common_wim,only: mwd                                      &
!    &      sdf_dir,freq_vec,wavdir

      implicit none
      include "wave_info.h"

      integer,intent(in)      :: REF_Hs_ICE_
      real,dimension(idm,jdm) :: mom0,mom_dir
      real,dimension(idm,jdm) :: Sfreq,mom_dir0,F
      integer                 :: w,wth
      real                    :: adv_dir,PI,wt_theta,om

      PI = 4*atan(1.0)
      if (ndir.eq.1) then
         wt_theta = 1.
      else
         wt_theta = 2.*PI/real(ndir)
      end if

      !! spectral moments;
      mom0     = 0.0
      mom_dir  = 0.0 !! for mwd

      do w=1,n_wave_freq

         om       = 2.*PI*freq_vec(w)
         Sfreq    = 0.  ! freq spec
         mom_dir0 = 0.  ! freq integral for mwd

         do wth=1,n_wavdir
            adv_dir  = -PI/180.0*(wavdir(wth)+90.0)
            Sfreq    = Sfreq     + wt_theta*sdf_dir(:,:,wth,w)
            mom_dir0 = mom_dir0  + wt_theta*sdf_dir(:,:,wth,w)*adv_dir
         end do! wth

         if (REF_Hs_ICE_.eq.1) then
            !! calc for ice
            F  = disp_ratio(:,:,w)
         else
            F  = 1.
         end if
         mom0     = mom0    + abs( wt_om(w)*F*F*Sfreq    )
         mom_dir  = mom_dir + abs( wt_om(w)*F*F*mom_dir0 )

      end do! w

      mwd   = 0.!problematic?
      where (mom0>0.)
         mwd         = -90. -180./PI*(mom_dir/mom0)!mwd (waves-from in degrees)
            !adv_dir  = -PI/180.0*(wavdir(wth)+90.0)
      end where

      end subroutine calc_mwd
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine check_grid()

      use mod_gridinfo
      use mod_file_utils

      implicit none
      character(len=80) :: afile,bfile
      integer,parameter :: Nrecs = 7

      afile = trim(outdir)//'/binaries/wim_grid.a'
      bfile = trim(outdir)//'/binaries/wim_grid.b'

      !! write variable names to .b file (description of stuff in .a file)
 111  format(a)
 112  format(i2.2,a)  
      call wrt_bfile_header(bfile,Nrecs)
      open(unit=2,file=trim(bfile),position = 'append')

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(2,111)    "Record number and name:"
      write(2,112)    1,    "       X"
      write(2,112)    2,    "       Y"
      write(2,112)    3,    "       scuy"
      write(2,112)    4,    "       scvx"
      write(2,112)    5,    "       scp2"
      write(2,112)    6,    "       scp2i"
      write(2,112)    7,    "       LANDMASK"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      close(2)

      !! dump arrays to .a file
      !! - X,Y,scuy,scvx,scp2,scp2i,LANDMASK
      call wrt_afile_2d(afile ,X          ,ii,jj  ,'replace',1)
      call wrt_afile_2d(afile ,Y          ,ii,jj  ,'old'    ,2)
      call wrt_afile_2d(afile ,scuy       ,ii,jj  ,'old'    ,3)
      call wrt_afile_2d(afile ,scvx       ,ii,jj  ,'old'    ,4)
      call wrt_afile_2d(afile ,scp2       ,ii,jj  ,'old'    ,5)
      call wrt_afile_2d(afile ,scp2i      ,ii,jj  ,'old'    ,6)
      call wrt_afile_2d(afile ,LANDMASK   ,ii,jj  ,'old'    ,7)

      end subroutine check_grid
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine check_init()

      use mod_gridinfo
      use mod_common_wim
      use mod_file_utils

      implicit none
      character(len=80) :: afile,bfile
      integer,parameter :: Nrecs = 6
      character(len=16) :: cts

      cts   = model_year_info%cdatetime
      afile = trim(outdir)//'/binaries/init/wim_init'
     &         //cts//'.a'
      bfile = trim(outdir)//'/binaries/init/wim_init'
     &         //cts//'.b'

      !! write variable names to .b file (description of stuff in .a file)
 111  format(a)
 112  format(i2.2,a)
      call wrt_bfile_header(bfile,Nrecs,.true.)
      open(unit=2,file=trim(bfile),position = 'append')

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(2,111)    "Record number and name:"
      write(2,112)    1,    "       icec"
      write(2,112)    2,    "       iceh"
      write(2,112)    3,    "       Dmax"
      write(2,112)    4,    "       Hs"
      write(2,112)    5,    "       Tp"
      write(2,112)    6,    "       mwd"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      close(2)

      !! dump arrays to .a file
      !! - X,Y,scuy,scvx,scp2,scp2i,LANDMASK
      call wrt_afile_2d(afile ,icec    ,ii,jj  ,'replace',1)
      call wrt_afile_2d(afile ,iceh    ,ii,jj  ,'old'    ,2)
      call wrt_afile_2d(afile ,dfloe   ,ii,jj  ,'old'    ,3)
      call wrt_afile_2d(afile ,Hs      ,ii,jj  ,'old'    ,4)
      call wrt_afile_2d(afile ,Tp      ,ii,jj  ,'old'    ,5)
      call wrt_afile_2d(afile ,mwd     ,ii,jj  ,'old'    ,6)

      end subroutine check_init
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine check_prog(OPT)

      use mod_gridinfo
      use mod_common_wim
      use mod_file_utils

      implicit none
      integer,intent(in)   :: OPT
      character(len=80)    :: afile,bfile
      character(len=16)    :: cts
      integer,parameter    :: Nrecs = 6

      if (OPT.eq.0) then
         !!dump only for progress report
         cts   = model_year_info%cdatetime
         afile = trim(outdir)//'/binaries/prog/wim_prog'                &
     &               //cts//'.a'
         bfile = trim(outdir)//'/binaries/prog/wim_prog'                &
     &               //cts//'.b'
      else
#if 1
         !!final output - use start time as output time in filename
         cts   = model_year_info_start%cdatetime
#else
         !!final output - use end time as output time in filename
         cts   = model_year_info%cdatetime
#endif
         afile = trim(outdir)//'/binaries/final/wim_out'
     &               //cts//'.a'
         bfile = trim(outdir)//'/binaries/final/wim_out'
     &               //cts//'.b'
      end if

      !! write variable names to .b file (description of stuff in .a file)
 111  format(a)
 112  format(i2.2,a)  
      call wrt_bfile_header(bfile,Nrecs,.true.)
      open(unit=2,file=trim(bfile),position = 'append')

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(2,111)    "Record number and name:"
      write(2,112)    1,    "       Dmax"
      write(2,112)    2,    "       tau_x"
      write(2,112)    3,    "       tau_y"
      write(2,112)    4,    "       Hs"
      write(2,112)    5,    "       Tp"
      write(2,112)    6,    "       mwd"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      close(2)

      !! dump arrays to .a file
      !! - X,Y,scuy,scvx,scp2,scp2i,LANDMASK
      call wrt_afile_2d(afile ,dfloe   ,ii,jj  ,'replace',1)
      call wrt_afile_2d(afile ,tau_x   ,ii,jj  ,'old'    ,2)
      call wrt_afile_2d(afile ,tau_y   ,ii,jj  ,'old'    ,3)
      call wrt_afile_2d(afile ,Hs      ,ii,jj  ,'old'    ,4)
      call wrt_afile_2d(afile ,Tp      ,ii,jj  ,'old'    ,5)
      call wrt_afile_2d(afile ,mwd     ,ii,jj  ,'old'    ,6)

      end subroutine check_prog
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine wrt_bfile_header(bfile,Nrecs,add_time)
      character(len=*),intent(in)   :: bfile
      integer,intent(in)            :: Nrecs
      logical,optional,intent(in)   :: add_time

 111  format(a)
 112  format(i2.2,a)  
 113  format(i3.3,a)  
 114  format(i4.4,a)  
 211  format(f9.1,a)  

      !! write .b file (description of stuff in .a file)
      open(unit=2,file=trim(bfile),status = 'replace')

      write(2,112) Nrecs,  "         Nrecs    # Number of records"
      write(2,112) 1,      "         Norder   # Storage order "//       &
     &               "[column-major (F/matlab) = 1, row-major (C) = 0]"
      write(2,114) ii,     "       nx       "//                         &
     &               "# Record length in x direction (elements)"
      write(2,114) jj,     "       ny       "//                         &
     &               "# Record length in y direction (elements)"
      write(2,111) ""

      !! if cts is given add a line with the model time at time of dump
      if (present(add_time)) then
         if (add_time) then
            write(2,'(a)') , model_year_info_start%cdatetime//
     &          "  t_start    # Model time of WIM call"
            write(2,'(a)') , model_year_info%cdatetime//
     &          "  t_out      # Model time of output"
         end if
         write(2,111) ""
      end if

      close(2)

      end subroutine wrt_bfile_header
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine floe_scaling(Dave,dmax,moment)

      implicit none

      real,    intent(out) :: Dave
      real,    intent(in)  :: dmax
      integer, intent(in)  :: moment
      integer              :: m,mm 
      real                 :: nm,nm1,dm,nsum,ndsum,r,ffac

      Dave  = max(dmax**moment,wim_info%dmin**moment)
      ffac  = wim_info%fragility*wim_info%xi**2

      !! determine no of breaks
      r     = dmax/wim_info%dmin
      mm    = 0  
      do while ( r .ge. wim_info%xi ) 
         !! if r<xi,no more breaking
         !! - don't change Dave
         r  = r/wim_info%xi 
         mm = mm+1 
      end do

      if ( mm .gt. 0 ) then 
         nm1   = 1.
         dm    = dmax   !floe length
         nsum  = 0      !eventually \sum_m=0^mm.n_m
         ndsum = 0      !eventually \sum_m=0^mm.n_m.d_m^moment

         do m = 0,mm-1
            nm    = (1.-wim_info%fragility)*nm1     !no of floes at dm
            nsum  = nsum +nm
            ndsum = ndsum+nm*dm**moment
            !!
            nm1   = nm1*ffac
            dm    = dm/wim_info%xi
         end do

         !m=M
         nsum  = nsum+nm1
         ndsum = ndsum+nm1*dm**moment
         

         Dave  = ndsum/nsum
      end if

!     !test subroutine (NB this will crash after 1st call, since status='new')
!     print*,'in floe_scaling'
!     open(unit=999,file='fsd_test_old.txt',                            &
!    &     action='write',status='new')
!     write(999,*) Dave
!     write(999,*) Dmax
!     write(999,*) fragility
!     write(999,*) xi
!     write(999,*) Dmin
!     close(999)

      end subroutine floe_scaling
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine floe_scaling_smooth(Dave,dmax,mom)
      implicit none

      real,    intent(out) :: Dave
      real,    intent(in)  :: dmax
      integer, intent(in)  :: mom
      integer              :: n
      real                 :: b,Dthresh,fsd_exp,A

      n        = mom
      fsd_exp  = 2+log(wim_info%fragility)/log(wim_info%xi)!!power law exponent: P(d>D)=(D_min/D)^fsd_exp;
      b        = n-fsd_exp
      Dthresh  = 200.

      !! calculate <D^n> from Dmax
      !! - uniform dist for larger floes
      Dave  = Dmax**n

      !! small floes
      if (Dmax.le.wim_info%Dmin) then
         Dave  = wim_info%Dmin**n
      elseif (Dmax.le.Dthresh) then
         !! bigger floes
         A     = (fsd_exp*exp(fsd_exp*(log(wim_info%Dmin)+log(Dmax))))
         A     = A/( exp(fsd_exp*log(Dmax))
     &               -exp(fsd_exp*log(wim_info%Dmin)) )
         Dave  = -(A/b)*(exp(b*log(wim_info%Dmin))-exp(b*log(Dmax)))
      end if

!!     !test subroutine (NB this will crash after 1st call, since status='new')
!      print*,'in floe_scaling_smooth'
!      open(unit=999,file='fsd_test_smooth.txt',                         &
!     &     action='write',status='new')
!      write(999,*) Dave
!      write(999,*) Dmax
!      write(999,*) mom
!      write(999,*) fragility
!      write(999,*) xi
!      write(999,*) Dmin
!      close(999)

      end subroutine floe_scaling_smooth


      subroutine do_breaking(brk_info)

      implicit none

      type(break_info),intent(inout)   :: brk_info
      logical :: BREAK_CRIT
      real    :: om,guess,inputs(2),lam,kice,kwtr,int_adm,              &
     &            Tp,modT,argR,argT,pi

      pi = 4*atan(1.)
      BREAK_CRIT   = (brk_info%conc.gt.0.)                            !ice present
     &               .and.(2*brk_info%var_strain.gt.wim_info%epsc**2) !big enough waves

      if (BREAK_CRIT) then
         ! =============================================================
         !calculate wavelength
         !inputs
         om    = sqrt(brk_info%mom2/brk_info%mom0) !omega=2*pi/T
         guess = om**2/wim_info%gravity
         
         inputs(1)   = brk_info%thick !absolute thickness
         inputs(2)   = om

         ! rest of arguments are outputs, but only need kice
         ! so use Tp_e,lam_e as dummy variables
         call RTparam_outer(lam,kice,kwtr,int_adm,                      &
     &            Tp,modT,argR,argT,                                    &
     &            inputs,guess)
         ! =============================================================

         lam = 2*pi/kice
         if (lam.lt.2*brk_info%dmax) then
            !breaking occurs
            brk_info%dmax     = max(wim_info%Dmin,lam/2.) !keep above Dmin
            brk_info%broken   = 1.                        !this records ice has been broken this time
            !print*,i,broken_e(i)
         end if
      end if

      end subroutine do_breaking
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module mod_wavesice
