      !mod_wavesice.F
      !Author: Timothy Williams
      !Date: 20141119, 05:22:32 CET
      module mod_wavesice

      use mod_gridinfo     !! grid info
!     use mod_wim_prams    !! ice parameters
!     use mod_wim_init     !! initial wave fields
      use mod_common_wim   !! initial ice fields
!     use mod_adv_atten    !! advection routine

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
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      public   :: wim_init,wim_step,do_wim_integrals
     &            ,check_grid,check_init,check_prog
     &            ,amin,amax,dt,WAVE_MASK2
     &            ,USE_ICE_VEL,REF_HS_ICE
      !,calc_mwd

      contains

      subroutine wim_init(STEADY,DO_ATTEN,CFL)
      !! set initial waves
      !! set attenuation coefficients
      !! - need to do this each time ice is updated

      use mod_gridinfo     !! grid info
      use mod_wim_prams,   only: young,visc_rp,gravity   !! ice parameters
      use mod_common_wim,  only: ICE_MASK,iceh           !! initial ice fields

      implicit none

      integer, intent(in) :: STEADY,DO_ATTEN
      real,    intent(in) :: CFL
      !!
      real     :: PI,tmp
      real     :: kwtr,kice,modT,int_adm,argT,argR
      real,dimension(nw)   :: wt_simp,wlng,ag,ap
      integer  :: i,j,w
      real     :: adv_dir,dtheta,om,dom,guess,inputs(4)
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
      wlng     = gravity/(2*PI*freq_vec**2)
      ap       = sqrt(gravity*wlng/(2*PI))   ! Phase speed
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
               guess = om**2/gravity!guess for root
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
            inputs(3)   = young
            inputs(4)   = visc_rp

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

      subroutine wim_step(wim_step_prams,nmesh_e,mesh_e)
      !! do 1 time step

      use mod_gridinfo     !! grid info
      use mod_wim_prams    !! ice parameters
!     use mod_wim_init     !! initial wave fields
      use mod_common_wim   !! initial ice fields
      use mod_adv_atten    !! advection routine
      use mod_interp, only: pwl_interp_2d_real

      implicit none

      integer,parameter :: n_wim_step_prams  = 7
      integer,dimension(n_wim_step_prams),intent(in) :: wim_step_prams  &
         !SCATMOD,ADV_DIM,ADV_OPT,BRK_OPT,STEADY,DO_ATTEN,FSD_OPT
      integer,intent(in)   :: nmesh_e

      real,dimension(nmesh_e,nmesh_vars),intent(inout),optional         &
     &   :: mesh_e
         ! this is [x_e,y_e,conc_e,thick_e,Nfloes_e,broken_e]
         ! -1st 4 are unchanged; last 2 get updated

      integer  :: SCATMOD,ADV_DIM,ADV_OPT                               &
     &            ,BRK_OPT,STEADY,DO_ATTEN,FSD_OPT

      real,dimension(idm,jdm)    :: mom0,mom2,var_strain                &
     &   ,S_freq,atten_dim,damp_dim                                     &
     &   ,taux_om,tauy_om,tmp1                                          &
     &   ,mom0w,mom2w
      integer  :: i,j,w

      !!to top up waves
      integer  :: wth
      real  :: adv_dir,PI

!     !!for testing
!     real     :: wt,Hs_test1,Hs_test2
      real  :: E_tot

      !!for breaking
      real     :: sig_strain,Pstrain,P_crit0,P_crit,wlng_crest,Dc
      real     :: om,ommin,ommax,om1,lam1,lam2,dom
      real     :: F,Dave,c1d,kice,tmp
      integer  :: jcrest
      logical  :: BREAK_CRIT

      logical              :: INTERP_MESH,crit_e
      real                 :: wim_break_damage,Tp_e,lam_e
      real,dimension(nmesh_e)  :: x_e,y_e                    !centres of elements
      real,dimension(nmesh_e)  :: mom0_e,mom2_e,var_strain_e !interpolated values onto elements
      real,dimension(nmesh_e)  :: conc_e,thick_e,Nfloes_e
     &                              ,broken_e,dmax_e
      real     :: kwtr,modT,int_adm,argT,argR,guess,inputs(4)
      !
      INTERP_MESH = .false.
      if (present(mesh_e)) then
         ! check if there is any ice
         INTERP_MESH = maxval(mesh_e).gt.0.
      end if

      !wim_step_prams   = (/
      !SCATMOD,ADV_DIM,ADV_OPT,BRK_OPT,STEADY,DO_ATTEN,FSD_OPT
      !/)
      SCATMOD  = wim_step_prams(1)
      ADV_DIM  = wim_step_prams(2)
      ADV_OPT  = wim_step_prams(3)
      BRK_OPT  = wim_step_prams(4)
      STEADY   = wim_step_prams(5)
      DO_ATTEN = wim_step_prams(6)
      FSD_OPT  = wim_step_prams(7)

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

      do w=1,nw

         !! CALC DIMENSIONAL ATTEN COEFF;
         atten_dim   = 0.0
         damp_dim    = 0.0
         do j = 1,jj
         do i = 1,ii

            if ((ICE_MASK(i,j).eq.1.0).and.(DO_ATTEN.eq.1)) then

               !! get expected no of floes met per unit
               !!  distance if travelling in a line;
               if (dfloe(i,j).le.200.0) then
                  if (FSD_OPT.eq.0) then
                     !! renormalisation group method (Dumont et al, 2011)
                     call floe_scaling(Dave,dfloe(i,j))
                  elseif (FSD_OPT.eq.1) then
                     !! smoothed power law method (Montiel et al, 2015)
                     call floe_scaling_smooth(Dave,dfloe(i,j),1)
                  end if
               else
                  Dave  = dfloe(i,j)
                  !! uniform lengths
               end if
               c1d = icec(i,j)/Dave!! floes per unit length;

               !! ENERGY attenuation coeff;
               atten_dim(i,j) = atten_nond(i,j,w)*c1d      !!scattering
               damp_dim(i,j)  = 2*damping(i,j,w)*icec(i,j) !!damping

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
         if (.true.) then!!new way (split advection from attenuation)

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
            tmp1  = rhowtr*gravity*taux_om/ap_eff(:,:,w)  !![Pa*s]
            tau_x = tau_x+wt_om(w)*tmp1                   !![Pa]
            tmp1  = rhowtr*gravity*tauy_om/ap_eff(:,:,w)  !![Pa*s]
            tau_y = tau_y+wt_om(w)*tmp1                   !![Pa]

         end if!!DO_ATTEN
         else!!old way

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !! Transport waves and attenuate them
         if (SCATMOD.eq.0) then
            !! Simple attenuation scheme - doesn't conserve scattered energy
            print*,'Simple atten 0'
            call adv_atten_simple(sdf_dir(:,:,:,w)                      &
     &         ,S_freq,taux_om,tauy_om                                  &
     &         ,atten_dim,damp_dim,ag_eff(:,:,w),dt,ADV_DIM,ADV_OPT)
            !print*,'Spectrum at',itest,jtest
            !print*,S_freq(itest,jtest)
            !print*,sdf_dir(itest,jtest,:,w)
         elseif (SCATMOD.eq.1) then
            !! same as SCATMOD==0, but scattered energy
            !! is distributed isotropically
            print*,'Isotropic atten 0'
            call adv_atten_isotropic(sdf_dir(:,:,:,w)                   &
     &         ,S_freq,taux_om,tauy_om                                  &
     &         ,atten_dim,damp_dim,ag_eff(:,:,w),dt,ADV_DIM,ADV_OPT)
         end if
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         !! integrate stress densities over frequency
         !! TODO: check if this is correct for ice-covered water
         tmp1  = rhowtr*gravity*taux_om/ap_eff(:,:,w)   !![Pa*s]
         tau_x = tau_x+wt_om(w)*tmp1                     !![Pa]
         tmp1  = rhowtr*gravity*tauy_om/ap_eff(:,:,w)  !![Pa*s]
         tau_y = tau_y+wt_om(w)*tmp1                     !![Pa]

         end if!!old way (end)

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

      if (.false.) then
         !!test: print atten coeff
         print*,'********************************************'
         print*,'Max atten coeff (per metre):'
         print*,maxval(atten_dim)
         print*,'********************************************'
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
      P_crit0  = 0.
      if (BRK_OPT.eq.0) then
         P_crit0  = 1.!!added to P_crit: Pstrain<=1 so now ice will never break
      end if

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
            P_crit   = P_crit0+exp(-1.0)
               !!this is critical prob if monochromatic wave
   
            !! FLOE BREAKING:
            BREAK_CRIT  = ( Pstrain.ge.P_crit ) !! breaks if larger than this

            !! Leave FSD alone if BRK_OPT==0
            BREAK_CRIT  = BREAK_CRIT.and.(BRK_OPT.gt.0)
   
            if (BREAK_CRIT) then
               !! use crest period to work out wavelength
               !! - interpolate using pre-calc'd wavelengths
               !! - half this is max poss floe length;
               om    = 2*PI/Tp(i,j)
               ommin = 2*PI*freq_vec(1)
               ommax = 2*PI*freq_vec(nw)
               if (om.le.ommin) then
                  wlng_crest  = wlng_ice(i,j,1)
               elseif (om.ge.ommax) then
                  wlng_crest  = wlng_ice(i,j,nw)
               else
                  jcrest      = floor((om-ommin+dom)/dom)
                  !print*,'crest lam',jcrest,nw
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

      if (INTERP_MESH) then
         x_e      = mesh_e(:,1)
         y_e      = mesh_e(:,2)
         conc_e   = mesh_e(:,3)
         thick_e  = mesh_e(:,4)
         Nfloes_e = mesh_e(:,5)
         broken_e = mesh_e(:,6)

         call pwl_interp_2d_real(ii,jj,X(:,1),Y(1,:),mom0
     &      ,nmesh_e,x_e,y_e,mom0_e)
         call pwl_interp_2d_real(ii,jj,X(:,1),Y(1,:),mom2
     &         ,nmesh_e,x_e,y_e,mom2_e)
         call pwl_interp_2d_real(ii,jj,X(:,1),Y(1,:),var_strain
     &      ,nmesh_e,x_e,y_e,var_strain_e)

         dmax_e   = 0.
         where (Nfloes_e.gt.0)
            dmax_e   = sqrt(conc_e/Nfloes_e)
         end where

         do i=1,nmesh_e
            crit_e   = (conc_e(i).gt.0.)                    !ice present
     &                 .and.(2*var_strain_e(i).gt.epsc**2)  !big enough waves

            if (crit_e) then
               ! =============================================================
               !calculate wavelength
               !inputs
               om    = sqrt(mom2_e(i)/mom0_e(i)) !omega=2*pi/T
               guess = om**2/gravity
               
               inputs(1)   = thick_e(i) !absolute thickness
               inputs(2)   = om
               inputs(3)   = young
               inputs(4)   = visc_rp

               ! rest of arguments are outputs, but only need kice
               ! so use Tp_e,lam_e as dummy variables
               call RTparam_outer(lam_e,kice,kwtr,int_adm,              &
     &                  Tp_e,modT,argR,argT,                            &
     &                  inputs,guess)
               ! =============================================================

               lam_e = 2*pi/kice
               if (lam_e.lt.2*dmax_e(i)) then
                  !breaking occurs
                  dmax_e(i)   = max(Dmin,lam_e/2.) !keep above Dmin
                  broken_e(i) = 1.                 !this records ice has been broken this time
                  !print*,i,broken_e(i)
               end if
            end if
         end do

         where (dmax_e.gt.0.)
            Nfloes_e   = conc_e/dmax_e**2
         end where
         mesh_e(:,5) = Nfloes_e
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
      subroutine check_grid(outdir)

      use mod_gridinfo
      use mod_file_utils

      implicit none
      character(len=*),intent(in)   :: outdir
      character(len=80)             :: afile,bfile
      integer,parameter             :: Nrecs = 7

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
      subroutine check_init(outdir)

      use mod_gridinfo
      use mod_common_wim
      use mod_file_utils

      implicit none
      character(len=*),intent(in)   :: outdir
      character(len=80)             :: afile,bfile
      integer,parameter             :: Nrecs = 6

      afile = trim(outdir)//'/binaries/wim_init.a'
      bfile = trim(outdir)//'/binaries/wim_init.b'

      !! write variable names to .b file (description of stuff in .a file)
 111  format(a)
 112  format(i2.2,a)
      call wrt_bfile_header(bfile,Nrecs)
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
      subroutine check_prog(nn,OPT,outdir)

      use mod_gridinfo
      use mod_common_wim
      use mod_file_utils

      implicit none
      integer,intent(in)            :: nn(2),OPT
      character(len=80),intent(in)  :: outdir
      character(len=80)             :: afile,bfile
      character(len=4)              :: cts
      integer,parameter             :: Nrecs = 5
      !!
      integer           :: n_step,n_total
      character(len=6)  :: frmt
      real              :: t_out

      n_step   = nn(1)
      n_total  = nn(2)
      if (n_total.le.10) then
         frmt  = '(i1.1)'
      elseif (n_total.le.100) then
         frmt  = '(i2.2)'
      elseif (n_total.le.1000) then
         frmt  = '(i3.3)'
      else
         frmt  = '(i4.4)'
      end if

      if (OPT.eq.0) then
         !!dump only for progress report
         write(cts,frmt) n_step
         afile = trim(outdir)//'/binaries/prog/wim_prog'                &
     &               //trim(cts)//'.a'
         bfile = trim(outdir)//'/binaries/prog/wim_prog'                &
     &               //trim(cts)//'.b'
         t_out = n_step*dt
      else
         !!final output
         afile = trim(outdir)//'/binaries/wim_out.a'
         bfile = trim(outdir)//'/binaries/wim_out.b'
         t_out = n_total*dt
      end if

      !! write variable names to .b file (description of stuff in .a file)
 111  format(a)
 112  format(i2.2,a)  
      call wrt_bfile_header(bfile,Nrecs,t_out)
      open(unit=2,file=trim(bfile),position = 'append')

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(2,111)    "Record number and name:"
      write(2,112)    1,    "       Dmax"
      write(2,112)    2,    "       tau_x"
      write(2,112)    3,    "       tau_y"
      write(2,112)    4,    "       Hs"
      write(2,112)    5,    "       Tp"
!     write(2,112)    6,    "       mwd"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      close(2)

      !! dump arrays to .a file
      !! - X,Y,scuy,scvx,scp2,scp2i,LANDMASK
      call wrt_afile_2d(afile ,dfloe   ,ii,jj  ,'replace',1)
      call wrt_afile_2d(afile ,tau_x   ,ii,jj  ,'old'    ,2)
      call wrt_afile_2d(afile ,tau_y   ,ii,jj  ,'old'    ,3)
      call wrt_afile_2d(afile ,Hs      ,ii,jj  ,'old'    ,4)
      call wrt_afile_2d(afile ,Tp      ,ii,jj  ,'old'    ,5)
!     call wrt_afile_2d(afile ,mwd     ,ii,jj  ,'old'    ,6)

      end subroutine check_prog
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine wrt_bfile_header(bfile,Nrecs,t_out)
      character(len=*),intent(in)   :: bfile
      integer,intent(in)            :: Nrecs
      real,optional,intent(in)      :: t_out

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

      !! if t_out is given add a line with the model time at time of dump
      if (present(t_out)) then
         write(2,211) t_out,   "  t_out    "//                          &
     &               "# Model time of output (s)"
      end if

      write(2,111) ""
      close(2)

      end subroutine wrt_bfile_header
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

      !print*,' '
      !print*,'dmax,dmin',dmax,dmin
      !print*,'r,mm',r,mm

      do while ( r .gt. xi ) 
         r  = r/xi 
         mm = mm+1 
         !print*,'r,mm',r,mm
         !print*,' '
      end do

      if ( mm .gt. 0 ) then 
         n     = 0.0
         nsum  = 0.0
         nd    = 0.0
         ndsum = 0.0

         do m = 0,mm
            n     = (1.0-fragility)*(fragility*xi**2)**m
            nd    = n/(xi**m)
            nsum  = nsum +n
            ndsum = ndsum+nd
            dfac  = ndsum/nsum
            !print*,'n,nd,dfac',n,nd,dfac
         end do

         Dave  = dfac*dmax

      else
         Dave = dmin
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
      use mod_wim_prams, only: fragility,xi,Dmin
      implicit none

      real,    intent(out) :: Dave
      real,    intent(in)  :: dmax
      integer, intent(in)  :: mom
      integer              :: n
      real                 :: b,Dthresh,fsd_exp,A

      n        = mom
      fsd_exp  = 2+log(fragility)/log(xi)!!power law exponent: P(d>D)=(D_min/D)^fsd_exp;
      b        = n-fsd_exp
      Dthresh  = 200.

      !! calculate <D^n> from Dmax
      !! - uniform dist for larger floes
      Dave  = Dmax**n

      !! small floes
      if (Dmax.le.Dmin) then
         Dave  = Dmin**n
      elseif (Dmax.le.Dthresh) then
         !! bigger floes
         A     = (fsd_exp*exp(fsd_exp*(log(Dmin)+log(Dmax))))
         A     = A/(exp(fsd_exp*log(Dmax))-exp(fsd_exp*log(Dmin)))
         Dave  = -(A/b)*(exp(b*log(Dmin))-exp(b*log(Dmax)))
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
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module mod_wavesice