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

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!variables

      !!TODO: make these inputs
      !integer,parameter :: STEADY      = 0
      !integer,parameter :: ADV_OPT     = 0
         !!0, not periodic; 2, periodic in y only; 1, periodic in both x,y;

      integer,parameter :: REF_Hs_ICE  = 0
      integer,parameter :: DO_ATTEN    = 1
      integer,parameter :: USE_ICE_VEL = 0
      !!
      real,parameter    :: CFL      = 0.7
      !!
      integer,parameter :: nw    = n_wave_freq
      integer,parameter :: ndir  = n_wavdir

      !!common to all subroutines
      real,dimension(idm,jdm)    :: WAVE_MASK2
      real,dimension(idm,jdm,nw) ::                                     &
     &   ag_eff,ap_eff,wlng_ice,atten_nond,damping,disp_ratio
      real,dimension(nw)   :: wt_om

      !!needed by outer routine
      real     :: amin,amax,dt
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      contains

      subroutine wim_init(STEADY)
      !! set initial waves
      !! set attenuation coefficients
      !! - need to do this each time ice is updated

      use mod_gridinfo     !! grid info
      use mod_wim_prams,   only: young,visc_rp,gravity   !! ice parameters
      use mod_common_wim,  only: ICE_MASK,iceh           !! initial ice fields

      implicit none

      integer, intent(in) :: STEADY
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

      subroutine wim_step(SCATMOD,ADV_DIM,ADV_OPT                       &
     &                     ,DO_BREAKING,STEADY)
      !! do 1 time step

      use mod_gridinfo     !! grid info
      use mod_wim_prams    !! ice parameters
!     use mod_wim_init     !! initial wave fields
      use mod_common_wim   !! initial ice fields
      use mod_adv_atten    !! advection routine

      implicit none

      integer,intent(in)   :: SCATMOD,ADV_DIM,ADV_OPT                   &
     &                        ,DO_BREAKING,STEADY

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
      real     :: sig_strain,Pstrain,P_crit,wlng_crest,Dc
      real     :: om,ommin,ommax,om1,lam1,lam2,dom
      real     :: F,Dave,c1d,kice,tmp
      integer  :: jcrest
      logical  :: BREAK_CRIT

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
               if (dfloe(i,j).lt.200.0) then
                  call floe_scaling(Dave,dfloe(i,j))
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

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !! Transport waves and attenuate them
         if (SCATMOD.eq.0) then
            !! Simple attenuation scheme - doesn't conserve scattered energy
            call adv_atten_simple(sdf_dir(:,:,:,w)                      &
     &         ,S_freq,taux_om,tauy_om                                  &
     &         ,atten_dim,damp_dim,ag_eff(:,:,w),dt,ADV_DIM,ADV_OPT)
            !print*,'Spectrum at',itest,jtest
            !print*,S_freq(itest,jtest)
            !print*,sdf_dir(itest,jtest,:,w)
         elseif (SCATMOD.eq.1) then
            !! same as SCATMOD==0, but scattered energy
            !! is distributed isotropically
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

      end do!! end spectral loop;

      if (REF_Hs_ICE.eq.1) then
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
      call calc_mwd(REF_Hs_ICE)

      if (.false.) then
         !!test: print atten coeff
         print*,'********************************************'
         print*,'Max atten coeff (per metre):'
         print*,maxval(atten_dim)
         print*,'********************************************'
      end if

      if ((STEADY.eq.0).and.(DO_BREAKING.eq.0)) then
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
            P_crit   = (1-DO_BREAKING)+exp(-1.0)
               !!this is critical prob if monochromatic wave
   
            !! FLOE BREAKING:
            BREAK_CRIT  = ( Pstrain.ge.P_crit ) !! breaks if larger than this

            !! Leave FSD alone if DO_BREAKING==0
            BREAK_CRIT  = BREAK_CRIT.and.(DO_BREAKING.eq.1)
   
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
      !use mod_common_wim
      use mod_file_utils

      implicit none
      character(len=80),intent(in)   :: outdir
      character(len=80)    :: afile,bfile
      character(len=10)    :: stat1,stat2

      afile = trim(outdir)//'/binaries/wim_grid.a'
      bfile = trim(outdir)//'/binaries/wim_grid.b'
      stat1 = 'replace'
      stat2 = 'old'

 111  format(a)
 112  format(i2.2,a)  
 113  format(i3.3,a)  
 114  format(i4.4,a)  

      !! write .b file (description of stuff in .a file)
      open(unit=2,file=trim(bfile),status = 'replace')

      write(2,112)                                                      &
     &  7,          "       Number of records"
      write(2,112)                                                      &
     &  1,          "       Storage order "//                           &
     &              "[column-major (F/matlab) = 1, row-major (C) = 0]"
      write(2,114)                                                      &
     &  ii,          "      Record length in x direction (elements)"
      write(2,114)                                                      &
     &  jj,          "      Record length in y direction (elements)"
      write(2,111) " "

      write(2,111)    "Record number and name:"
      write(2,112)                                                      &
     &  1,    "       X"
      write(2,112)                                                      &
     &  2,    "       Y"
      write(2,112)                                                      &
     &  3,    "       scuy"
      write(2,112)                                                      &
     &  4,    "       scvx"
      write(2,112)                                                      &
     &  5,    "       scp2"
      write(2,112)                                                      &
     &  6,    "       scp2i"
      write(2,112)                                                      &
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
      subroutine check_init(outdir)

      use mod_gridinfo
      use mod_common_wim
      use mod_file_utils

      implicit none
      character(len=80),intent(in)   :: outdir
      character(len=80)    :: afile,bfile
      character(len=10)    :: stat1,stat2

      afile = trim(outdir)//'/binaries/wim_init.a'
      bfile = trim(outdir)//'/binaries/wim_init.b'
      stat1 = 'replace'
      stat2 = 'old'

 111  format(a)
 112  format(i2.2,a)  
 113  format(i3.3,a)  
 114  format(i4.4,a)  

      !! write .b file (description of stuff in .a file)
      open(unit=2,file=trim(bfile),status = 'replace')

      write(2,112)                                                      &
     &  6,          "       Number of records"
      write(2,112)                                                      &
     &  1,          "       Storage order "//                           &
     &              "[column-major (F/matlab) = 1, row-major (C) = 0]"
      write(2,114)                                                      &
     &  ii,          "      Record length in x direction (elements)"
      write(2,114)                                                      &
     &  jj,          "      Record length in y direction (elements)"
      write(2,111) ""

      write(2,111)    "Record number and name:"
      write(2,112)                                                      &
     &  1,    "       icec"
      write(2,112)                                                      &
     &  2,    "       iceh"
      write(2,112)                                                      &
     &  3,    "       Dmax"
      write(2,112)                                                      &
     &  4,    "       Hs"
      write(2,112)                                                      &
     &  5,    "       Tp"
      write(2,112)                                                      &
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
      subroutine check_prog(SCATMOD,nn,OPT,outdir)

      use mod_gridinfo
      use mod_common_wim
      use mod_file_utils

      implicit none
      integer,intent(in)            :: SCATMOD,nn(2),OPT
      character(len=80),intent(in)  :: outdir
      character(len=80)             :: afile,bfile
      character(len=10)             :: stat1,stat2
      character(len=4)              :: cts
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
      stat1 = 'replace'
      stat2 = 'old'

 111  format(a)
 112  format(i2.2,a)  
 113  format(i3.3,a)  
 114  format(i4.4,a)  
 211  format(f8.1,a)  

      !! write .b file (description of stuff in .a file)
      open(unit=2,file=trim(bfile),status = 'replace')

      write(2,112)                                                      &
     &  5,           "         Number of records"
      write(2,112)                                                      &
     &  1,          "         Storage order "//                         &
     &              "[column-major (F/matlab) = 1, row-major (C) = 0]"
      write(2,114)                                                      &
     &  ii,           "        Record length in x direction (elements)"
      write(2,114)                                                      &
     &  jj,           "        Record length in y direction (elements)"
      write(2,211)                                                      &
     &  t_out,           "   Model time of output (s)"
      write(2,111) ""

      write(2,111)    "Record number and name:"
      write(2,112)                                                      &
     &  1,    "       Dmax"
      write(2,112)                                                      &
     &  2,    "       tau_x"
      write(2,112)                                                      &
     &  3,    "       tau_y"
      write(2,112)                                                      &
     &  4,    "       Hs"
      write(2,112)                                                      &
     &  5,    "       Tp"
!     write(2,112)                                                      &
!    &  6,    "       mwd"

      close(2)

      !! dump arrays to .a file
      !! - X,Y,scuy,scvx,scp2,scp2i,LANDMASK
      call wrt_afile_2d(afile ,dfloe   ,ii,jj  ,stat1,1)
      call wrt_afile_2d(afile ,tau_x   ,ii,jj  ,stat2,2)
      call wrt_afile_2d(afile ,tau_y   ,ii,jj  ,stat2,3)
      call wrt_afile_2d(afile ,Hs      ,ii,jj  ,stat2,4)
      call wrt_afile_2d(afile ,Tp      ,ii,jj  ,stat2,5)
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
      !print*,'dave (in f-s)',Dave

      end subroutine floe_scaling
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module mod_wavesice
