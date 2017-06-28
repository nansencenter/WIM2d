      module mod_atten_young

      implicit none
      
      private
      public RTparam_outer

      contains

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!function [damping,kice,kwtr,int_adm,...
      !!            alp_scat,modT,argR,argT] =...
      !!               RT_param_outer(h,om,E,drag_rp,guess)
      subroutine RTparam_outer(outputs,inputs,guess)

      use mod_common_wim

      type(atten_outputs),intent(out)  :: outputs

      real,intent(in)    :: inputs(2),guess
      real               :: h,om
      
      real   :: alp_nd,hnd,zeta_nd
      real   :: Hw_nd,varpi,L
      real   :: damping_nd,drag_rp_nd,vews_nd,coeff_D
      real   :: H_nd = 4.0d0!!infinite depth used in scattering calculation;
      !!
      real   :: ki,kw,coeff_del,BG1,BG2
      !real   :: E,g,rhow,rhoi,nu
      real   :: D,rho

      h        = inputs(1)
      om       = inputs(2)
      
      rho      = wim_info%rhoice/wim_info%rhowtr
      D        = wim_info%young*(h*h*h)/12.0/(1.0-wim_info%poisson**2)
      L        = exp( 0.2*log(D/wim_info%rhowtr/om/om) )
      alp_nd   = om*om/wim_info%gravity*L
      hnd      = h/L
      zeta_nd  = rho*hnd
      
      !!get wavenumber for ice;
      varpi = 1/alp_nd-zeta_nd
      !print*,'DLvH',D,L,varpi,H_nd

      !![ki,BG2,coeff_del]   = gen_root_ice(varpi,H_nd,guess*L);
      call gen_root_ice(ki,BG2,coeff_del,coeff_D,varpi,H_nd,guess*L)
      outputs%kice = ki/L
      !print*,'vHkGB_i',varpi,H_nd,kice,guess*L,BG1
         
      !!get wavenumber for water;
      varpi = 1.0/alp_nd
      Hw_nd = H_nd+zeta_nd
      !![kw,BG1] = gen_root_wtr(varpi,Hw_nd,alp_nd);
      call gen_root_wtr(kw,BG1,varpi,Hw_nd,alp_nd)
      outputs%kwtr  = kw/L
      !print*,'vHkGB_w',varpi,Hw_nd,kwtr,alp_nd,BG1
         
      !!get intrinsic admittance;
      !!|R|^2+int_adm*|T|^2=1
      !print*,'BG1   = ', BG1
      !print*,'BG2   = ', BG2
      outputs%int_adm  = BG1/BG2
         
      !!get dissipative attenuation;
      drag_rp_nd        = wim_info%drag_rp/wim_info%rhowtr/om/L !!RP drag
      vews_nd           = h**3*rho*wim_info%visc_ws/
     &                     (6.*(1.+wim_info%poisson)*om*L**5)
      damping_nd        = -coeff_del*drag_rp_nd -coeff_D*vews_nd
      outputs%damping   = damping_nd/L
      
      !!get interpolated attenuation coefficient;
      call RTparam_fast(outputs%atten_nond
     &                  ,outputs%modT
     &                  ,outputs%argR
     &                  ,outputs%argT
     &                  ,alp_nd,hnd,outputs%int_adm)

      end subroutine RTparam_outer
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine RTparam_fast(ac,modT,argR,argT,
     &    alp_nd,hnd,int_adm)

      real, intent(out)  :: ac,modT,argR,argT
      real, intent(out)  :: alp_nd,hnd,int_adm

      real   :: t_a,t_h,modR            !!inputs to RTparam_inner;
      real   :: alp_lin3,alp_lin4
      integer  :: OPT,LOW
      !!
      real   :: alp_nd_lims(6) =
     &   (/1.0e-6,0.005e0,0.3e0,1.5e0,2.5e0,3.5e0/)
      real   :: mc_alplin(2)   =
     &   (/-3.323529252398524e0,3.119943407349375e0/)
      real   :: y0_ll          = 40.0e0
      real   :: dy_ll          = 120.0e0
      integer  :: n_ll           = 3
      real   :: h1_ll          = .4e0
      
      real   :: hnd_lims(3)    = (/1.0e-2,0.2e0,0.4e0/)
      integer  :: LOG_A_vec(5)   = (/1,1,1,0,1/)
      integer  :: LOG_A
      real   :: a0,a1,h0,h1,l0,l1,dtmp,PI
      
      !!start of approximately linear regime;
      alp_lin3       = mc_alplin(2)+mc_alplin(1)*log(hnd)
      alp_nd_lims(5) = alp_lin3
      
      !!end of calculatable results;
      PI       = 4.0*atan(1.0)
      dtmp     = cos(hnd/h1_ll*PI/2.0)
      alp_lin4 = y0_ll+dy_ll*exp(n_ll*log(dtmp))
      
      alp_nd_lims(6) = alp_lin4
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!determine if low/high thickness regime;
      if (hnd<=hnd_lims(1)) then
         hnd    = hnd_lims(1) 
         LOW    = 1
      
      else if (hnd>=hnd_lims(3)) then
         hnd    = hnd_lims(3)
         LOW    = 0
      
      !!now between hnd_lims[0] and hnd_lims[1];
      else if (hnd<hnd_lims(2)) then
         LOW = 1
      else
         LOW = 0
      end if
      
      !!Interpolation value (-1<=t<=1) to use for hnd;
      h0  = hnd_lims(2-LOW)
      h1  = hnd_lims(3-LOW)
      
      if (LOW==0) then
         !!linear interpolation wrt hnd;
         t_h   = -1.0+2*(hnd-h0)/(h1-h0)
      else
         !!interpolation wrt log(hnd);
         !printf("\n(h0,h1)=(%f,%f)\n",h0,h1);
         l0  = log(h0)
         l1  = log(h1)
         t_h = -1.0+2.0*(log(hnd)-l0)/(l1-l0)
      end if
      !printf("\n(t_a,t_h)=(%f,%f)\n",t_a,t_h);
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!determine frequency regime;
      if (alp_nd<=alp_nd_lims(1)) then
         alp_nd = alp_nd_lims(1)
         OPT    = 1
      else if (alp_nd>=alp_nd_lims(6)) then
         alp_nd = alp_nd_lims(6)
         OPT    = 5
      !!
      else if (alp_nd<alp_nd_lims(2)) then
         OPT = 1
      else if (alp_nd<alp_nd_lims(3)) then
         OPT = 2
      else if (alp_nd<alp_nd_lims(4)) then
         OPT = 3
      else if (alp_nd<alp_nd_lims(5)) then
         OPT = 4
      else
         OPT = 5
      end if
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      a0     = alp_nd_lims(OPT)
      a1     = alp_nd_lims(OPT+1)
      LOG_A  = LOG_A_vec(OPT)
      if (LOG_A==1) then
         !!interpolation wrt log(alp_nd);
         l0  = log(a0)
         l1  = log(a1)
         t_a = -1.0+2.0*(log(alp_nd)-l0)/(l1-l0)
      else
         !!interpolation wrt alp_nd;
         t_a = -1.0+2.0*(alp_nd-a0)/(a1-a0)
      end if
        
      !!call main interpolation routine;
      call RTparam_inner(ac,modT,argR,argT,
     &              t_a,t_h,OPT,LOW,int_adm)

      end subroutine RTparam_fast
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine RTparam_inner(ac,modT,argR,argT,
     &                  tx,ty,OPT,LOW,int_adm)
      
      use mod_atten_young_hardcoded

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!declare "input/output variables";
      real, intent(out)   :: ac,modT,argR,argT
      real, intent(in)    :: tx,ty,int_adm
      integer, intent(in)   :: OPT,LOW
      !!
      
      !!auxiliary variables;
      real, allocatable  :: chebys(:)!!use allocate/deallocate later
      !real, allocatable  :: Zthings(:)!!use deallocate/deallocate later
      real               :: Zthings(4)!!use deallocate/deallocate later
      !!
      real    :: z,Rr,Ri,Tr,Ti,modR
      integer   :: Ncx,Ncy,ncol
      integer   :: i,j,r
      integer   :: Nthings
      integer   :: INTERP_MODE(5) = (/1,1,3,2,1/)
      integer   :: IM
      
      !!Order of chebyshev expansion for lower thicknesses;
      integer  :: Ncx_L(5) = (/10,10,10,10,3/) !!for alp_nd;
      integer  :: Ncy_L(5) = (/10,10,10,10,10/)!!for hnd;
      
      !!Order of chebyshev expansion for higher thicknesses;
      integer  :: Ncx_H(5) = (/10,10,10,10,4/) !!for alp_nd;
      integer  :: Ncy_H(5) = (/10,10,10,10,10/)!!for hnd;
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      !!Mode of interpolation;
      IM  =INTERP_MODE(OPT)
      if (IM==1) then
         !!Interpolate log(ac), Arg[R,T];
         Nthings   = 3
      else if (IM==2) then
         !!Interpolate (ac), Arg[R,T];
         Nthings   = 3
      else if (IM==3) then
         !!Interpolate Re[R,T], Im[R,T];
         Nthings   = 4
      end if
      
      !print*,'LOW,OPT=',LOW,OPT
      
      if (LOW==1) then
         !!Lower non-dimensional thicknesses;
         Ncx = Ncx_L(OPT)
         Ncy = Ncy_L(OPT)
         call get_Zthings(Zthings,tx,ty,Nthings,Ncx,Ncy,LOW,OPT)
      else
         !!Higher non-dimensional thicknesses;
         Ncx = Ncx_H(OPT)
         Ncy = Ncy_H(OPT)
         call get_Zthings(Zthings,tx,ty,Nthings,Ncx,Ncy,LOW,OPT)
      end if
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!interpret interpolated quantities;
      call RTparam_get_ac(ac,modT,argR,argT,Zthings,IM,int_adm)
      
      end subroutine
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine RTparam_get_ac(ac,modT,argR,argT,
     &      Zthings,IM,int_adm)
      
      real ,intent(out)  :: ac,modT,argR,argT
      real ,intent(in)   :: Zthings(4),int_adm
      integer,intent(in)   :: IM
      !!
      real   ::Rr,Ri,Tr,Ti,modR
      
      if (IM==1) then
         ac     = exp(Zthings(1))
         argR   = Zthings(2)
         argT   = Zthings(3)
         modT   = sqrt(exp(-ac/2.0)/int_adm)
      else if (IM==2) then
         ac     = Zthings(1)
         argR   = Zthings(2)
         argT   = Zthings(3)
         modT   = sqrt(exp(-ac/2.0)/int_adm)
      else if (IM==3) then
         Rr     = Zthings(1)
         Ri     = Zthings(2)
         Tr     = Zthings(3)
         Ti     = Zthings(4)
         !!
         argR   = atan2(Ri,Rr)
         argT   = atan2(Ti,Tr)
         modR   = sqrt(Rr*Rr+Ri*Ri)
         modT   = sqrt(Tr*Tr+Ti*Ti)
         ac     = -2*log(1-modR*modR)
      end if
      
      end subroutine RTparam_get_ac
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      !!%%function [ki,BG2,coeff_del,coeff_D]=gen_root_ice(del,H,guess)
      !!finds the root of the ice dispersion relation nearest to 'guess'.
      subroutine gen_root_ice(ki,BG2,coeff_del,coeff_D
     &                ,del, H,guess)

      real,intent(out)  :: ki,BG2,coeff_del,coeff_D
      real,intent(in)   :: del,H,guess
      
      real     :: fac,k0,dk,tol
      real     :: res,Lam,Lampr,denom
      integer  :: its,MAXITS
      
      fac   = 1.0
      k0    = guess
      tol   = 1.0e-7!relative tolerance
      
      !!Call dispersion relation function;
      !!dk    = NR_corr_term(k0,del,H,fac);
      call NR_corr_term(dk,Lam,Lampr,k0,del,H,fac)
      ki = k0-dk
      MAXITS   = 15
      its   = 1
      do while((abs(dk/ki).ge.tol).and.(its.lt.MAXITS))
         k0 = ki
        
         !!Call dispersion relation function;
         !!dk = NR_corr_term(k0,del,H,fac);
         call NR_corr_term(dk,Lam,Lampr,k0,del,H,fac)
         ki    = k0-dk
         its   = its+1
      end do
        
      !!Call dispersion relation function;
      !![dk,Lam,Lampr] = NR_corr_term(ki,del,H,fac);
      call NR_corr_term(dk,Lam,Lampr,ki,del,H,fac)
      !!
      denom = H*(Lam*Lam*ki*ki-1.0)+Lampr
      res   = -ki/denom
      
      !!Outputs;
      BG2  = Lam*Lam*res
      coeff_del  = -ki/denom   !! dk/d(del) - non-dimensional
      coeff_D    = -ki**5/denom!! dk/d(D)   - non-dimensional
      end subroutine gen_root_ice
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      !!function [kw,BG1]=gen_root_wtr(del,H,guess)
      !!%% finds the root of the water dispersion relation nearest to 'guess'.
      subroutine gen_root_wtr(kw,BG1,
     &                 del, H,guess)
      
      real,intent(out)   :: kw,BG1
      real,intent(in)    :: del,H,guess
      !!
      real     :: fac,k0,dk,tol
      real     :: res,Lam,Lampr,denom
      integer  :: its,MAXITS
      
      tol   = 1.0e-7!relative tolerance
      fac   = 0.0
      k0    = guess
      
      !!Call dispersion relation function;
      !!dk    = NR_corr_term(k0,del,H,fac);
      call NR_corr_term(dk,Lam,Lampr,k0,del,H,fac)
      kw = k0-dk

      MAXITS   = 15
      its      = 1
      
      do while((abs(dk/kw).gt.tol).and.(its.le.MAXITS))
         k0 = kw
      
         !!Call dispersion relation function;
         !!dk       = NR_corr_term(k0,del,H,fac);
         call NR_corr_term(dk,Lam,Lampr,k0,del,H,fac)
         kw    = k0-dk
         its   = its+1
      end do
      
      !!Call dispersion relation function;
      !![dk,Lam,Lampr] = NR_corr_term(kw,del,H,fac);
      call NR_corr_term(dk,Lam,Lampr,kw,del,H,fac)
      denom          = H*(Lam*Lam*kw*kw-1.0)+Lampr
      res            = -kw/denom
      !!printf("(denom,res,Bg2)=(%.9e,%.9e,%.9e)",res,denom,Lam*Lam*res);
      
      !!Outputs;
      BG1  = Lam*Lam*res
      end subroutine gen_root_wtr
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      !! function [dk,Lam,Lampr] = NR_corr_term(k,del,H,fac)
      !! %% dk=f/f_k, where f has the same zeros as of the dispersion function, 
      !! %% is the correction term in the Newton-Rhapson method for finding zeros in f.
      subroutine NR_corr_term(dk,Lam,Lampr,
     &                 k, del, H,fac)
      
      real,intent(out)   :: dk,Lam,Lampr
      real,intent(in)    :: k,del,H,fac
      !!
      real   :: f,df,x,k4
      
      k4    = k*k*k*k
      Lam   = fac*k4+del
      Lampr = 5.0d0*fac*k4+del
      x     = 7.5d0
      
      if (abs(k*H)<=x) then
         f  = Lam*k*sinh(k*H)-cosh(k*H)
         df = Lam*(k*H)*cosh(k*H)+(Lampr-H)*sinh(k*H)
      else
         f  = Lam*k*tanh(k*H)-1
         df = Lam*k*H+(Lampr-H)*tanh(k*H)
      end if
      
      !!Outputs;
      dk      = f/df
      
      end subroutine NR_corr_term
      end module mod_atten_young
