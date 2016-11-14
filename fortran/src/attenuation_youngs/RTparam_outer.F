      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine atten_youngs_f2py(outputs,inputs)
      use mod_common_wim, only: wim_info
      !wrapper for f2py (make sure doubles are changed to reals if necessary)

      real*8,intent(in)   :: inputs(5)  !(/h,om,young,drag_rp,visc_ws/)
      real*8,intent(out)  :: outputs(8) !(/damping,kice,kwtr,int_adm,ac,modT,argR,argT/)

      real  :: damping,kice,kwtr,int_adm,ac,modT,argR,argT

      real  :: guess,om

      !!set some params manually
      wim_info%poisson  = 0.3
      wim_info%rhowtr   = 1025.
      wim_info%gravity  = 9.81
      om                = real(inputs(2))
      wim_info%young    = real(inputs(3))
      wim_info%drag_rp  = real(inputs(4))
      wim_info%visc_ws  = real(inputs(5))
      guess             = om**2/wim_info%gravity

      call RTparam_outer(damping,kice,kwtr,int_adm,
     &                  ac,modT,argR,argT,
     &                  real(inputs(1:2)),guess)

      outputs  = dble((/damping,kice,kwtr,int_adm
     &                 ,ac,modT,argR,argT/))

      end subroutine atten_youngs_f2py
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!function [damping,kice,kwtr,int_adm,...
      !!            alp_scat,modT,argR,argT] =...
      !!               RT_param_outer(h,om,E,drag_rp,guess)
      subroutine RTparam_outer(damping,kice,kwtr,int_adm,
     &                  ac,modT,argR,argT,
     &                  inputs,guess)

      use mod_common_wim, only:wim_info

      real,intent(out)   :: damping,kice,kwtr,int_adm,
     &                           ac,modT,argR,argT
      real,intent(in)    :: inputs(2),guess
      real               :: h,om
      
      real   :: alp_nd,hnd,zeta_nd
      real   :: Hw_nd,varpi,L
      real   :: damping_nd,drag_rp_nd
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
      kice = ki/L
      !print*,'vHkGB_i',varpi,H_nd,kice,guess*L,BG1
         
      !!get wavenumber for water;
      varpi = 1.0/alp_nd
      Hw_nd = H_nd+zeta_nd
      !![kw,BG1] = gen_root_wtr(varpi,Hw_nd,alp_nd);
      call gen_root_wtr(kw,BG1,varpi,Hw_nd,alp_nd)
      kwtr  = kw/L
      !print*,'vHkGB_w',varpi,Hw_nd,kwtr,alp_nd,BG1
         
      !!get intrinsic admittance;
      !!|R|^2+int_adm*|T|^2=1
      !print*,'BG1   = ', BG1
      !print*,'BG2   = ', BG2
      int_adm  = BG1/BG2
         
      !!get dissipative attenuation;
      drag_rp_nd  = wim_info%drag_rp/wim_info%rhowtr/om/L !!RP drag
      vews_nd     = h**3*rho*wim_info%visc_ws/
     &               (6.*(1.+wim_info%poisson)*om*L**5)
      damping_nd  = -coeff_del*drag_rp_nd -coeff_D*vews_nd
      damping     = damping_nd/L
      
      !!get interpolated attenuation coefficient;
      call RTparam_fast(ac,modT,argR,argT,
     &             alp_nd,hnd,int_adm)

      end subroutine RTparam_outer
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
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
