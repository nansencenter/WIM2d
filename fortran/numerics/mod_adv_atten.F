      module mod_adv_atten
      implicit none
      contains

      subroutine adv_atten_simple(Sdir                                  &
     &   ,Sfreq,taux_om,tauy_om                                         &
     &   ,atten_dim,damp_dim,ag_eff,dt)

      use mod_gridinfo
      use mod_common_wim
      use mod_advect

      implicit none

      integer,parameter :: ndir     = n_wavdir
      integer,parameter :: ADV_OPT  = 2
         !1: periodic in x,y; 2: periodic in y only; 0: not periodic
      !!
      real,dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,ndir)
     &   ,intent(inout)  :: Sdir
      real,dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)                   &
     &   ,intent(out)    :: Sfreq,taux_om,tauy_om
      real,dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)                   &
     &   ,intent(in)     :: ag_eff,atten_dim,damp_dim
      real,intent(in)   :: dt
      !!
      integer  :: i,j,wth
      real,dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::                &
     &   uwave,vwave
      real  :: PI,adv_dir,wt_theta(ndir),S_th,tmp,alp_dim,source

      PI          = 4.0*atan(1.0)

      !!do advection
      print*,'Spectrum at',itest,jtest
      print*,'Before advection',Sdir(itest,jtest,1)
      do wth  = 1,ndir
         !!define wave direction & speeds
         adv_dir  = -PI/180.0*(90.0+wavdir(wth))
         uwave    = ag_eff*cos(adv_dir)
         vwave    = ag_eff*sin(adv_dir)

         !!advect
         call waveadv_weno(Sdir(:,:,wth),uwave,vwave,dt,ADV_OPT)
      end do
      !print*,'ag,adv_dir      ',ag_eff(itest,jtest),adv_dir
      print*,'u,v of wave     ',uwave(itest,jtest),vwave(itest,jtest)
      print*,'After  advection',Sdir(itest,jtest,1)
      print*,'Attenuation     ',atten_dim(itest,jtest)                  &
     &                         ,damp_dim(itest,jtest)

      !!weights for integral over direction
      if (ndir.eq.1) then
         wt_theta = 1.0
      else
         wt_theta = 2*PI/(1.0*ndir)
      end if

      Sfreq    = 0.0
      taux_om  = 0.0
      tauy_om  = 0.0

      do j = 1,jj
      do i = 1,ii
         if (ICE_MASK(i,j).gt.0.0) then
            !! atten_dim = ENERGY attenuation coeff [m^{-1}]
            do wth=1,ndir
               adv_dir  = -PI/180.0*(90.0+wavdir(wth))
               S_th     = Sdir(i,j,wth)
               alp_dim  = atten_dim(i,j)+damp_dim(i,j)

               !!stress calculation
               source         = -alp_dim*ag_eff(i,j)*S_th!! m^{-1}*[m/s]*[m^2s] = m^2
               tmp            = -cos(adv_dir)*wt_theta(wth)*source
               taux_om(i,j)  = taux_om(i,j)+tmp
               tmp            = -sin(adv_dir)*wt_theta(wth)*source
               tauy_om(i,j)  = tauy_om(i,j)+tmp
               !! tau_x,tau_y need to be multiplied by rho_wtr*g/phase_vel
               !!  and integrated over frequency as well;
               !! units: [m^2]*[kg/m^3]*[m/s^2]*[s/m]*s^{-1}
               !!         = kg/m/s^2 = Pa
      
               !! NB take '-' because here we have calc'd the stress
               !! ON the waves (cf Donelan et al, 2012, JGR)
               !! - we want the stress on the ice

               !!do attenuation
               Sdir(i,j,wth) = S_th*exp(-alp_dim*ag_eff(i,j)*dt)
            end do
      
         end if!!end ice check (attenuation and stress calc)
      
         !! INTEGRATE SPECTRUM OVER DIRECTION;
         !! NB do this in water also
         do wth=1,ndir
            Sfreq(i,j)  = Sfreq(i,j)+wt_theta(wth)*Sdir(i,j,wth)
         end do
      end do!i
      end do!j
      
      end subroutine adv_atten_simple
      end module mod_adv_atten
