      module mod_adv_atten
      implicit none
      contains

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine adv_atten_simple(Sdir                                  &
     &   ,Sfreq,taux_om,tauy_om                                         &
     &   ,atten_dim,damp_dim,ag_eff,dt,ADV_DIM,ADV_OPT)

      use mod_gridinfo
      use mod_common_wim, only: dfloe,dfloe_pack_init,ICE_MASK          &
     &   ,wavdir,n_wavdir
      use mod_advect

      implicit none

      integer,parameter :: ndir     = n_wavdir
      !!
      real,dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,ndir)
     &   ,intent(inout)    :: Sdir
      real,dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)                   &
     &   ,intent(out)      :: Sfreq,taux_om,tauy_om
      real,dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)                   &
     &   ,intent(in)       :: ag_eff,atten_dim,damp_dim
      real,intent(in)      :: dt
      integer,intent(in)   :: ADV_DIM,ADV_OPT
      !!
      integer  :: i,j,wth
      real,dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::                &
     &   uwave,vwave
      real  :: PI,adv_dir,wt_theta(ndir),S_th,tmp,alp_dim,source

      PI          = 4.0*atan(1.0)

      !!do advection
      !print*,'Spectrum at',itest,jtest
      !print*,'Before advection',Sdir(itest,jtest,1)
      uwave = 0.0
      vwave = 0.0

      !print*,'max= ',maxval(ag_eff)
      !print*,'min= ',minval(ag_eff)

      do wth  = 1,ndir
         !!define wave direction & speeds
         adv_dir  = -PI/180.0*(90.0+wavdir(wth))
         !print*,'adv_dir= ',cos(adv_dir)
         uwave    = ag_eff*cos(adv_dir)

         if (ADV_DIM.eq.2) then
            !! normal 2d advection
            vwave = ag_eff*sin(adv_dir)

            !if (wth.eq.4) then
             !  print*,'min 1= ',minval(Sdir(:,:,wth))
              ! print*,'max 1= ',maxval(Sdir(:,:,wth))
            !end if

            !!advect
            call waveadv_weno(Sdir(:,:,wth),uwave,vwave,dt,ADV_OPT)

            !if (wth.eq.4) then
            !   print*,'min 2= ',minval(Sdir(:,:,wth))
            !   print*,'max 2= ',maxval(Sdir(:,:,wth))
            !end if

         else
            !! * in some situations where no variation is expected
            !! in the 'y' dirn, we can ignore the speed in this
            !! dirn
            !! * this can also be used to do multiple 1d results at once
            !!   - ie each row can be a different 1d expt

            !!advect
            call waveadv_weno_1d(Sdir(:,:,wth),uwave,dt,ADV_OPT)
         end if
      end do
!     print*,'ag,adv_dir      ',ag_eff(itest,jtest),adv_dir
!     print*,'u,v of wave     ',uwave(itest,jtest),vwave(itest,jtest)
!     print*,'After  advection',Sdir(itest,jtest,1)
!     print*,'Attenuation     ',atten_dim(itest,jtest)                  &
!    &                         ,damp_dim(itest,jtest)

      !!weights for integral over direction
      if (ndir.eq.1) then
         wt_theta = 1.0
      else
         wt_theta = 2*PI/(1.0*ndir)
      end if

      Sfreq    = 0.0
      taux_om  = 0.0
      tauy_om  = 0.0

      print*,'min= ',minval(damp_dim)
      print*,'max= ',maxval(damp_dim)


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
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine adv_atten_isotropic(Sdir                               &
     &   ,Sfreq,taux_om,tauy_om                                         &
     &   ,atten_dim,damp_dim,ag_eff,dt,ADV_DIM,ADV_OPT)

      use mod_gridinfo
      use mod_common_wim, only: dfloe,dfloe_pack_init,ICE_MASK          &
     &   ,wavdir,n_wavdir
      use mod_advect

      implicit none

      integer,parameter :: ndir     = n_wavdir
      integer,parameter :: ncs      = nint(ndir/2.0)!!number of sine/cosine coeff's
      !!
      real,dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,ndir)
     &   ,intent(inout)    :: Sdir
      real,dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)                   &
     &   ,intent(out)      :: Sfreq,taux_om,tauy_om
      real,dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)                   &
     &   ,intent(in)       :: ag_eff,atten_dim,damp_dim
      real,intent(in)      :: dt
      integer,intent(in)   :: ADV_DIM,ADV_OPT
      !!
      integer                 :: i,j,wth,nth,n,jp1,jm1,ss
      integer,dimension(ndir) :: nvec
      real,dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::                &
     &   uwave,vwave
      real  :: PI,adv_dir,tmp,alp_dim,source

      real,dimension(ndir)    :: K_fou,S_th,theta_vec,wt_theta,         &
     &                              tmp1,evals_x
      complex,dimension(ndir) :: S_fou
      complex                 :: zi,src_fou_p1,src_fou_m1
      real,dimension(ncs)     :: S_cos,S_sin
      real  :: cg,q_scat,q_abs,q_tot,src_cos_1,src_sin_1
      integer :: cpt = 0

      PI = 4.0*atan(1.0)
      !zi = cmplx(1,0)!sqrt(-1)
      zi = cmplx(0,1)!sqrt(-1)

      !!weights for integral over direction
      wt_theta = 2*PI/(1.0*ndir) !!\Delta\theta

      !!do advection
      !print*,'Spectrum at',itest,jtest
      !print*,'Before advection',Sdir(itest,jtest,1)
      uwave = 0.0
      vwave = 0.0
      do wth  = 1,ndir
         !!define wave direction & speeds
         adv_dir  = -PI/180.0*(90.0+wavdir(wth))
         uwave    = ag_eff*cos(adv_dir)

         if (ADV_DIM.eq.2) then
            !! normal 2d advection
            vwave = ag_eff*sin(adv_dir)

!!advect

            !print*,'Min 1= ',minval(Sdir(:,:,wth))
            !print*,'Max 1= ',maxval(Sdir(:,:,wth))

            call waveadv_weno(Sdir(:,:,wth),uwave,vwave,dt,ADV_OPT)

            !print*,'Min 2= ',minval(Sdir(:,:,wth))
            !print*,'Max 2= ',maxval(Sdir(:,:,wth))


         else
            !! * in some situations where no variation is expected
            !! in the 'y' dirn, we can ignore the speed in this
            !! dirn
            !! * this can also be used to do multiple 1d results at once
            !!   - ie each row can be a different 1d expt

            !!advect
            call waveadv_weno_1d(Sdir(:,:,wth),uwave,dt,ADV_OPT)
         end if

         !!useful later
         theta_vec(wth) = adv_dir

         print*,'theta_vec[',wth-1,']= ',theta_vec(wth)

         nvec(wth)      = wth-1
      end do
!     print*,'ag,adv_dir      ',ag_eff(itest,jtest),adv_dir
!     print*,'u,v of wave     ',uwave(itest,jtest),vwave(itest,jtest)
!     print*,'After  advection',Sdir(itest,jtest,1)
!     print*,'Attenuation     ',atten_dim(itest,jtest)                  &
!    &                         ,damp_dim(itest,jtest)

      Sfreq    = 0.0
      taux_om  = 0.0
      tauy_om  = 0.0

      cpt = 0

      do j = 1,jj
      do i = 1,ii

         do wth=1,ndir
            S_th(wth)      = Sdir(i,j,wth)
         end do
         S_fou    = zi

         !! S_fou(1) = \int_0^{2\pi}S(\theta)d\theta
         S_fou(1) = cmplx( sum(wt_theta*S_th) ) !n=0 Fourier coefficient

         !if ((i.eq.1).and.(j.eq.1)) then
         !if (cpt.eq.0) then
         !   print*,'S_fou R= ',realpart(S_fou(1))
         !   print*,'S_fou I= ',imagpart(S_fou(1))
         !end if

         !cpt = cpt+1

         if (ICE_MASK(i,j).gt.0.0) then

            cpt = cpt+1

            !! atten_dim,damp_dim = ENERGY attenuation coeffs [m^{-1}]
            if (dfloe(i,j).lt.dfloe_pack_init) then
               !! only apply scattering in the MIZ
               !! - otherwise ice is completely broken
               q_scat   = atten_dim(i,j)
               q_abs    = damp_dim(i,j)
            else
               !! in pack all energy is dissipated
               q_scat   = 0.0
               q_abs    = atten_dim(i,j)+damp_dim(i,j)
            end if
            q_tot    = q_scat+q_abs
            cg       = ag_eff(i,j)

!           !!uncomment below for testing against WIM2d.m
!           !!NB !!DANGER!! !!overwriting actual values!!
!           !!test case = waves(12s,2m),ice(2m,70%),1 freq, 16 dirs
!           if ((i.eq.24).and.(j.eq.1)) then
!              print*,'NB !!DANGER!! !!'//
!    &                'overwriting actual values for testing!!'
!              q_scat   = 8.041601e-06
!              q_tot    = 2.736895e-05
!              q_abs    = q_tot-q_scat
!              !!
!              S_th(1 ) = 1.348636e-07
!              S_th(2 ) = 1.443994e-07
!              S_th(3 ) = 1.650405e-07
!              S_th(4 ) = 1.985781e-07
!              S_th(5 ) = 2.397257e-07
!              S_th(6 ) = 2.287021e-07
!              S_th(7 ) = 1.746545e-05
!              S_th(8 ) = 2.584725e-04
!              S_th(9 ) = 5.674894e-04
!              S_th(10) = 2.584722e-04
!              S_th(11) = 1.746496e-05
!              S_th(12) = 2.276650e-07
!              S_th(13) = 2.384907e-07
!              S_th(14) = 1.976913e-07
!              S_th(15) = 1.645170e-07
!              S_th(16) = 1.439858e-07
!              !!
!              S_fou(1) = cmplx( sum(wt_theta*S_th) )
!           end if

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !! - kernel of Boltzmann eqn
            !! - eigenvalues

            !! In general K_fou will be a cosine series
            !! (with coefficients K_cos)
            !! since K will be of the form K(\theta-\theta')
            !! => K_fou will be real (with K_fou(1) = K_cos(1),
            !! &  K_fou(n)=K_fou(ndir+1-n)=K_cos(n)/2.0, n>1 )
            !! For Meylan & Masson (2006) model, pass in K_cos
            K_fou    = 0.0
            K_fou(1) = q_scat

            !! evals_x  = K_fou-q_tot
            !! NB these solve \pa_x.S=evals_x.S
            !! possibly need to be careful to avoid rounding errors
            !! making eigenvalues go >0
            evals_x     = -q_tot
            evals_x(1)  = -q_abs
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            !!S_fou(n) = s_n = \int_0^{2\pi}S(\theta)e^{i*n*\theta}d\theta
            do nth=1,ncs

               !wt_theta = 1.0
               !theta_vec = 1.0
               S_cos(nth)     = sum(wt_theta*S_th*cos(nth*theta_vec))
               S_sin(nth)     = sum(wt_theta*S_th*sin(nth*theta_vec))

               !if ((cpt.eq.1).and.(nth.eq.1)) then
               if ((cpt.eq.1)) then
                  !print*,'ncs  = ',ncs
                  !print*,'S_th= ',sum(S_th)
                  !print*,'Real= ',S_cos(nth)
!print*,'Imag= ',S_sin(nth)
                  print*,'S_th= ',S_th(1)
                  print*,'S_th= ',S_th(2)
               end if


               !! 1<=n<=N/2: e^{i*n*\theta}=\cos(n\theta)+i*\sin(n\theta)
               S_fou(nth+1)   = cmplx(S_cos(nth),S_sin(nth))

               !!also need n=-nth index
               !! N/2<n<=N (n<0): e^{i*n*\theta}=\cos(n\theta)-i*\sin(n\theta)
               n           = ndir+1-nth
               S_fou(n)    = cmplx(S_cos(nth),-S_sin(nth))
                  !! NB don't need nth=ncs here,
                  !! but if \theta_m=2\pi(m/N),
                  !! e^{iN/2\theta_m}=e^{im\pi}=e^{-im\pi}=e^{-iN/2\theta_m}
            end do

            !!stresses
            jp1   = 2
            jm1   = ndir
            !!
            src_fou_p1  = cg*( -q_tot*S_fou(jp1)+K_fou(jp1)*S_fou(jp1) )
            src_fou_m1  = cg*( -q_tot*S_fou(jm1)+K_fou(jm1)*S_fou(jm1) )
            src_cos_1   = real(     .5*(src_fou_p1+src_fou_m1) )
            src_sin_1   = real( -zi*.5*(src_fou_p1-src_fou_m1) )
            !!
            taux_om(i,j)   = -src_cos_1
            tauy_om(i,j)   = -src_sin_1

            !if ((cpt.eq.1)) then
               !print*,'ncs  = ',ncs
               !print*,'S_th= ',sum(S_th)
            !   print*,'Real= ',-src_cos_1
            !   print*,'Imag= ',-src_sin_1
            !end if


!           !!uncomment below for testing against WIM2d.m
!           if ((i.eq.24).and.(j.eq.1)) then
!              print*,'Hs',4.0*sqrt(real(S_fou(1)))                     &
!    &               ,4.0*sqrt(sum(wt_theta*S_th))
!              print*,'q_scat,q_abs,q_tot',q_scat,q_abs,q_tot
!              print*,'theta,S_th,S_fou,K_fou'
!              do wth=1,ndir
!                 print*,theta_vec(wth),S_th(wth),S_fou(wth),K_fou(wth)
!              end do
!              print*,'tau_x,tau_y',taux_om(i,j),tauy_om(i,j)
!           end if

            !!do scattering/attenuation and transform back
            S_fou = S_fou*exp(evals_x*cg*dt)

            !if ((cpt.eq.1)) then
            !   !print*,'S_th= ',sum(S_th)
            !   print*,'***S_fou R= ',realpart(S_fou(1))
            !   print*,'***S_fou I= ',imagpart(S_fou(1))
            !end if


            !!S(\theta)=(1/2/\pi)*\sum_n.S_n.e^{-i*n*\theta}
            do wth=1,ndir
!              tmp1  = real(    S_fou)*cos(nvec*theta_vec(wth))         &
!    &                -real(-zi*S_fou)*sin(nvec*theta_vec(wth))
               tmp1  = realpart(S_fou)*cos(nvec*theta_vec(wth))         &
     &                -imagpart(S_fou)*sin(nvec*theta_vec(wth))
!     !
               !tmp1=1
               Sdir(i,j,wth)  = 1.0/(2*PI)*sum(tmp1)

             !  if ((cpt.eq.1)) then
             !     do ss=1,ndir
              !       print*,'tmp1[',ss-1,']= ',tmp1(ss)
              !    end do
              ! end if


               !if ((cpt.eq.1).and.(wth.eq.1)) then
               if ((cpt.eq.1).and.(wth.eq.2)) then

                  !do ss=1,ndir
                  !   print*,'tmp1[',ss-1,']= ',tmp1(ss)
                  !end do

                  !print*,'***Temp R= ',realpart(S_fou(1))
!     print*,'***Temp I= ',imagpart(S_fou(1))
                  !print*,'theta[',wth-1,']=',theta_vec(wth)
                  !print*,'Sdir[',wth-1,']=',Sdir(i,j,wth)
               end if

               !if ((cpt.eq.1).and.(wth.eq.2)) then
               !   print*,'Min OUT= ',minval(Sdir)
               !   print*,'Max OUT= ',maxval(Sdir)
               !end if



               S_th(wth)      = Sdir(i,j,wth)
            end do

!           !!uncomment below for testing against WIM2d.m
!           if ((i.eq.24).and.(j.eq.1)) then
!              print*,'Hs after',4.0*sqrt(real(S_fou(1)))               &
!    &                          ,4.0*sqrt(sum(wt_theta*S_th))
!              print*,'eigenvalues, S_th, S_fou after atten'
!              do wth=1,ndir
!                 print*,dd(wth),S_th(wth),S_fou(wth)
!              end do
!           end if

         end if!!end ice check (attenuation and stress calc)

         !! INTEGRATE SPECTRUM OVER DIRECTION;
         !! already when getting the Fourier transform
         !! NB do this in water also
         Sfreq(i,j)  = real(S_fou(1))

         !if ((cpt.eq.1)) then
         !   print*,'Sfreq[',i,j,']=',Sfreq(i,j)
!end if

         !print*,'Sfreq[',i,j,']=',Sfreq(i,j)

      end do!i
      end do                    !j

!print*,'cpt= ',cpt

      !print*,'Min OUT= ',minval(Sdir)
      !print*,'Max OUT= ',maxval(Sdir)

      print*,'Min OUT= ',minval(Sfreq)
      print*,'Max OUT= ',maxval(Sfreq)


      end subroutine adv_atten_isotropic

      end module mod_adv_atten
