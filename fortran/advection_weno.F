      subroutine iceadv_v2(h,u,v,scuy,scvx,scp2i,scp2,dt)
!
! --- ------------------------------------------------------------------
! --- Advection is done with flux limited 3rd order WENO in space and
! --- 2nd order Runge-Kutta in time
! --- ------------------------------------------------------------------
!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! modified [TW 28.1.2014] for wave advection
      !! - iceadv builds up ice (& thus waves) at coast;
      !! - iceadv_v2 lets waves go into/out of land,
      !!   then cancels the waves on land at end of routine;
      !! - just adds the wave cancellation to Luke's original fix
      !!   for waves on boundaries;
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! ARGUMENTS:
      !! in/out:  h is qty to be advected;
      !!
      !! in:      u,v,dt are speed components and time step;
      !! in:      scuy,scvx is mesh size at u (v) points in y (x) direction
      !!           - see common_blocks.h;
      !! in:      scp2, scp2i are grid box area at p points, and its inverse;
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use mod_xc
!
      implicit none
!
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::               &
     &  h,u,v,scuy,scvx,scp2i,scp2
      real dt
!
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: sao,hp
      real dtm
      integer i,j,l

       sao(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)=0
!
! --- Use a modified time step since velocities are in m/s while scale
! --- factors are in cm
!      dtm=dt*1.e2
       dtm=dt
! --- Prediction step
      call weno3pd_v2(h,sao,u,v,scuy,scvx,scp2i,scp2,dtm)
!$OMP PARALLEL DO
       margin=nbdy 
      do j=1-margin,jj+margin
!       do l=1,isp(j)
!       do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
        do i=1-margin,ii+margin!![TW 28.1.2014]
          hp(i,j)=h(i,j)+dtm*sao(i,j)
        enddo
!       enddo
!       enddo
      enddo
!$OMP END PARALLEL DO
!
! --- Correction step
      call weno3pd_v2(hp,sao,u,v,scuy,scvx,scp2i,scp2,dtm)
!$OMP PARALLEL DO
      do j=1-margin,jj+margin
!       do l=1,isp(j)
!       do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
        do i=1-margin,ii+margin!![TW 28.1.2014]
          if (LANDMASK(i,j).eq.0.0) then
             !!if not on land keep as it is;
             h(i,j)  = .5*(h(i,j)+hp(i,j)+dtm*sao(i,j))
          else
             !!if on land set quantitity (waves) to zero
             h(i,j)  = 0.0
          end if
        enddo
!       enddo
!       enddo
      enddo
!$OMP END PARALLEL DO
      return
      end subroutine iceadv_v2


       subroutine weno3pd_v2(g,sao,u,v,scuy,scvx,scp2i,scp2,dt)
!
! --- ------------------------------------------------------------------
! --- By a weighted essentially non-oscillatory scheme with up to 3th
! --- order accuracy, obtain the spatial advective operator of a
! --- 2-dimensional field defined at the scalar points of a C-grid. The
! --- fluxes are limited to make the scheme positive definite.
! --- Advective velocities in the i- and j-direction are defined at u-
! --- and v-points, respectively.
! --- ------------------------------------------------------------------
!
      use mod_xc
!
      implicit none
!
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::               &
     &  g,sao,u,v,scuy,scvx,scp2i,scp2
      real dt
!
      real cq00,cq01,cq10,cq11,ca0,ca1,eps
      parameter (cq00=-1./2.,cq01= 3./2.,                               &
     &           cq10= 1./2.,cq11= 1./2.,                               &
     &           ca0=1./3.,ca1=2./3.,                                   &
     &           eps=1.e-12)
!
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::               &
     &  ful,fuh,fvl,fvh,gt
      real q0,q1,a0,a1,q
      integer i,j,l,im1,im2,ip1,jm1,jm2,jp1
!
! --- Compute grid cell boundary fluxes. Split in a low order flux
! --- (donor cell) and a high order correction flux.
!
!$OMP PARALLEL DO
      do j=0,jj+2
        do i=0,ii+2
          ful(i,j)=0.
          fuh(i,j)=0.
          fvl(i,j)=0.
          fvh(i,j)=0.
        enddo
      enddo
!$OMP END PARALLEL DO
!
!      call xctilr(g, 1,1, 3,3, halo_ps)
!
!$OMP PARALLEL DO PRIVATE(im1,im2,q0,q1,a0,a1,ip1)
      do j=0,jj+1
!       do l=1,isu(j)
!       do i=max(0,ifu(j,l)),min(ii+2,ilu(j,l))
        do i=0,ii+2!![TW 28.1.2014]
          im1=i-1
!
          if (u(i,j).gt.0.) then
             im2=im1-iu(im1,j)
!
            q0=cq00*g(im2,j)+cq01*g(im1,j)
            q1=cq10*g(im1,j)+cq11*g(i  ,j)
!
            a0=ca0
            a1=ca1*(abs(g(im2,j)-g(im1,j))+eps)                         &
     &            /(abs(g(im1,j)-g(i  ,j))+eps)
!
            ful(i,j)=u(i,j)*g(im1,j)*scuy(i,j)
!

          else
            ip1=i+iu(i+1,j)
!
            q0=cq11*g(im1,j)+cq10*g(i  ,j)
            q1=cq01*g(i  ,j)+cq00*g(ip1,j)
!
            a0=ca1
            a1=ca0*(abs(g(im1,j)-g(i  ,j))+eps)                         &
     &            /(abs(g(i  ,j)-g(ip1,j))+eps)
!
            ful(i,j)=u(i,j)*g(i  ,j)*scuy(i,j)
!
          endif
!
          fuh(i,j)=u(i,j)*(a0*q0+a1*q1)/(a0+a1)*scuy(i,j)-ful(i,j)
!
        enddo
!       enddo
!       enddo
      enddo
!$OMP END PARALLEL DO
!
!$OMP PARALLEL DO PRIVATE(jm1,q0,q1,a0,a1,jm2,jp1)
      do j=0,jj+2
        jm1=j-1
!       do l=1,isv(j)
!       do i=max(0,ifv(j,l)),min(ii+1,ilv(j,l))
        do i=0,ii+1!![TW 28.1.2014]
!
          if (v(i,j).gt.0.) then
            jm2=jm1-iv(i,jm1)
!
            q0=cq00*g(i,jm2)+cq01*g(i,jm1)
            q1=cq10*g(i,jm1)+cq11*g(i,j  )
!
            a0=ca0
            a1=ca1*(abs(g(i,jm2)-g(i,jm1))+eps)                         &
     &            /(abs(g(i,jm1)-g(i,j  ))+eps)
!
            fvl(i,j)=v(i,j)*g(i,jm1)*scvx(i,j)
!
          else
            jp1=j+iv(i,j+1)
!
            q0=cq11*g(i,jm1)+cq10*g(i,j  )
            q1=cq01*g(i,j  )+cq00*g(i,jp1)
!
            a0=ca1
            a1=ca0*(abs(g(i,jm1)-g(i,j  ))+eps)                         &
     &            /(abs(g(i,j  )-g(i,jp1))+eps)
!
            fvl(i,j)=v(i,j)*g(i,j  )*scvx(i,j)
!
          endif
!
          fvh(i,j)=v(i,j)*(a0*q0+a1*q1)/(a0+a1)*scvx(i,j)-fvl(i,j)
!        
        enddo
!       enddo
!       enddo
      enddo
!$OMP END PARALLEL DO
!
! --- Update field with low order fluxes.
!$OMP PARALLEL DO
      do j=0,jj+1
!       do l=1,isp(j)
!       do i=max(0,ifp(j,l)),min(ii+1,ilp(j,l))
        do i=0,ii+1!![TW 28.1.2014]
          gt(i,j)=g(i,j)-dt*(ful(i+1,j)-ful(i,j)                        &
     &                      +fvl(i,j+1)-fvl(i,j))*scp2i(i,j)
        enddo
!       enddo
!       enddo
      enddo
!$OMP END PARALLEL DO
!
! --- Obtain fluxes with limited high order correction fluxes.
      q=.25/dt
!$OMP PARALLEL DO
      do j=1,jj
!       do l=1,isu(j)
!       do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l))
        do i=1,ii+1!![TW 28.1.2014]
          fuh(i,j)=ful(i,j)+max(-q*gt(i  ,j)*scp2(i  ,j),               &
     &                      min( q*gt(i-1,j)*scp2(i-1,j),fuh(i,j)))
        enddo
!       enddo
!       enddo
      enddo
!$OMP END PARALLEL DO
!$OMP PARALLEL DO
      do j=1,jj+1
!       do l=1,isv(j)
!       do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
        do i=1,ii!![TW 28.1.2014]
          fvh(i,j)=fvl(i,j)+max(-q*gt(i,j  )*scp2(i,j  ),               &
     &                      min( q*gt(i,j-1)*scp2(i,j-1),fvh(i,j)))
        enddo
!       enddo
!       enddo
      enddo
!$OMP END PARALLEL DO
!
! --- Compute the spatial advective operator.
!$OMP PARALLEL DO
      do j=1,jj
!       do l=1,isp(j)
!       do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
        do i=1,ii!![TW 28.1.2014]
          sao(i,j)=-(fuh(i+1,j)-fuh(i,j)+fvh(i,j+1)-fvh(i,j))*scp2i(i,j)
        enddo
!       enddo
!       enddo
      enddo
!$OMP END PARALLEL DO
!
      return
      end subroutine weno3pd_v2
