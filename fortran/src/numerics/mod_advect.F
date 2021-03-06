      module mod_advect

      implicit none

      private

      !!public subroutines:
      public   :: waveadv_weno
     &           ,waveadv_weno_1d
     &           ,waveadv_lax_1d

      contains

      ! ---------------------------------------------------------------
      ! Wrappers (called from outside)
      ! ---------------------------------------------------------------

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine waveadv_weno(h,u,v,dt,ADV_OPT)
      !! WENO 2d wrapper
!
! --- ------------------------------------------------------------------
! --- Advection is done with flux limited 3rd order WENO in space and
! --- 2nd order Runge-Kutta in time
! --- ------------------------------------------------------------------
!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! modified [TW 28.1.2014] for wave advection
      !! - iceadv builds up ice (& thus waves) at coast;
      !! - iceadv_v2 lets waves go into/out of land,
      !!   then cancels the waves on land at end of routine;
      !! - just adds the wave cancellation to Luke's original fix
      !!   for waves on boundaries;
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! ARGUMENTS:
      !! in/out:  h is qty to be advected;
      !!
      !! in:      u,v,dt are speed components and time step;
      !!
      !! Now in mod_gridinfo.F
      !! in: scuy,scvx is mesh size at u (v) points in y (x) direction
      !!      - see common_blocks.h;
      !! in: scp2, scp2i are grid box area at p points, and its inverse;
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use mod_gridinfo
      use mod_file_utils
      use mod_common_advect, only: pad_var

      implicit none

      integer, parameter   :: n_bdy = 4
         !n_bdy=6: don't need to apply boundary condition between
         ! prediction/correction steps (WENO advection)
         !n_bdy=3: do need to apply boundary condition between
         ! prediction/correction steps (WENO advection)

      real,intent(in),dimension(idm,jdm) :: u,v
      !,scuy,scvx,scp2i,scp2,LANDMASK
      real,intent(inout),dimension(idm,jdm)  :: h
      real,intent(in)                        :: dt
      integer,intent(in)                     :: ADV_OPT
!
      real, dimension(1-n_bdy:idm+n_bdy,1-n_bdy:jdm+n_bdy) :: sao,hp    &
     &   ,u_pad,v_pad,scp2_pad,scp2i_pad,scuy_pad,scvx_pad,h_pad
      integer  :: i,j,l
      character(len=80) :: outfile,stat1,stat2

      sao   = 0.0

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!pad out variables
      !! - apply boundary conditions via "ghost cells"

      !!make u,v,sc* periodic
      call pad_var(u_pad      ,u       ,1,n_bdy)
      call pad_var(v_pad      ,v       ,1,n_bdy)
      call pad_var(scp2_pad   ,scp2    ,1,n_bdy)
      call pad_var(scp2i_pad  ,scp2i   ,1,n_bdy)
      call pad_var(scuy_pad   ,scuy    ,1,n_bdy)
      call pad_var(scvx_pad   ,scvx    ,1,n_bdy)

      !!ADV_OPT applies boundary conditions
      !!0: make h just be zero on boundaries
      !!1: make h be periodic in x,y
      !!2: make h be periodic in y only
      call pad_var(h_pad,h,ADV_OPT,n_bdy)

      if (.false.) then
         !! confirm in correct subroutine
         print*,'******************************************'
         print*,'In waveadv_weno (2d advection)'
         print*,'nx,ny = ',ii,jj
         print*,'******************************************'
      end if

      if (.false.) then
      !if (.true.) then
         !!test padding with zero on boundaries
         stat1    = 'replace'
         outfile  = 'test_out/h.a'

         print*,'******************************************'
         print*,'testing waveadv_weno:'
         print*,'check',outfile
         print*,'******************************************'
         print*,''

         i        = idm+2*n_bdy
         j        = jdm+2*n_bdy
         call wrt_afile_2d(outfile,h_pad,i,j,stat1,1)
      elseif (.false.) then
         !!test padding with periodic boundary conditions
         stat1    = 'replace'
         stat2    = 'old'
         outfile  = 'test_out/all.a'

         print*,'******************************************'
         print*,'testing waveadv_weno:'
         print*,'check: ',outfile
         print*,'******************************************'
         print*,''

         i     = idm+2*n_bdy
         j     = jdm+2*n_bdy
         call wrt_afile_2d(outfile,u_pad,    i,j,stat1,1)
         call wrt_afile_2d(outfile,v_pad,    i,j,stat2,2)
         call wrt_afile_2d(outfile,scp2_pad, i,j,stat2,3)
         call wrt_afile_2d(outfile,scp2i_pad,i,j,stat2,4)
         call wrt_afile_2d(outfile,scuy_pad, i,j,stat2,5)
         call wrt_afile_2d(outfile,scvx_pad, i,j,stat2,6)
      end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
! --- Prediction step
      !print*,'prediction'
      call weno3pd_v2(h_pad,sao,u_pad,v_pad,                            &
     &   scuy_pad,scvx_pad,scp2i_pad,scp2_pad,dt,n_bdy)

      if (n_bdy.gt.3) then
         !! original version
         do j=1-n_bdy,jj+n_bdy
         do i=1-n_bdy,ii+n_bdy!![TW 28.1.2014]
             hp(i,j) = h_pad(i,j)+dt*sao(i,j)
         enddo
         enddo
      elseif (n_bdy.eq.3) then
         !! apply boundary condition between prediction and correction step
         do j=1,jj
         do i=1,ii!![TW 28.1.2014]
             hp(i,j) = h_pad(i,j)+dt*sao(i,j)
         enddo
         enddo
         call pad_var(hp,hp(1:ii,1:jj),ADV_OPT,n_bdy)
      else
         print*,'Choose n_bdy>=3'
         stop
      end if

! --- Correction step
      !print*,'correction'
      call weno3pd_v2(hp,sao,u_pad,v_pad,                               &
     &   scuy_pad,scvx_pad,scp2i_pad,scp2_pad,dt,n_bdy)

      !! final output
      do j=1,jj
      do i=1,ii!![TW 28.1.2014]
           h(i,j) = .5*(h_pad(i,j)+hp(i,j)+dt*sao(i,j))
      enddo
      enddo

      !!mask land (no waves on land)
      h  = h*(1-LANDMASK)

      end subroutine waveadv_weno
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine waveadv_weno_1d(h,u,dt,ADV_OPT)
      !! WENO 1d wrapper

!
! --- ------------------------------------------------------------------
! --- Advection is done with flux limited 3rd order WENO in space and
! --- 2nd order Runge-Kutta in time
! --- ------------------------------------------------------------------
!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! modified [TW 28.1.2014] for wave advection
      !! - iceadv builds up ice (& thus waves) at coast;
      !! - iceadv_v2 lets waves go into/out of land,
      !!   then cancels the waves on land at end of routine;
      !! - just adds the wave cancellation to Luke's original fix
      !!   for waves on boundaries;
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! ARGUMENTS:
      !! in/out:  h is qty to be advected;
      !!
      !! in:      u,v,dt are speed components and time step;
      !!
      !! Now in mod_gridinfo.F
      !! in: scuy,scvx is mesh size at u (v) points in y (x) direction
      !!      - see common_blocks.h;
      !! in: scp2, scp2i are grid box area at p points, and its inverse;
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use mod_gridinfo
      use mod_file_utils
      use mod_common_advect, only: pad_var_1d

      implicit none

      integer, parameter   :: n_bdy = 4
         !n_bdy=6: don't need to apply boundary condition between
         ! prediction/correction steps (WENO advection)
         !n_bdy=3: do need to apply boundary condition between
         ! prediction/correction steps (WENO advection)

      real,intent(in),dimension(idm,jdm)     :: u
      real,intent(inout),dimension(idm,jdm)  :: h
      real,intent(in)                        :: dt
      integer,intent(in)                     :: ADV_OPT
!
      real, dimension(1-n_bdy:idm+n_bdy,1:jdm) :: sao,hp                &
     &   ,u_pad,scp2_pad,scp2i_pad,scuy_pad,h_pad
      integer  :: i,j,l
      character(len=80) :: outfile,stat1,stat2

      sao   = 0.0

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!pad out variables
      !! - apply boundary conditions via "ghost cells"

      !!make u,v,sc* periodic
      call pad_var_1d(u_pad    ,u    ,1,n_bdy)
      call pad_var_1d(scp2_pad ,scp2 ,1,n_bdy)
      call pad_var_1d(scp2i_pad,scp2i,1,n_bdy)
      call pad_var_1d(scuy_pad ,scuy ,1,n_bdy)

      !!ADV_OPT applies boundary conditions
      !!0: make h just be zero on boundaries
      !!1: make h be periodic in x
      call pad_var_1d(h_pad,h,ADV_OPT,n_bdy)

      if (.false.) then
         !! confirm in correct subroutine
         print*,'******************************************'
         print*,'In waveadv_weno_1d'
         print*,'******************************************'
      endif

      if (.false.) then
      !if (.true.) then
         !!test padding of h
         stat1    = 'replace'
         outfile  = 'test_out/h.a'

         print*,'******************************************'
         print*,'testing waveadv_weno:'
         print*,'check',outfile

         i  = idm+2*n_bdy
         j  = jdm
         call wrt_afile_2d(outfile,h_pad,i,j,stat1,1)

         print*,'******************************************'
         print*,''
      end if

      if (.false.) then
      !if (.true.) then
         !!test padding of other things
         !!(periodic bc's)
         stat1    = 'replace'
         stat2    = 'old'
         outfile  = 'test_out/all.a'

         print*,'******************************************'
         print*,'testing waveadv_weno:'
         print*,'check: ',outfile

         i  = idm+2*n_bdy
         j  = jdm
         call wrt_afile_2d(outfile,u_pad,    i,j,stat1,1)
         call wrt_afile_2d(outfile,scp2_pad, i,j,stat2,2)
         call wrt_afile_2d(outfile,scp2i_pad,i,j,stat2,3)
         call wrt_afile_2d(outfile,scuy_pad, i,j,stat2,4)

         print*,'******************************************'
         print*,''
      end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
! --- Prediction step
      call weno3pd_v2_1d(h_pad,sao,u_pad,                               &
     &   scuy_pad,scp2i_pad,scp2_pad,dt,n_bdy)

      if (n_bdy.gt.3) then
         !! original version - but I think this needed xctilr also
         do j=1,jj
         do i=1-n_bdy,ii+n_bdy
             hp(i,j) = h_pad(i,j)+dt*sao(i,j)
         enddo
         enddo
      elseif (n_bdy.eq.3) then
         !! apply boundary condition between prediction and correction step
         do j=1,jj
         do i=1,ii
             hp(i,j) = h_pad(i,j)+dt*sao(i,j)
         enddo
         enddo
         call pad_var_1d(hp,hp(1:ii,1:jj),ADV_OPT,n_bdy)
      else
         print*,'choose n_bdy>=3'
         stop
      end if

! --- Correction step
      call weno3pd_v2_1d(hp,sao,u_pad,                                  &
     &   scuy_pad,scp2i_pad,scp2_pad,dt,n_bdy)

      !final output
      do j=1,jj
      do i=1,ii
           h(i,j)  = .5*(h_pad(i,j)+hp(i,j)+dt*sao(i,j))
      enddo
      enddo

      !!mask land (no waves on land)
      h  = h*(1-LANDMASK)

      end subroutine waveadv_weno_1d
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine waveadv_lax_1d(h,u,dt,ADV_OPT)
      !! Lax-Wendroff 1d wrapper
!
! --- ------------------------------------------------------------------
! --- Advection is done with flux limited 3rd order WENO in space and
! --- 2nd order Runge-Kutta in time
! --- ------------------------------------------------------------------
!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! modified [TW 28.1.2014] for wave advection
      !! - iceadv builds up ice (& thus waves) at coast;
      !! - iceadv_v2 lets waves go into/out of land,
      !!   then cancels the waves on land at end of routine;
      !! - just adds the wave cancellation to Luke's original fix
      !!   for waves on boundaries;
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! ARGUMENTS:
      !! in/out:  h is qty to be advected;
      !!
      !! in:      u,v,dt are speed components and time step;
      !!
      !! Now in mod_gridinfo.F
      !! in: scuy,scvx is mesh size at u (v) points in y (x) direction
      !!      - see common_blocks.h;
      !! in: scp2, scp2i are grid box area at p points, and its inverse;
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use mod_gridinfo
      use mod_file_utils
      use mod_common_advect, only: pad_var_1d

      implicit none

      integer, parameter   :: n_bdy = 4!TODO this could probably be reduced for this adv scheme

      real,intent(in),dimension(idm,jdm) :: u
      !,scuy,scvx,scp2i,scp2,LANDMASK
      real,intent(inout),dimension(idm,jdm)  :: h
      real,intent(in)                        :: dt
      integer,intent(in)                     :: ADV_OPT
!
      real, dimension(1-n_bdy:idm+n_bdy,1:jdm) :: sao,hp                &
     &   ,u_pad,scp2_pad,scp2i_pad,scuy_pad,h_pad
      integer  :: i,j,l
      character(len=80) :: outfile,stat1,stat2

      sao   = 0.0

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!pad out variables
      !! - apply boundary conditions via "ghost cells"

      !!make u,v,sc* periodic
      call pad_var_1d(u_pad    ,u    ,1,n_bdy)
      call pad_var_1d(scp2_pad ,scp2 ,1,n_bdy)
      call pad_var_1d(scp2i_pad,scp2i,1,n_bdy)
      call pad_var_1d(scuy_pad ,scuy ,1,n_bdy)

      !!ADV_OPT applies boundary conditions
      !!0: make h just be zero on boundaries
      !!1: make h be periodic in x
      call pad_var_1d(h_pad,h,ADV_OPT,n_bdy)

      if (.false.) then
      !if (.true.) then
         !!test padding of h
         stat1    = 'replace'
         outfile  = 'test_out/h.a'

         print*,'******************************************'
         print*,'testing waveadv_weno:'
         print*,'check',outfile

         i  = idm+2*n_bdy
         j  = jdm
         call wrt_afile_2d(outfile,h_pad,i,j,stat1,1)

         print*,'******************************************'
         print*,''
      end if

      if (.false.) then
      !if (.true.) then
         !!test padding of other things
         !!(periodic bc's)
         stat1    = 'replace'
         stat2    = 'old'
         outfile  = 'test_out/all.a'

         print*,'******************************************'
         print*,'testing waveadv_weno:'
         print*,'check: ',outfile

         i  = idm+2*n_bdy
         j  = jdm
         call wrt_afile_2d(outfile,u_pad,    i,j,stat1,1)
         call wrt_afile_2d(outfile,scp2_pad, i,j,stat2,2)
         call wrt_afile_2d(outfile,scp2i_pad,i,j,stat2,3)
         call wrt_afile_2d(outfile,scuy_pad, i,j,stat2,4)

         print*,'******************************************'
         print*,''
      end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
! --- Single call to get Lax-Wendroff spatial advection operator
      call sao_lax_1d(h_pad,sao,u_pad,                                  &
     &   scuy_pad,scp2i_pad,scp2_pad,dt,n_bdy)

      do j=1,jj
      do i=1,ii
          h(i,j)  = h_pad(i,j)+dt*sao(i,j)
      enddo
      enddo

      !!mask land (no waves on land)
      h  = h*(1-LANDMASK)

      end subroutine waveadv_lax_1d
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! ---------------------------------------------------------------
      ! End of wrappers (called from outside)
      ! ---------------------------------------------------------------


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! ---------------------------------------------------------------
      ! End of boundary conditions
      ! ---------------------------------------------------------------
      

      ! ---------------------------------------------------------------
      ! Main algorithms
      ! ---------------------------------------------------------------

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine weno3pd_v2(g,sao,u,v,scuy,scvx,scp2i,scp2,dt,n_bdy)
!
! --- ------------------------------------------------------------------
! --- By a weighted essentially non-oscillatory scheme with up to 3rd
! --- order accuracy (2/3 in time/space), obtain the spatial advective operator of a
! --- 2-dimensional field defined at the scalar points of a C-grid. The
! --- fluxes are limited to make the scheme positive definite.
! --- Advective velocities in the i- and j-direction are defined at u-
! --- and v-points, respectively.
! --- ------------------------------------------------------------------
!

      use mod_gridinfo, only: ii,jj,idm,jdm

      implicit none

      integer,intent(in)   :: n_bdy

      real,intent(in),                                                  &
     &  dimension(1-n_bdy:idm+n_bdy,1-n_bdy:jdm+n_bdy) :: g,u,v         &
     &    ,scp2,scp2i,scuy,scvx
      real,intent(out),                                                 &
     &  dimension(1-n_bdy:idm+n_bdy,1-n_bdy:jdm+n_bdy) :: sao
      real,intent(in) :: dt
!
      real   ::  cq00,cq01,cq10,cq11,ca0,ca1,eps
      parameter (cq00=-1./2.,cq01= 3./2.,                               &
     &           cq10= 1./2.,cq11= 1./2.,                               &
     &           ca0=1./3.,ca1=2./3.,                                   &
     &           eps=1.e-12)
!
      real, dimension(1-n_bdy:idm+n_bdy,1-n_bdy:jdm+n_bdy) ::           &
     &  ful,fuh,fvl,fvh,gt
      real   :: q0,q1,a0,a1,q
      integer  :: i,j,l,im1,im2,ip1,jm1,jm2,jp1
!
! --- Compute grid cell boundary fluxes. Split in a low order flux
! --- (donor cell) and a high order correction flux.
!
      ful  = 0.
      fuh  = 0.
      fvl  = 0.
      fvh  = 0.

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!fluxes in x dirn:
      do j=1-n_bdy,jj+n_bdy
      do i=1-n_bdy+2,ii+n_bdy-1
          im1=i-1
!
          if (u(i,j).gt.0.) then
            im2=im1-1
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
            ip1=i+1
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
      enddo
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!fluxes in y dirn:
      do j=1-n_bdy+2,jj+n_bdy-1
      do i=1-n_bdy,ii+n_bdy

          jm1=j-1
          if (v(i,j).gt.0.) then
            jm2=jm1-1
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
            jp1=j+1
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
      enddo
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! --- Update field with low order fluxes.
      do j = 1-n_bdy,jj+n_bdy-1
      do i = 1-n_bdy,ii+n_bdy-1
          gt(i,j)=g(i,j)-dt*(ful(i+1,j)-ful(i,j)                        &
     &                      +fvl(i,j+1)-fvl(i,j))*scp2i(i,j)
      enddo
      enddo

! --- Obtain fluxes with limited high order correction fluxes.
      q=.25/dt
      do j = 1-n_bdy,jj+n_bdy
      do i = 1-n_bdy+1,ii+n_bdy
          fuh(i,j)=ful(i,j)+max(-q*gt(i  ,j)*scp2(i  ,j),               &
     &                      min( q*gt(i-1,j)*scp2(i-1,j),fuh(i,j)))
      enddo
      enddo

      do j = 1-n_bdy+1,jj+n_bdy
      do i = 1-n_bdy,ii+n_bdy
          fvh(i,j)=fvl(i,j)+max(-q*gt(i,j  )*scp2(i,j  ),               &
     &                      min( q*gt(i,j-1)*scp2(i,j-1),fvh(i,j)))
      enddo
      enddo

! --- Compute the spatial advective operator.
      do j = 1-n_bdy,jj+n_bdy-1
      do i = 1-n_bdy,ii+n_bdy-1
          sao(i,j)=-(fuh(i+1,j)-fuh(i,j)+fvh(i,j+1)-fvh(i,j))*scp2i(i,j)
      enddo
      enddo

      end subroutine weno3pd_v2
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine weno3pd_v2_1d(g,sao,u,scuy,scp2i,scp2,dt,n_bdy)
!
! --- ------------------------------------------------------------------
! --- By a weighted essentially non-oscillatory scheme with up to 3rd
! --- order accuracy (2/3 in time/space), obtain the spatial advective operator of a
! --- 2-dimensional field defined at the scalar points of a C-grid. The
! --- fluxes are limited to make the scheme positive definite.
! --- Advective velocities in the i- and j-direction are defined at u-
! --- and v-points, respectively.
! --- ------------------------------------------------------------------
!

      use mod_gridinfo, only: ii,jj,idm,jdm

      implicit none

      integer,intent(in)   :: n_bdy

      real,intent(in),                                                  &
     &  dimension(1-n_bdy:idm+n_bdy,1:jdm) :: g,u                       &
     &    ,scp2,scp2i,scuy
      real,intent(out),                                                 &
     &  dimension(1-n_bdy:idm+n_bdy,1:jdm) :: sao
      real,intent(in) :: dt
!
      real   ::  cq00,cq01,cq10,cq11,ca0,ca1,eps
      parameter (cq00=-1./2.,cq01= 3./2.,                               &
     &           cq10= 1./2.,cq11= 1./2.,                               &
     &           ca0=1./3.,ca1=2./3.,                                   &
     &           eps=1.e-12)
!
      real, dimension(1-n_bdy:idm+n_bdy,1:jdm) ::                       &
     &  ful,fuh,gt
      real   :: q0,q1,a0,a1,q
      integer  :: i,j,l,im1,im2,ip1,jm1,jm2,jp1
!
! --- Compute grid cell boundary fluxes. Split in a low order flux
! --- (donor cell) and a high order correction flux.
!
      ful  = 0.
      fuh  = 0.

      !!fluxes in x dirn:
      do j=1,jj
      do i=1-n_bdy+2,ii+n_bdy-1
         im1=i-1
!
         if (u(i,j).gt.0.) then
            !! lower-order fluxes
            im2      = im1-1
            ful(i,j) = u(i,j)*g(im1,j)*scuy(i,j)
!
            !! coefficents to calc higher-order fluxes
            q0 = cq00*g(im2,j)+cq01*g(im1,j)
            q1 = cq10*g(im1,j)+cq11*g(i  ,j)
            a0 = ca0
            a1 = ca1*(abs(g(im2,j)-g(im1,j))+eps)                       &
     &            /(abs(g(im1,j)-g(i  ,j))+eps)
!
          else
            !! lower-order fluxes
            ip1      = i+1
            ful(i,j) = u(i,j)*g(i,j)*scuy(i,j)

            !! coefficents to calc higher-order fluxes
            q0 = cq11*g(im1,j)+cq10*g(i  ,j)
            q1 = cq01*g(i  ,j)+cq00*g(ip1,j)
            a0 = ca1
            a1 = ca0*(abs(g(im1,j)-g(i  ,j))+eps)                       &
     &            /(abs(g(i  ,j)-g(ip1,j))+eps)
          endif

          !! higher-order fluxes
          fuh(i,j)=u(i,j)*(a0*q0+a1*q1)/(a0+a1)*scuy(i,j)-ful(i,j)
!
        enddo
      enddo

! --- Update field with low order fluxes.
      do j=1,jj
      do i = 1-n_bdy,ii+n_bdy-1
         gt(i,j)  = g(i,j)-dt*(ful(i+1,j)-ful(i,j))*scp2i(i,j)
      enddo
      enddo

! --- Obtain fluxes with limited high order correction fluxes.
      q=.25/dt
      do j=1,jj
      do i = 1-n_bdy+1,ii+n_bdy
          fuh(i,j)=ful(i,j)+max(-q*gt(i  ,j)*scp2(i  ,j),               &
     &                      min( q*gt(i-1,j)*scp2(i-1,j),fuh(i,j)))
      enddo
      enddo

! --- Compute the spatial advective operator.
      do j=1,jj
      do i = 1-n_bdy,ii+n_bdy-1
          sao(i,j)=-(fuh(i+1,j)-fuh(i,j))*scp2i(i,j)
      enddo
      enddo

      end subroutine weno3pd_v2_1d
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine sao_lax_1d(g,sao,u,scuy,scp2i,scp2,dt,n_bdy)
!
! --- ------------------------------------------------------------------
! --- By a weighted essentially non-oscillatory scheme with up to 3rd
! --- order accuracy (2/3 in time/space), obtain the spatial advective operator of a
! --- 2-dimensional field defined at the scalar points of a C-grid. The
! --- fluxes are limited to make the scheme positive definite.
! --- Advective velocities in the i- and j-direction are defined at u-
! --- and v-points, respectively.
! --- ------------------------------------------------------------------
!

      use mod_gridinfo, only: ii,jj,idm,jdm

      implicit none

      !flux limiter method
      integer,parameter :: FL_method   = 2

      integer, intent(in)  :: n_bdy

      real,intent(in),                                                  &
     &  dimension(1-n_bdy:idm+n_bdy,1:jdm) :: g,u                       &
     &    ,scp2,scp2i,scuy
      real,intent(out),                                                 &
     &  dimension(1-n_bdy:idm+n_bdy,1:jdm) :: sao
      real,intent(in) :: dt
!
!
      real, dimension(1-n_bdy:idm+n_bdy,1:jdm) :: flux_ip1
      real     :: q_flux,r_flux,dx_i,dx_ip1,up,um,ful,fuh
      integer  :: i,j,im1,im2,ip1,ip2
!
! --- Compute grid cell boundary fluxes. Split in a low order flux
! --- (donor cell) and a high order correction flux.
!
      do j=1,jj
      do i=0,ii+2
         flux_ip1(i,j)  = 0.0 !!F_{i+1/2,j}
      enddo
      enddo

      !!fluxes in x dirn:


      do j=1,jj
      do i=0,ii+2!![TW 28.1.2014]

         ip1   = i+1
         ip2   = i+2
         im1   = i-1
         im2   = i-2
         !!
         up    = ( u(i,j)   + abs(u(i,j))   )/2.0
         um    = ( u(ip1,j) - abs(u(ip1,j)) )/2.0
         !!
         dx_i     = scuy(i,j)*scp2i(i,j)     !!TODO: check this is x(i+1,j)-x(i,j)   ??
         dx_ip1   = scuy(ip1,j)*scp2i(ip1,j) !!TODO: check this is x(i+2,j)-x(i+1,j) ??

         !! low-order flux
         !! - 1st-order upwind
         !! - NB *scuy because *scp2i later
         ful   = up*g(i,j)+um*g(ip1,j)

         !! 2nd-order flux
         !! - Lax_Wendroff
         !! - NB *scuy because *scp2i later
         fuh   = ful + 0.5*abs(u(ip1,j))*(1-dt*dx_ip1)*g(ip1,j)         &
     &               - 0.5*abs(u(i,j))*(1-dt*dx_i)*g(i,j)

         if (u(i,j).gt.0.0) then
            r_flux   = (g(i,j)-g(im1,j))/(g(ip1,j)-g(i,j))
         else
            r_flux   = (g(ip2,j)-g(ip1,j))/(g(ip1,j)-g(i,j))
         end if
         !!
         call flux_limiter(q_flux,r_flux,FL_method)
         flux_ip1(i,j)  = ful+q_flux*(fuh-ful)
      enddo
      enddo

! --- Compute the spatial advective operator.
      do j=1,jj
      do i=1,ii
         dx_i     = scuy(i,j)*scp2i(i,j)  !!1/dx
         sao(i,j) = -(flux_ip1(i,j)-flux_ip1(i-1,j))*dx_i
      enddo
      enddo

      return
      end subroutine sao_lax_1d
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! ---------------------------------------------------------------
      ! End of main algorithms
      ! ---------------------------------------------------------------


      ! ---------------------------------------------------------------
      ! Additional code
      ! ---------------------------------------------------------------

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine flux_limiter(r,phi,FL_method)
      implicit none
      integer, intent(in)  :: FL_method
      real, intent(in)  :: r
      real, intent(out) :: phi
      real  :: beta,phi0,phi1

      if (FL_method.eq.1) then
         phi   = (abs(r)+r)/(abs(r)+1.0)                 !!van Leer
      elseif (FL_method.eq.2) then
         phi   = max(0.0,max(min(1.0,2*r),min(r,2.0)))   !!Superbee
      elseif (FL_method.eq.3) then
         phi = (r**2+r)/(r**2+1.0)                       !! van Albada 1
      elseif (FL_method.eq.4) then
         beta  = 1.5 !! choose/tune 1<=beta<=2;
         phi0  = min(beta*r,1.0)
         phi1  = min(r,beta)
         phi   = max(0.0,max(phi0,phi1))                 !! Sweby
      elseif (FL_method.eq.5) then
         phi = 1.5*(r**2+r)/(r**2+r+1.0)                 !! Ospre
      elseif (FL_method.eq.6) then
         beta  = 1.5 !! choose/tune 1<=beta<=2;
         phi   = max(0.0,min(r,beta))                    !! Osher
      elseif (FL_method.eq.7) then
         phi = max(0.0,min(1.0,r))                       !! minmod
      elseif (FL_method.eq.8) then
         phi0  = min(2*r,(1+r)/2.0)
         phi   = max(0.0,min(phi0,2.0))                  !! MC (monotonised central)

      !!Not symmetric:
      elseif (FL_method.eq.9) then
         phi0  = min(2*r,(1+2*r)/3.0)
         phi   = max(0.0,min(phi0,2.0))                  !! Koren
         !! 3rd order accurate for smooth-enough data
      elseif (FL_method.eq.10) then
         phi0  = min(2*r,(1+3*r)/4.0)
         phi1  = min(2.0,(3+r)/4.0)
         phi   = max(0.0,min(phi0,phi1))                 !! UMIST

      !! non 2nd order TVD:   
      elseif (FL_method.eq.11) then
         phi = (2*r)/(r**2+1.0)                          !! van Albada 2
      elseif (FL_method.eq.12) then
         phi0  = min(2*r,(1+3*r)/4.0)
         phi   = max(0.0,min(phi0,4.0))                  !! smart
      elseif (FL_method.eq.13) then
         phi = 2*(r+abs(r))/(r+3.0)                      !! HQUICK
      elseif (FL_method.eq.14) then
         phi = 1.5*(r+abs(r))/(r+2.0)                    !! HCUS
      elseif (FL_method.eq.15) then
         phi      = 0.0
         if (r.gt.0.0) then
            phi   = r*(3*r+1.0)/(r+1.0)**2               !! CHARM
         end if
      end if

      end subroutine flux_limiter
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module mod_advect
