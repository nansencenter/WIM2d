      module mod_advect_ice

      use mod_common_advect

      implicit none

      private

      !!public subroutines:
      public   :: init_iceadv
     &           ,iceadv_remap

      contains

      ! ---------------------------------------------------------------
      ! Wrappers (called from outside)
      ! ---------------------------------------------------------------

c ==============================================================================
      subroutine init_iceadv()

      use mod_gridinfo

      implicit none

      integer :: i,j

      ! ===================================================================
      ! make masks according to bigrid.F from HYCOM

      ip = 1-nint(LANDMASK) ! 1 if water at the p-point
      iu = 0
      iv = 0

      do j=2,jj
      do i=2,ii
         if (ip(i-1,j).gt.0.and.ip(i,j).gt.0) then
           iu(i,j)=1
         endif
         if (ip(i,j-1).gt.0.and.ip(i,j).gt.0) then
           iv(i,j)=1
         endif
      enddo
      enddo

      call indxi(ip,ifp,ilp,isp)
      call indxi(iu,ifu,ilu,isu)
      call indxi(iv,ifv,ilv,isv)
      ! ===================================================================


      end subroutine init_iceadv
c ==============================================================================


c ==============================================================================
      subroutine iceadv_remap(dp,u,v,dt,ADV_OPT)
      !! wrapper for incremental remapping advection scheme
!
! --- ------------------------------------------------------------------
! --- Advection is done with flux limited 3rd order WENO in space and
! --- 2nd order Runge-Kutta in time
! --- ------------------------------------------------------------------
!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! ARGUMENTS:
      !! in/out:  dp is qty to be advected;
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
      use mod_remap
!
      implicit none
!
      integer, parameter   :: n_bdy = 2 !no of ghost cells

      real,intent(in),dimension(idm,jdm) :: u,v
      !,scuy,scvx,scp2i,scp2,LANDMASK
      real,intent(inout),dimension(idm,jdm)  :: dp
      real,intent(in)                        :: dt
      integer,intent(in)                     :: ADV_OPT
!
      real, dimension(1-n_bdy:idm+n_bdy,1-n_bdy:jdm+n_bdy) ::           &
     &    u_pad,v_pad,scp2_pad,scp2i_pad,scuy_pad,scvx_pad,dp_pad       &
     &   ,pbmin,plo,pbu,pbv,temp_pad,saln_pad                           &
     &   ,uflx,vflx,utflx,vtflx,usflx,vsflx

      integer  :: i,j,l
      character(len=80) :: outfile,stat1,stat2

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
      !!0: make dp just be zero on boundaries
      !!1: make dp be periodic in x,y
      !!2: make dp be periodic in y only
      call pad_var(dp_pad,dp,ADV_OPT,n_bdy)

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
         call wrt_afile_2d(outfile,dp_pad,i,j,stat1,1)
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

      ! inputs for ocean bottom topography (TODO remove these)
      pbmin = 1000.
      pbu   = 1000.
      pbv   = 1000.

      ! pressure at lower interface of ocean layer
      ! - make plo << pbmin,pbu,pbv so topography plays no role
      plo   = 900.
      
      !temp and salinity
      temp_pad = 10.
      saln_pad = 30.

      call remap_eitvel(scuy_pad,scvx_pad,scp2i_pad,scp2_pad
     &                        ,pbmin,pbu,pbv,plo
     &                        ,u,v,dt
     &                        ,n_bdy
     &                        ,dp_pad,temp_pad,saln_pad
     &                        ,uflx,vflx,utflx,vtflx,usflx,vsflx
#ifdef TRC
     &                       ,k,trc
#endif
     &                       )

      !! final output
      do j=1,jj
      do i=1,ii
           dp(i,j) = dp_pad(i,j)
      enddo
      enddo

      end subroutine iceadv_remap
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module mod_advect_ice
