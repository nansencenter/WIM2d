      module mod_common_advect

      use mod_gridinfo

      implicit none

      integer,parameter :: ms = 100 !max no of sections

      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     &   ip,iu,iv !wet masks for p,u,v points

      integer, dimension (1-nbdy:jdm+nbdy) ::
     &   isp,isu,isv !no of wet sections in each column

      integer, dimension (1-nbdy:jdm+nbdy,ms) ::
     &   ifp,ifu,ifv
     &  ,ilp,ilu,ilv
           ! ifp(j,k) is first point of the k-th section in column j (according to the ip mask)
           ! ilp(j,k) is last  point of the k-th section in column j (according to the ip mask)

      contains


      ! ---------------------------------------------------------------
      ! Boundary conditions
      ! ---------------------------------------------------------------

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine pad_var(u_pad,u,OPT,n_bdy)

      use mod_gridinfo, only: ii,jj,idm,jdm

      implicit none

      integer,intent(in)   :: OPT,n_bdy
      real,intent(out),                                                 &
     &  dimension(1-n_bdy:idm+n_bdy,1-n_bdy:jdm+n_bdy) :: u_pad
      real,intent(in),dimension(idm,jdm)               :: u

      u_pad             = 0.0
      u_pad(1:ii,1:jj)  = u

      if (OPT.eq.1) then
         !!make things periodic

         !!make it periodic in i
         u_pad(1-n_bdy:0,1:jj)       = u(ii-n_bdy:ii-1,:)
         u_pad(ii+1:ii+n_bdy,1:jj)   = u(1:n_bdy,:)

         !!make it periodic in j
         u_pad(1:ii,1-n_bdy:0)       = u(1:ii,jj-n_bdy:jj-1)
         u_pad(1:ii,jj+1:jj+n_bdy)   = u(1:ii,1:n_bdy)

         !!BR,TL
         u_pad(ii+1:ii+n_bdy,jj+1:jj+n_bdy)                             &
     &     = u(1:n_bdy,1:n_bdy)
         u_pad(1-n_bdy:0,1-n_bdy:0)                                     &
     &     = u(ii-n_bdy:ii-1,jj-n_bdy:jj-1)

         !!BL,TR
         u_pad(ii+1:ii+n_bdy,1-n_bdy:0)                                 &
     &     = u(1:n_bdy,jj-n_bdy:jj-1)
         u_pad(1-n_bdy:0,jj+1:jj+n_bdy)                                 &
     &     = u(ii-n_bdy:ii-1,1:n_bdy)

      else if (OPT.eq.2) then
         !!make it periodic in j only
         u_pad(1:ii,1-n_bdy:0)       = u(1:ii,jj-n_bdy:jj-1)
         u_pad(1:ii,jj+1:jj+n_bdy)   = u(1:ii,1:n_bdy)

      end if

      end subroutine pad_var
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine pad_var_1d(u_pad,u,OPT,n_bdy)

      use mod_gridinfo, only: ii,jj,idm,jdm

      implicit none

      integer,intent(in)   :: OPT,n_bdy
      real,intent(out),                                                 &
     &  dimension(1-n_bdy:idm+n_bdy,1:jdm)   :: u_pad
      real,intent(in),dimension(idm,jdm)     :: u

      u_pad             = 0.0
      u_pad(1:ii,1:jj)  = u

      if (OPT.eq.1) then
         !!make it periodic in i
         u_pad(1-n_bdy:0,1:jj)       = u(ii-n_bdy:ii-1,:)
         u_pad(ii+1:ii+n_bdy,1:jj)   = u(1:n_bdy,:)
      end if

      end subroutine pad_var_1d
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine indxi(ipt,if,il,is)
      implicit none
c
      integer, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     &         ipt
      integer, dimension (1-nbdy:jdm+nbdy,ms) ::
     &         if,il
      integer, dimension (1-nbdy:jdm+nbdy) ::
     &         is
c
c --- input array ipt contains 1 at grid point locations, 0 elsewhere
c --- output is arrays if, il, is  where
c --- if(j,k) gives row index of first point in column j for k-th section
c --- il(j,k) gives row index of last point
c --- is(j) gives number of sections in column j (maximum: ms)
c
      integer i,j,k,last
c
      do j=1-nbdy,jj+nbdy
        is(j) = 0
        do k=1,ms
          if(j,k) = 0
          il(j,k) = 0
        end do
c
        k=1
        last = ipt(1-nbdy,j)
        if     (last .eq. 1) then
          if(j,k) = 1-nbdy
        endif
        do i=2-nbdy,ii+nbdy
          if      (last .eq. 1 .and. ipt(i,j) .eq. 0) then
            il(j,k) = i-1
            k = k+1
          elseif (last .eq. 0 .and. ipt(i,j) .eq. 1) then
            if     (k .gt. ms) then
              print*,
     &          ' error in indxi -- ms too small at i,j =',i,j
              stop '(indxi)'
            endif
            if(j,k) = i
          endif
          last = ipt(i,j)
        enddo
        if     (last .eq. 1) then
          il(j,k) = ii+nbdy
          is(j) = k
        else
          is(j) = k-1
        endif
      enddo

      return
      end subroutine indxi
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module mod_common_advect
