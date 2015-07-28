      !p_save_grid.F
      !Author: Timothy Williams
      !Date:   20141128, 16:06:18 CET

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      program p_save_grid

      use mod_save_grid, only: save_grid_info_hdr,save_grid
      use mod_grid
      implicit none

      character(len=80)       :: outdir,outdir2
      integer,parameter       :: testing  = 1

      if (testing.eq.1) then
         !!for testing
         outdir   = 'test/out'
         outdir2  = 'test/out'
      else
         !!proper locations of outputs:
         outdir   = '../run/inputs'    !.a,.b
         outdir2  = '../header_files'  !.h
      end if

      !! main subroutine
      !! - define the grid
      call set_grid()

      !! save main parameters to grid_info.h:
      call save_grid_info_hdr(outdir2,ii,jj,dx,dy)

      !! save grid as binary
      call save_grid(outdir,grid_arrays,ii,jj)

      end program p_save_grid
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module mod_grid
      !!change ii,jj,dx,dy,GRID_OPT here
      !!remaining stuff set in set_grid

      implicit none

      integer, parameter   :: GRID_OPT = 1
      integer,parameter    :: ii       = 150
      integer,parameter    :: jj       = 50
      real,parameter       :: dx       = 4.0e3!m
      real,parameter       :: dy       = 4.0e3!m
      integer, parameter   :: nbdy     = 0
      integer, parameter   :: idm      = ii
      integer, parameter   :: jdm      = jj
      !!
      real  ,dimension(ii,jj) ::                                        &
     &   X,Y,scuy,scvx,scp2,scp2i,LANDMASK
      !!
      real,dimension(ii,jj,7) :: grid_arrays

      end module mod_grid
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine set_grid()
      use mod_grid

      implicit none

      integer              :: i,j
      real                 :: xx(ii),yy(jj),xm,ym,Rsq
      real                 :: x0,y0,x1,y1,R0

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      x0 = 0.
      y0 = 0.
      xm = (ii*dx)/2.0
      ym = (jj*dy)/2.0

      X  = 0.0
      Y  = 0.0
      do i=1,ii
         xx(i)    = dx*(i-1)
         X(i,:)   = xx(i)
      end do
      do j=1,jj
         yy(j)    = dy*(j-1)
         Y(:,j)   = yy(j)
      end do
      !!
      LANDMASK = 0.0
      scuy     = dy
      scvx     = dx
      scp2     = scuy*scvx
      do j=1,jj
      do i=1,ii
         scp2i(i,j)  = 1/scp2(i,j)

         !!set LANDMASK
         if (GRID_OPT.eq.0) then

            !!make an island;
            x1    = x0+2*xm-.25*(2*xm) !!centre
            y1    = y0+.25*(2*ym)      !!centre
            R0    = 30.e3              !!radius
            Rsq   = (X(i,j)-x1)**2+(Y(i,j)-y1)**2
            if (Rsq.lt.R0**2) then
               LANDMASK(i,j)  = 1
            end if

         elseif (GRID_OPT.eq.1) then

            !!column of land to right
            x1 = x0+1.5*xm
            if (X(i,j).gt.x1) then
               LANDMASK(i,j)  = 1
            end if

         elseif (GRID_OPT.eq.2) then

            !!land in bottom right corner
            x1 = x0+1.5*xm
            y1 = y0+0.5*ym
            if ( (X(i,j).gt.x1).and.(Y(i,j).lt.y1) ) then
               LANDMASK(i,j)  = 1
            end if

         end if

      end do
      end do

      grid_arrays(:,:,1)   = X
      grid_arrays(:,:,2)   = Y
      grid_arrays(:,:,3)   = scuy
      grid_arrays(:,:,4)   = scvx
      grid_arrays(:,:,5)   = scp2
      grid_arrays(:,:,6)   = scp2i
      grid_arrays(:,:,7)   = LANDMASK

      end subroutine set_grid
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
