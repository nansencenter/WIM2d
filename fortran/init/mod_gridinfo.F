      !mod_gridinfo.F
      !Author: Timothy Williams
      !Date:   20141128, 16:06:18 CET

      module mod_gridinfo

      implicit none

      include "grid_info.h"   !ii,jj,dx,dy (grid dimensions and resolution)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!global variables
      integer, parameter   :: nbdy  = 0
      !!
      integer, parameter   :: idm   = ii
      integer, parameter   :: jdm   = jj
      !!
      real  ,dimension(idm,jdm)   ::                                    &
     &   X,Y,scuy,scvx,scp2,scp2i,LANDMASK
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      contains

      subroutine get_grid(GRID_OPT)

      implicit none

      integer, intent(in)  :: GRID_OPT
      integer              :: i,j
      real                 :: xx(ii),yy(jj),xm,ym,Rsq
      real                 :: x0,y0,R0

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!initialise global variables
      xm = (ii+1)/2.0*dx
      ym = (jj+1)/2.0*dy

      X  = 0.0
      Y  = 0.0
      do i=1,ii
         xx(i)    = -xm+dx*i
         X(i,:)   = xx(i)
      end do
      do j=1,jj
         yy(j)    = -ym+dy*j
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
            x0    = xm-.25*(2*xm)   !!centre
            y0    = -ym+.25*(2*ym)  !!centre
            R0    = 30e3            !!radius
            Rsq   = (X(i,j)-x0)**2+(Y(i,j)-y0)**2
            if (Rsq.lt.R0**2) then
               LANDMASK(i,j)  = 1
            end if

         elseif (GRID_OPT.eq.1) then

            !!column of land to right
            x0 = .5*xm
            if (X(i,j).gt.x0) then
               LANDMASK(i,j)  = 1
            end if

         elseif (GRID_OPT.eq.2) then

            !!land in bottom right corner
            x0 = .5*xm
            y0 = -.5*ym
            if ( (X(i,j).gt.x0).and.(Y(i,j).lt.y0) ) then
               LANDMASK(i,j)  = 1
            end if

         end if

      end do
      end do

      end subroutine get_grid
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_grid_file()

      use mod_file_utils, only: rd_afile_2d

      implicit none

      character(len=80) :: afile

      afile = 'inputs/wim_grid.a'
      call rd_afile_2d(afile ,X          ,ii,jj,1)
      call rd_afile_2d(afile ,Y          ,ii,jj,2)
      call rd_afile_2d(afile ,scuy       ,ii,jj,3)
      call rd_afile_2d(afile ,scvx       ,ii,jj,4)
      call rd_afile_2d(afile ,scp2       ,ii,jj,5)
      call rd_afile_2d(afile ,scp2i      ,ii,jj,6)
      call rd_afile_2d(afile ,LANDMASK   ,ii,jj,7)

      end subroutine get_grid_file
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module mod_gridinfo
