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

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_grid_from_file(indir)
      !! read in the grid from file

      use mod_file_utils, only: rd_afile_2d

      implicit none

      character(len=*),intent(in)   :: indir
      character(len=80)             :: afile
      logical                       :: fex

      afile = trim(indir)//'/wim_grid.a'
      inquire(file=afile,exist=fex)
      if (.not.fex) then
         print*,'grid file = '//afile
         stop('grid file not present')
      else
         call rd_afile_2d(afile ,X          ,ii,jj,1)
         call rd_afile_2d(afile ,Y          ,ii,jj,2)
         call rd_afile_2d(afile ,scuy       ,ii,jj,3)
         call rd_afile_2d(afile ,scvx       ,ii,jj,4)
         call rd_afile_2d(afile ,scp2       ,ii,jj,5)
         call rd_afile_2d(afile ,scp2i      ,ii,jj,6)
         call rd_afile_2d(afile ,LANDMASK   ,ii,jj,7)
      end if

      end subroutine get_grid_from_file
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_grid_auto()
      !! make grid automatically

      implicit none

      integer           :: i,j

      do j=1,jj
      do i=1,ii
         X(i,j)      = i*dx ! x increases along i
         Y(i,j)      = j*dy ! y increases along j
      end do
      end do

      scuy     = dy
      scvx     = dx
      scp2     = dx*dy
      scp2i    = 1./(dx*dy)
      LANDMASK = 0.

      end subroutine get_grid_auto
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module mod_gridinfo
