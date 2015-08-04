      !save_grid_f2py.F
      !Author: Timothy Williams
      !Date:   20141128, 16:06:18 CET

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine save_grid_info_hdr_f2py(outdir,ii,jj,dx,dy)

      use mod_save_grid, only: save_grid_info_hdr
      implicit none

      character(len=*),intent(in)   :: outdir
      integer, intent(in)           :: ii,jj
      real, intent(in)              :: dx,dy

      call save_grid_info_hdr(outdir,ii,jj,dx,dy)

      end subroutine save_grid_info_hdr_f2py
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine save_grid_f2py(outdir,grid_arrays,ii,jj)

      use mod_save_grid, only: save_grid

      implicit none
      integer,parameter                   :: nz = 7
      integer,intent(in)                  :: ii,jj
      character(len=*),intent(in)         :: outdir
      real,dimension(ii,jj,nz),intent(in) :: grid_arrays

      !comments for f2py
      !don't need to pass these inputs in with python
!f2py intent(hide)   :: ii,jj

      call save_grid(outdir,grid_arrays,ii,jj)

      end subroutine save_grid_f2py
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
