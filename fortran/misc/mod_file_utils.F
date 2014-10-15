      !! test_advection_weno.F
      !! Author: Timothy Williams
      !! Date:   20140821, 12:22:17 CEST
      module mod_file_utils

      contains

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine wrt_afile_1d(outfile,x,nx,stat)

      implicit none

      character(len=80), intent(in)       :: outfile
      real*8, dimension(nx), intent(in)   :: x  !!1d array to write in
      integer, intent(in)                 :: nx !!length of x
      character(len=80), intent(in)       :: stat
      !! stat can be:
      !! 'replace' - overwrite file
      !! 'new'     - file doesn't exist already
      !! 'old'     - append


      integer  :: rec_len

      inquire(IOLENGTH=rec_len) x !should be 8*nx

      open(unit=1,file=trim(outfile),status = trim(stat),               &
     &   form='unformatted', access='direct',recl=rec_len)
      write(1,rec=1) x
      close(1)

      end subroutine wrt_afile_1d
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine wrt_afile_2d(outfile,x,nx,ny,stat,recno)

      implicit none

      character(len=80), intent(in)          :: outfile
      real*8, dimension(nx,ny), intent(in)   :: x     !!2d array to write in
      integer, intent(in)                    :: nx,ny !!size of x
      integer, intent(in)                    :: recno !!record no
      character(len=10), intent(in)          :: stat
      !! stat can be:
      !! 'replace' - overwrite file
      !! 'new'     - file doesn't exist already
      !! 'old'     - append

      integer  :: rec_len

      inquire(IOLENGTH=rec_len) x   !should be 8*nx*ny

      open(unit=1,file=trim(outfile),status = trim(stat),               &
     &   form='unformatted', access='direct',recl=rec_len)
      write(1,rec=recno) x
      close(1)

      end subroutine wrt_afile_2d
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module mod_file_utils
