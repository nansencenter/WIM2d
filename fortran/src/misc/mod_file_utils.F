      module mod_file_utils

      contains

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine wrt_afile_1d(outfile,x,nx,stat)

      implicit none

      character(len=*), intent(in)    :: outfile
      real, dimension(nx), intent(in)  :: x  !!1d array to write in
      integer, intent(in)              :: nx !!length of x
      character(len=*), intent(in)    :: stat
      !! stat can be:
      !! 'replace' - overwrite file
      !! 'new'     - file doesn't exist already
      !! 'old'     - append


      integer  :: rec_len

      inquire(IOLENGTH=rec_len) x !should be 4*nx

      open(unit=1,file=trim(outfile),status = trim(stat),               &
     &   form='unformatted', access='direct',recl=rec_len)
      write(1,rec=1) x
      close(1)

      end subroutine wrt_afile_1d
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine wrt_afile_2d(outfile,x,nx,ny,stat,recno)

      implicit none

      character(len=*), intent(in)       :: outfile
      real, dimension(nx,ny), intent(in)  :: x     !!2d array to write in
      integer, intent(in)                 :: nx,ny !!size of x
      integer, intent(in)                 :: recno !!record no
      character(len=*), intent(in)       :: stat
      !! stat can be:
      !! 'replace' - overwrite file
      !! 'new'     - file doesn't exist already
      !! 'old'     - append

      integer  :: rec_len

      inquire(IOLENGTH=rec_len) x   !should be 4*nx*ny

      open(unit=1,file=trim(outfile),status = trim(stat),               &
     &   form='unformatted', access='direct',recl=rec_len)
      write(1,rec=recno) x
      close(1)

      end subroutine wrt_afile_2d
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine wrt_afile_3d(outfile,x,nx,ny,nz,stat,recno)

      implicit none

      character(len=*), intent(in)       :: outfile
      real, dimension(nx,ny), intent(in)  :: x        !!3d array to write in
      integer, intent(in)                 :: nx,ny,nz !!size of x
      integer, intent(in)                 :: recno    !!record no
      character(len=*), intent(in)       :: stat
      !! stat can be:
      !! 'replace' - overwrite file
      !! 'new'     - file doesn't exist already
      !! 'old'     - append

      integer  :: rec_len

      inquire(IOLENGTH=rec_len) x   !should be 4*nx*ny*nz

      open(unit=1,file=trim(outfile),status = trim(stat),               &
     &   form='unformatted', access='direct',recl=rec_len)
      write(1,rec=recno) x
      close(1)

      end subroutine wrt_afile_3d
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine rd_afile_1d(infile,x,nx,recno)

      implicit none

      character(len=*), intent(in)       :: infile
      real, dimension(nx), intent(out)    :: x     !!1d array to read in
      integer, intent(in)                 :: nx    !!size of x
      integer, intent(in)                 :: recno !!record no
      integer                             :: rec_len

      x  = 0.0
      inquire(IOLENGTH=rec_len) x   !should be 4*nx

      open(unit=1,file=trim(infile),status = 'old',                     &
     &   form='unformatted', access='direct',recl=rec_len)
      read(1,rec=recno) x
      close(1)

      end subroutine rd_afile_1d
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine rd_afile_2d(infile,x,nx,ny,recno)

      implicit none

      character(len=*), intent(in)        :: infile
      real, dimension(nx,ny), intent(out) :: x     !!2d array to read in
      integer, intent(in)                 :: nx,ny !!size of x
      integer, intent(in)                 :: recno !!record no
      real*4, dimension(nx,ny)            :: x4 !if wim_grid.[a,b] are single precision arrays
      real*8, dimension(nx,ny)            :: x8 !if wim_grid.[a,b] are double precision arrays

      character(len=200)   :: bfile,ctmp
      integer              :: fid,rec_len,filesize,filename_len,nrecs
      logical              :: fex

      filename_len   = len(trim(infile))
      bfile          =  infile(1:filename_len-2)//'.b'

      x     = 0.0
      x4    = 0.0
      fid   = 99
      inquire(IOLENGTH=rec_len) x4   !should be 4*nx*ny
      inquire(file=infile,exist=fex)

      if (fex) then
         inquire(file=infile,SIZE=filesize)

         ! get number of records from bfile
         open(unit=fid,file=trim(bfile),action = 'read')
         read(fid,'(i2.2,a)') nrecs,ctmp
         close(fid)

         print*,recno,rec_len,nrecs
         print*,filesize

         if (filesize.eq.rec_len*nrecs) then
            !file is in single precision
            open(unit=fid,file=trim(infile),status = 'old',             &
     &         form='unformatted', access='direct',recl=rec_len)
            read(fid,rec=recno) x4
            close(fid)
            x  = real(x4)
         else
            !file is in double precision
            open(unit=fid,file=trim(infile),status = 'old',             &
     &         form='unformatted', access='direct',recl=2*rec_len)
            read(fid,rec=recno) x8
            close(fid)
            x  = real(x8)
         end if

      else
         print*,'file not present: ',infile
         stop
      end if

      end subroutine rd_afile_2d
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine rd_afile_3d(infile,x,nx,ny,nz,recno)

      implicit none

      character(len=*), intent(in)       :: infile
      real, dimension(nx,ny), intent(out) :: x        !!3d array to read in
      integer, intent(in)                 :: nx,ny,nz !!size of x
      integer, intent(in)                 :: recno    !!record no
      integer                             :: rec_len

      inquire(IOLENGTH=rec_len) x   !should be 4*nx*ny*nz

      open(unit=1,file=trim(infile),status = 'old',                     &
     &   form='unformatted', access='direct',recl=rec_len)
      read(1,rec=recno) x
      close(1)

      end subroutine rd_afile_3d
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module mod_file_utils
