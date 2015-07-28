      !save_grid_f2py.F
      !Author: Timothy Williams
      !Date:   20141128, 16:06:18 CET

      module mod_save_grid
      contains
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine save_grid_info_hdr(outdir,ii,jj,dx,dy)

      implicit none

      integer, intent(in)           :: ii,jj
      real, intent(in)              :: dx,dy
      character(len=*),intent(in)  :: outdir
      character(len=*),parameter    :: blk   = '      '
      !!
      character(len=80) :: hfile,ss
      character(len=80) :: ssi,ssj,ssx,ssy

      !! write .h file
      hfile = trim(outdir)//'/grid_info.h'
      print*,'Making header file : ',hfile

      open(unit=2,file=trim(hfile),status = 'replace')

      ss = blk//"!! header file with grid dimensions and resolution"
      write(2,'(a)') trim(ss)
      ss = blk//"!! - file generated automatically in "//               &
     &          "fortran/grid_setup"
      write(2,'(a)') trim(ss)
      ss = blk//"!! - do not edit !!"
      write(2,'(a)') trim(ss)
      write(2,'(a)') ''

      write(ssi,'(a,i3.3)'),  'integer,parameter   :: ii = ',ii
      write(ssj,'(a,i3.3)'),  'integer,parameter   :: jj = ',jj
      write(ssx,'(a,e10.4)'), 'real,parameter   :: dx = ',dx
      write(ssy,'(a,e10.4)'), 'real,parameter   :: dy = ',dy

      write(2,'(a,a,a)')                                                &
     &  blk,trim(ssi),   "   !! Record length in x direction (elements)"
      write(2,'(a,a,a)')                                                &
     &  blk,trim(ssj),   "   !! Record length in y direction (elements)"
      write(2,'(a)') " "
      write(2,'(a,a,a)')                                                &
     &  blk,trim(ssx),   "   !! Resolution in x direction"
      write(2,'(a,a,a)')                                                &
     &  blk,trim(ssy),   "   !! Resolution in y direction"

      close(2)

      end subroutine save_grid_info_hdr
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine save_grid(outdir,grid_arrays,ii,jj)

      use mod_file_utils

      implicit none
      integer,parameter                   :: nz = 7
      integer,intent(in)                  :: ii,jj
      character(len=*),intent(in)        :: outdir
      real,dimension(ii,jj,nz),intent(in) :: grid_arrays
      character(len=80)    :: afile,bfile
      character(len=10)    :: stat1,stat2
      integer              :: j

      !comments for f2py
      !don't need to pass these inputs in with python
!f2py intent(hide)   :: ii,jj

      afile = trim(outdir)//'/wim_grid.a'
      bfile = trim(outdir)//'/wim_grid.b'
      stat1 = 'replace'
      stat2 = 'old'

 111  format(i2.2,a)  
 112  format(i3.3,a)  
 113  format(a)

      !! write .b file (description of stuff in .a file)
      print*,'Making .b file     : ',bfile
      open(unit=2,file=trim(bfile),status = 'replace')

      write(2,111)                                                      &
     &  7,          "       Number of records"
      write(2,111)                                                      &
     &  1,          "       Storage order "//                           &
     &              "[column-major (F/matlab) = 1, row-major (C) = 0]"
      write(2,112)                                                      &
     &  ii,          "      Record length in x direction (elements)"
      write(2,112)                                                      &
     &  jj,          "      Record length in y direction (elements)"
      write(2,113) ""

      write(2,113)    "Record number and name:"
      write(2,111)                                                      &
     &  1,    "       X"
      write(2,111)                                                      &
     &  2,    "       Y"
      write(2,111)                                                      &
     &  3,    "       scuy"
      write(2,111)                                                      &
     &  4,    "       scvx"
      write(2,111)                                                      &
     &  5,    "       scp2"
      write(2,111)                                                      &
     &  6,    "       scp2i"
      write(2,111)                                                      &
     &  7,    "       LANDMASK"

      close(2)

      !! dump arrays to .a file
      !! - X,Y,scuy,scvx,scp2,scp2i,LANDMASK
      print*,'Making .a file     : ',afile

      if (.false.) then
         !! write 3d array directly
         call wrt_afile_3d(afile,grid_arrays,ii,jj,nz,stat1,1)
      else
         !! write 7 2d arrays
         j  = 1
         call wrt_afile_2d(afile,grid_arrays(:,:,j),ii,jj,stat1,j)
         do j=2,7
            call wrt_afile_2d(afile,grid_arrays(:,:,j),ii,jj,stat2,j)
         end do
      end if

      end subroutine save_grid
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module mod_save_grid
