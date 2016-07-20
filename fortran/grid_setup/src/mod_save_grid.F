      !save_grid_f2py.F
      !Author: Timothy Williams
      !Date:   20141128, 16:06:18 CET

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     module mod_save_grid_vars
!     implicit none

!     integer  :: GRID_OPT
!     integer  :: ii,jj
!     real     :: dx,dy,x0,y0
!     !!
!     real  ,dimension(ii,jj) ::                                        &
!    &   X,Y,scuy,scvx,scp2,scp2i,LANDMASK
!     !!
!     real,dimension(ii,jj,7) :: grid_arrays

!     end module mod_save_grid_vars
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module mod_save_grid

      implicit none

      contains
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine save_wave_info_hdr(outdir,nfreq,ndir,Tmin,Tmax)

      implicit none

      character(len=*),intent(in)   :: outdir
      integer,intent(in)            :: nfreq,ndir
      real,intent(in)               :: Tmin,Tmax
      character(len=*),parameter    :: blk   = '      '
      !!
      character(len=80) :: hfile,ss
      character(len=80) :: ssi,ssj,ssx,ssy

      !! write .h file
      hfile = trim(outdir)//'/wave_info.h'
      print*,'Making header file : ',hfile
      print*,' '

      open(unit=2,file=trim(hfile),status = 'replace')

 310  format(a)
 311  format(a,i4.4)
 312  format(a,i4.4,a)
 411  format(a,f6.2)

      write(2,311),  blk//'integer,parameter   :: n_wavdir    = ',ndir
      write(2,312),  blk//'integer,parameter   :: n_wave_freq = ',nfreq &
     &                  ," !should be odd (Simpson's rule)"
      write(2,310),  ''
      write(2,411), blk//'real,parameter   :: Tmin = ',Tmin
      write(2,411), blk//'real,parameter   :: Tmax = ',Tmax

      close(2)

      end subroutine save_wave_info_hdr
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine save_grid_info_hdr(outdir,ii,jj,dx,dy)

      implicit none

      character(len=*),intent(in)   :: outdir
      integer,intent(in)            :: ii,jj
      real,intent(in)               :: dx,dy
      character(len=*),parameter    :: blk   = '      '
      !!
      character(len=80) :: hfile,ss
      character(len=80) :: ssi,ssj,ssx,ssy

      !! write .h file
      hfile = trim(outdir)//'/grid_info.h'
      print*,'Making header file : ',hfile
      print*,' '

      open(unit=2,file=trim(hfile),status = 'replace')

      ss = blk//"!! header file with grid dimensions and resolution"
      write(2,'(a)') trim(ss)
      ss = blk//"!! - file generated automatically in "//               &
     &          "fortran/grid_setup"
      write(2,'(a)') trim(ss)
      ss = blk//"!! - do not edit !!"
      write(2,'(a)') trim(ss)
      write(2,'(a)') ''

      write(ssi,'(a,i4.4)'),  'integer,parameter   :: ii = ',ii
      write(ssj,'(a,i4.4)'),  'integer,parameter   :: jj = ',jj
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
      character(len=*),intent(in)         :: outdir
      real,dimension(ii,jj,nz),intent(in) :: grid_arrays
      !!
      character(len=80)    :: afile,bfile
      integer              :: j

      !comments for f2py
      !don't need to pass these inputs in with python
!f2py intent(hide)   :: ii,jj

      afile = trim(outdir)//'/wim_grid.a'
      bfile = trim(outdir)//'/wim_grid.b'

      !! write .b file (description of stuff in .a file)
      print*,'Making .b file     : ',bfile
 111  format(a)
 112  format(i2.2,a)  
      call wrt_bfile_header(bfile,nz,ii,jj)
      open(unit=2,file=trim(bfile),position='append')

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(2,111)    "Record number and name:"
      write(2,112)    1,    "       X"
      write(2,112)    2,    "       Y"
      write(2,112)    3,    "       scuy"
      write(2,112)    4,    "       scvx"
      write(2,112)    5,    "       scp2"
      write(2,112)    6,    "       scp2i"
      write(2,112)    7,    "       LANDMASK"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      close(2)

      !! dump arrays to .a file
      print*,'Making .a file     : ',afile

      !! write 7 2d arrays
      j  = 1
      call wrt_afile_2d(afile,grid_arrays(:,:,j),ii,jj,'replace',j)
      do j=2,nz
         call wrt_afile_2d(afile,grid_arrays(:,:,j),ii,jj,'old',j)
      end do

      end subroutine save_grid
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine wrt_bfile_header(bfile,Nrecs,ii,jj)
      character(len=*),intent(in)   :: bfile
      integer,intent(in)            :: Nrecs,ii,jj

 111  format(a)
 112  format(i2.2,a)  
 113  format(i3.3,a)  
 114  format(i4.4,a)  

      !! write .b file (description of stuff in .a file)
      open(unit=2,file=trim(bfile),status = 'replace')

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(2,112) Nrecs,  "       Nrecs    # Number of records"
      write(2,112) 1,      "       Norder   # Storage order "//         &
     &               "[column-major (F/matlab) = 1, row-major (C) = 0]"
      write(2,114) ii,     "     nx       "//                           &
     &               "# Record length in x direction (elements)"
      write(2,114) jj,     "     ny       "//                           &
     &               "# Record length in y direction (elements)"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      write(2,111) ""
      close(2)

      end subroutine wrt_bfile_header
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module mod_save_grid
