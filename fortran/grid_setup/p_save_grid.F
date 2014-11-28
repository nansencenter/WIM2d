      !mod_gridinfo.F
      !Author: Timothy Williams
      !Date:   20141128, 16:06:18 CET

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      program p_save_grid

      implicit none

      character(len=80) :: outdir,outdir2

      if (.false.) then
         !!for testing
         outdir   = 'test'
         outdir2  = 'test'
      else
         !!proper locations of outputs:
         outdir   = '../run/inputs'    !.a,.b
         outdir2  = '../header_files'  !.h
      end if

      !! save main parameters to grid_info.h:
      call save_grid_info_hdr(outdir2)

      !! main subroutine
      !! - define the grid
      call set_grid()

      !! save grid as binary
      call save_grid(outdir)

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
      real  ,dimension(ii,jj)   ::                                      &
     &   X,Y,scuy,scvx,scp2,scp2i,LANDMASK

      end module mod_grid
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine save_grid_info_hdr(outdir)

      use mod_grid, only: ii,jj,dx,dy

      implicit none
      character(len=80),intent(in)  :: outdir
      character(len=*),parameter    :: blk   = '      '
      !!
      character(len=80) :: hfile,ss
      character(len=10) :: ssi,ssj

      !! write .h file
      hfile = trim(outdir)//'/grid_info.h'
      print*,'Making header file : ',hfile

      open(unit=2,file=trim(hfile),status = 'replace')

      ss = blk//"!! header file with grid dimensions and resolution"
      write(2,'(a)') trim(ss)
      ss = blk//"!! - file generated automatically with setup_grid.sh"  &
     &        //"/p_save_grid.F"
      write(2,'(a)') trim(ss)
      ss = blk//"!! - do not edit !!"
      write(2,'(a)') trim(ss)
      write(2,'(a)') ''

      write(ssi,'(i3.3)'),ii
      write(ssj,'(i3.3)'),jj
      write(2,'(a,a,a)')                                                &
     &  blk,ssi,   "   !! Record length in x direction (elements)"
      write(2,'(a,a,a)')                                                &
     &  blk,ssj,   "   !! Record length in y direction (elements)"
      write(2,'(a,e10.4,a)')                                            &
     &  blk,dx,   "   !! Resolution in x direction"
      write(2,'(a,e10.4,a)')                                            &
     &  blk,dy,   "   !! Resolution in y direction"

      close(2)

      end subroutine save_grid_info_hdr
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine set_grid()
      use mod_grid

      implicit none

      integer              :: i,j
      real                 :: xx(ii),yy(jj),xm,ym,Rsq
      real                 :: x0,y0,R0

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

      end subroutine set_grid
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine save_grid(outdir)

      use mod_grid
      use mod_file_utils

      implicit none
      character(len=80),intent(in)   :: outdir
      character(len=80)    :: afile,bfile
      character(len=10)    :: stat1,stat2

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
     &  5,          "       Number of records"
      write(2,112)                                                      &
     &  ii,          "      Record length in x direction (elements)"
      write(2,112)                                                      &
     &  jj,          "      Record length in y direction (elements)"
      write(2,111)                                                      &
     &  GRID_OPT,   "       Option number for initial conditions"
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
      call wrt_afile_2d(afile ,X          ,ii,jj  ,stat1,1)
      call wrt_afile_2d(afile ,Y          ,ii,jj  ,stat2,2)
      call wrt_afile_2d(afile ,scuy       ,ii,jj  ,stat2,3)
      call wrt_afile_2d(afile ,scvx       ,ii,jj  ,stat2,4)
      call wrt_afile_2d(afile ,scp2       ,ii,jj  ,stat2,5)
      call wrt_afile_2d(afile ,scp2i      ,ii,jj  ,stat2,6)
      call wrt_afile_2d(afile ,LANDMASK   ,ii,jj  ,stat2,7)

      end subroutine save_grid
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
