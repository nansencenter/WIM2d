      !p_save_grid.F
      !Author: Timothy Williams
      !Date:   20141128, 16:06:18 CET

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      program p_save_grid

      use mod_save_grid
      !use mod_save_grid_vars
      implicit none

      character(len=80) :: outdir,outdir2
      character(len=80) :: infile
      integer,parameter :: testing  = 0

      !! grid prams: set here or read in from infile_grid.txt
      integer  :: ii,jj,GRID_OPT
      real     :: reals(4)
      logical  :: fex

      integer :: nfreq,ndir
      real    :: Tmin,Tmax

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
      infile   = 'infile_grid.txt'
      inquire(file=trim(infile),exist=fex) 

      if (fex) then
         !! read from infile
         print*,' '
         print*,'Reading ',infile
         print*,' '
         call read_infile_grid(infile,ii,jj,GRID_OPT,reals)
      else
         !! else set here
         ii       = 150
         jj       = 20
         GRID_OPT = 0
         !!
         reals(1) = 4.e3   !dx
         reals(2) = 4.e3   !dy
         reals(3) = 0.     !x0
         reals(4) = 0.     !y0
      end if

      print*,'nx            : ',ii
      print*,'ny            : ',jj
      print*,'dx       (km) : ',reals(1)/1.e3
      print*,'dy       (km) : ',reals(2)/1.e3
      print*,'x0       (km) : ',reals(3)/1.e3
      print*,'y0       (km) : ',reals(4)/1.e3
      print*,'x extent (km) : ',ii*reals(1)/1.e3
      print*,'y extent (km) : ',jj*reals(2)/1.e3
      print*,' '

      call set_grid(outdir,ii,jj,GRID_OPT,reals)

      !! save main parameters to grid_info.h:
      call save_grid_info_hdr(outdir2,ii,jj,reals(1),reals(2))

      !! save main parameters to wave_info.h:
      infile   = 'infile_waves.txt'
      inquire(file=trim(infile),exist=fex) 
      if (fex) then
         !! read from infile
         print*,'Reading ',infile
         print*,' '
         call read_infile_waves(infile,nfreq,ndir,Tmin,Tmax)
      else
         nfreq = 1
         ndir  = 16
         Tmin  = 2.5
         Tmax  = 25.
      end if

      print*,'nfreq    : ',nfreq
      print*,'ndir     : ',ndir
      print*,'Tmin (s) : ',Tmin
      print*,'Tmax (s) : ',Tmax
      print*,' '

      call save_wave_info_hdr(outdir2,nfreq,ndir,Tmin,Tmax)

      end program p_save_grid
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine read_infile_grid(infile,ii,jj,GRID_OPT,reals)
      !use mod_save_grid_vars
      implicit none

      character(len=*),intent(in)   :: infile
      integer,intent(out)           :: ii,jj,GRID_OPT
      real,intent(out)              :: reals(4)
      character(len=80)             :: ctmp
      integer,parameter             :: vers_no0 = 1
      integer                       :: fid,vers_no

 111  format(i1.1,a)
 114  format(i4.4,a)  
 212  format(f9.2,a)  

      open(unit=fid,file=infile,action='read')
      read(fid,111), vers_no,ctmp

      if (vers_no.ne.vers_no0) then
         print*,'Wrong version number for '//trim(infile)
         print*,'Version number : ',vers_no
         print*,'Should be      : ',vers_no0
         close(fid)
         stop
      end if

      read(fid,114)  ,ii         ,ctmp !nx
      read(fid,114)  ,jj         ,ctmp !ny
      read(fid,212)  ,reals(1)   ,ctmp !dx
      read(fid,212)  ,reals(2)   ,ctmp !dy
      read(fid,212)  ,reals(3)   ,ctmp !x0
      read(fid,212)  ,reals(4)   ,ctmp !y0
      read(fid,111)  ,GRID_OPT   ,ctmp

      close(fid)
      end subroutine read_infile_grid
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine read_infile_waves(infile,nfreq,ndir,Tmin,Tmax)
      !use mod_save_grid_vars
      implicit none

      character(len=*),intent(in)   :: infile
      integer,intent(out)           :: nfreq,ndir
      real,intent(out)              :: Tmin,Tmax
      character(len=80)             :: ctmp
      integer,parameter             :: vers_no0 = 1
      integer                       :: fid,vers_no

 111  format(i1.1,a)
 112  format(i2.2,a)  
 113  format(i3.3,a)  
 212  format(f5.2,a)  

      open(unit=fid,file=infile,action='read')
      read(fid,111), vers_no,ctmp

      if (vers_no.ne.vers_no0) then
         print*,'Wrong version number for '//trim(infile)
         print*,'Version number : ',vers_no
         print*,'Should be      : ',vers_no0
         close(fid)
         stop
      end if

      read(fid,112)  ,nfreq      ,ctmp !nx
      read(fid,113)  ,ndir       ,ctmp !ny
      read(fid,212)  ,Tmin       ,ctmp !dx
      read(fid,212)  ,Tmax       ,ctmp !dy

      close(fid)
      end subroutine read_infile_waves
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine set_grid(outdir,ii,jj,GRID_OPT,reals)
      use mod_save_grid

      implicit none

      integer,intent(in)            :: ii,jj,GRID_OPT
      real,intent(in)               :: reals(4)
      character(len=*),intent(in)   :: outdir

      integer                 :: i,j
      real                    :: x0,y0,dx,dy
      real                    :: xx(ii),yy(jj),xm,ym,Rsq
      real                    :: x1,y1,R0
      real,dimension(ii,jj)   :: X,Y,scuy,scvx,scp2,scp2i,LANDMASK
      real,dimension(ii,jj,7) :: grid_arrays

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      dx = reals(1)
      dy = reals(2)
      x0 = reals(3)
      y0 = reals(4)

      xm = (ii*dx)/2.0
      ym = (jj*dy)/2.0

      X  = 0.0
      Y  = 0.0
      do i=1,ii
         xx(i)    = x0+dx*(i-.5) !centres of cells
         X(i,:)   = xx(i)
      end do
      do j=1,jj
         yy(j)    = y0+dy*(j-.5) !centres of cells
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
            !!no land
         elseif (GRID_OPT.eq.1) then

            !!column of land to right
            x1 = x0+1.8*xm
            if (X(i,j).gt.x1) then
               LANDMASK(i,j)  = 1
            end if

         elseif (GRID_OPT.eq.2) then

            !!make an island;
            x1    = x0+2*xm-.25*(2*xm) !!centre
            y1    = y0+.25*(2*ym)      !!centre
            R0    = .2*ym              !!radius
            Rsq   = (X(i,j)-x1)**2+(Y(i,j)-y1)**2
            if (Rsq.lt.R0**2) then
               LANDMASK(i,j)  = 1
            end if

         elseif (GRID_OPT.eq.3) then

            !!land in top,bottom and right cells
            LANDMASK(ii,:) = 1   !far right cells
            LANDMASK(:,jj) = 1   !top cells
            LANDMASK(:,1)  = 1   !bottom cells

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

      !! save grid as binary
      call save_grid(outdir,grid_arrays,ii,jj)

      end subroutine set_grid
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
