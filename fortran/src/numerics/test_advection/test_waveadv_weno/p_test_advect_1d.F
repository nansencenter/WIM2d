      !! test_waveadv_weno.F
      !! Author: Timothy Williams
      !! Date:   20140821, 12:22:17 CEST
      program main

      use mod_advect
      use mod_gridinfo
      use mod_file_utils

      implicit none

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!variables
      integer              :: OPT,ADV_opt
      integer, parameter   :: ADV_METHOD = 1
      real                 :: CFL
      !!
      real   :: x0,y0,xm,ym
      real   :: xx(ii),yy(jj)
      real,dimension(idm,jdm)   :: R,Theta,u,h,v,source

      !! can try to have attenuation (set to 0.0 if don't want it)
      !real,parameter :: alp   = 4.0e-4 !needs alp*dt<<1
      real,parameter :: alp   = 0.0

      real   :: adv_dir,xc,yc,uc,Rc,angrot,dt,dtheta,PI
      integer  :: nt,n,i,j
      logical  :: critter,fex2

      character(len=*),parameter :: outfile = "ADVweno"
      character(len=*),parameter :: outdir = "out/"
      character(len=3)           :: cts
      character(len=80)          :: afile,bfile,stat1,stat2,ss
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      PI = 4.0*atan(1.0)

      inquire(file='infile.txt',exist=fex2)
      if (fex2) then
         call read_infile(OPT,ADV_opt,CFL,adv_dir)
      else
         OPT      = 1 ! geometry
         ADV_OPT  = 2
            !1: periodic in x,y; 2: periodic in y only; 0: not periodic
         CFL      = .4
         adv_dir  = 135.0!deg
      end if

      v  = 0.0 !!place holder for test functions

      source   = 0.0

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!set grid here
      x0 = 0.
      y0 = 0.
      !!
      X  = 0.0
      Y  = 0.0
      do i=1,ii
         xx(i)    = x0+dx*(i-1)
         X(i,:)   = xx(i)
      end do
      xm = .5*(maxval(X)-x0)
      !!
      do j=1,jj
         yy(j)    = y0+dy*(j-1)
         Y(:,j)   = yy(j)
      end do
      ym = .5*(maxval(Y)-y0)
      !!
      LANDMASK = 0.0
      scuy     = dy
      scvx     = dx
      scp2     = scuy*scvx
      do j=1,jj
      do i=1,ii
         scp2i(i,j)  = 1.0/scp2(i,j)
      end do
      end do
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do j=1,jj
      do i=1,ii
         R(i,j)      = sqrt(X(i,j)**2+Y(i,j)**2)
         Theta(i,j)  = atan2(Y(i,j),X(i,j))
      end do
      end do

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!Test situations
      if (OPT.eq.1) then
         !!constant speed
         print*,'******************************************'
         print*,'using constant speed: OPT = ',OPT
         print*,'******************************************'
         print*,''
         uc       = 30!const speed m/s
         xc       = x0+5*xm/3.0
         u        = uc*cos(pi/180.0*adv_dir)
         dt       = CFL*dx/uc
         if (ADV_opt.eq.1) then
            nt = 2*aint(2*xm/(uc*dt))  !test periodic boundary condition
         elseif (ADV_opt.eq.2) then
            nt = aint(3*xm/(uc*dt))  !test periodic boundary condition
         else
            nt = aint(2*xm/(uc*dt))
         end if
         !!
         h  = 0.0
         do j=1,jj
         do i=1,ii
            if (X(i,j).gt.xc) then
               h(i,j)   = 1.0
            end if
         end do
         end do

      elseif (OPT.eq.2) then
         !!variable speed
         uc = 30!max speed m/s
         xc = x0+5*xm/3.0
         u  = -uc*(X-xm)/xm!speed variable with x
         dt = CFL*dx/uc
         nt = 2*xm/(uc*dt)
         !!
         h  = 0.0
         do j=1,jj
         do i=1,ii
            if (X(i,j).gt.xc) then
               h(i,j)   = 1.0
            end if
         end do
         end do
      end if

      !! Print info
      print*,' '
      print*,'*****************************************'
      print*,'Grid size : ',ii,jj
      if (ADV_METHOD.eq.1) then
         print*,'Scheme    :   WENO'
      else
         print*,'Scheme    :   Lax-Wendroff'
      end if
      print*,'CFL       : ',CFL
      print*,'Time step : ',dt
      print*,'*****************************************'
      print*,' '


      !if (.false.) then
      if (.true.) then
         stat1 = 'replace'
         stat2 = 'old'
         afile = 'test_out/all0.a'

         print*,'******************************************'
         print*,'testing initial conditions:'
         print*,'check: ',afile
         print*,'******************************************'
         print*,''

         call wrt_afile_2d(afile,u,    idm,jdm,stat1,1)
         call wrt_afile_2d(afile,v,    idm,jdm,stat2,2)
         call wrt_afile_2d(afile,scp2, idm,jdm,stat2,3)
         call wrt_afile_2d(afile,scp2i,idm,jdm,stat2,4)
         call wrt_afile_2d(afile,scuy, idm,jdm,stat2,5)
         call wrt_afile_2d(afile,scvx, idm,jdm,stat2,6)
      end if

      !!output initial conditions for reading by matlab/python
      n  = 0
      stat1 = 'replace'
      stat2 = 'old'
      !!
      write(cts,'(i3.3)') n
      afile = outdir//outfile//cts//'.a'
      bfile = outdir//outfile//cts//'.b'
      call wrt_bfile(bfile,ii,jj,OPT)
      call wrt_afile_2d(afile,x,ii,jj,stat1,1)
      call wrt_afile_2d(afile,y,ii,jj,stat2,2)
      call wrt_afile_2d(afile,h,ii,jj,stat2,3)
      !!
      print*,'******************************************'
      print*,'finished writing initial conditions to:'
      print*,afile
      print*,bfile
      print*,'******************************************'
      print*,''

      do n = 1,nt
      !do n = 1,1
         !print*,'timestep = ',n,' out of ',nt

         source   = -alp*h!!attenuation source (sink)
         if (1.eq.0) then
            ss = 'Source before advection'
            h  = h+source*dt
         end if

         if (ADV_METHOD.eq.1) then
            !!advect for one time step
            call waveadv_weno_1d(h,u,dt,ADV_opt)
         else
            !!advect for one time step
            if (.false.) then
               print*,'min/max source:',minval(source),maxval(source)
            endif
            call waveadv_lax_1d(h,u,dt,ADV_opt)
         end if

         if (1.eq.1) then
            if (n.eq.1) ss = '*Split-step attenuation'
            h  = h*exp(-alp*dt)
         elseif (1.eq.1) then
            if (n.eq.1) ss = '*Aprox. split-step attenuation'
            !! gives an idea of effect of time step
            h  = h*(1-alp*dt)
         elseif (1.eq.1) then
            if (n.eq.1) ss = '*Explicit source in advection scheme'
            h  = h+source*dt
         end if

         if (n.eq.1) then
            print*,ss
            print*,'Atten. coeff :',alp
         end if

         !!output initial conditions for reading by matlab/python
         write(cts,'(i3.3)') n
         afile = outdir//outfile//cts//'.a'
         bfile = outdir//outfile//cts//'.b'
         call wrt_bfile(bfile,ii,jj,OPT)
         call wrt_afile_2d(afile,x,ii,jj,stat1,1)
         call wrt_afile_2d(afile,y,ii,jj,stat2,2)
         call wrt_afile_2d(afile,h,ii,jj,stat2,3)
      end do

      print*,' '
      print*,'******************************************'
      print*,'finished writing results to:'
      print*,afile
      print*,bfile
      print*,'******************************************'
      print*,''

      end program main
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine wrt_bfile(outfile,nx,ny,OPT)

      implicit none
      character(len=80), intent(in) :: outfile
      integer, intent(in)           :: nx,ny,OPT
      
 111  format(i2.2,a)  
 112  format(i3.3,a)  
 113  format(a)

      open(unit=2,file=trim(outfile),status = 'replace')

      write(2,111)                                                      &
     &  2,     "       Number of records"
      write(2,112)                                                      &
     &  nx,     "      Record length in x direction (elements)"
      write(2,112)                                                      &
     &  ny,     "      Record length in y direction (elements)"
      write(2,111)                                                      &
     &  OPT,    "      Option number for initial conditions"
      write(2,113) ""
      write(2,113)    "Record number and name:"
      write(2,111)                                                      &
     &  1,    "       X (m)"
      write(2,111)                                                      &
     &  2,    "       Y (m)"
      write(2,111)                                                      &
     &  3,    "       h (advected quantity)"

      close(2)

      end subroutine wrt_bfile
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine read_infile(OPT,ADV_OPT,CFL,adv_dir)
      !!determines where to read inputs/outputs from
      
      implicit none
      
      integer,intent(out)        :: OPT,ADV_OPT
      real,intent(out)           :: CFL,adv_dir
      integer,parameter          :: fid=99
      character(len=80)          :: ctmp
      character(len=*),parameter :: infile = 'infile.txt'

      print*,' '
      print*,'******************************************************'
      print*,'Opening infile: infile.txt'

      open(unit=fid,file=infile,action='read')

      !! integer prams
      read(fid,'(a)'), ctmp
      read(fid,'(i2.2,a)'), OPT     ,ctmp
      print*,'OPT: ',OPT
      read(fid,'(i2.2,a)'), ADV_OPT     ,ctmp
      print*,'ADV_OPT: ',ADV_OPT

      !! real prams
      read(fid,'(a)'), ctmp
      read(fid,'(a)'), ctmp
      read(fid,'(e4.2,a)'),  CFL ,ctmp
      print*,"CFL: ",CFL
      read(fid,'(e4.2,a)'),  adv_dir ,ctmp
      print*,"adv_dir: ",adv_dir
      print*,'******************************************************'
      print*,' '
      !!
      close(fid)

      end subroutine read_infile
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
