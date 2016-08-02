#include "fintrf.h"
C
#if 0
C
C     CALL:
C     [Sdir_io,out_arrays,mesh_e]=...
C        WIM2d_run_io_mex_vSdir(Sdir_io,ice_arrays,int_prams,...
C                               real_prams,Tp_io,mwd_io,mesh_e);
C
#endif
C
C     This is a MEX-file for MATLAB.

      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
C-----------------------------------------------------------------------
C
      use mod_WIM2d_run, only: mex_io_gate_vSdir_mesh
      implicit none

      include "grid_info.h"
      include "wave_info.h"

      mwpointer plhs(*), prhs(*)
      mwpointer mxCreateDoubleMatrix, mxGetPr
      integer nlhs, nrhs, mxIsNumeric
      mwsize mxGetM, mxGetN
      integer*4   :: mexPrintf,j,k,nmeshi

C-----INPUTS------------------------------------------------------------
      !!pointers to inputs
      mwpointer x3d_pr     !collect all 2d input arrays into 1 3d array (input as vector)
      mwpointer x4d_pr     !pointer to input of Sdir
      mwpointer xmesh_pr   !pointer to inputs on mesh
      mwpointer ip_pr      !pointer to int_prams
      mwpointer rp_pr      !pointer to real_prams
      mwpointer Tp_pr      !pointer to Tp  (used if only single freq)
      mwpointer MD_pr      !pointer to mwd (used if only single dir)
      mwpointer nmesh_e_pr !pointer to no of mesh elements

      integer,parameter :: nxy4d      = n_wavdir*n_wave_freq   ! size(xy4d) = ii*jj*nxy4d
      integer,parameter :: Nin        = 8                      ! no of inputs
      integer,parameter :: nx3d       = 3                      ! size(x3d) = ii*jj*nx3d -- ice_arrays rearranged to a vector
      integer,parameter :: ni         = 9                      ! length of int_prams
      integer,parameter :: nr         = 4                      ! length of real_prams

      mwsize,dimension (Nin)  :: mvec_in,nvec_in,size_vec_in
      integer,dimension(Nin)  :: mvec_in_0,nvec_in_0,size_vec_in_0

      !!inputs from matlab - need to be double precision real
      real*8,dimension(ii*jj*nx3d)  :: x3d
      real*8,dimension(ii*jj*nxy4d) :: x4d
      real*8, allocatable           :: xmesh(:)
      real*8,dimension(ni)          :: ip
      real*8,dimension(nr)          :: rp
      real*8                        :: Tp_io
      real*8                        :: mwd_io
      real*8                        :: nmesh_e
C-----------------------------------------------------------------------


C-----OUTPUTS-----------------------------------------------------------
      integer,parameter :: Nout  = 3   ! no of outputs

      !!pointers to outputs
      mwpointer y3d_pr!collect all 2d output arrays into 1 3d array (output as vector)
      mwpointer y4d_pr!pointer to output of Sdir
      mwpointer ymesh_pr!pointer to output on mesh

      integer,parameter       :: ny3d  = 6   ! size(y3d) = ii*jj*ny3d
      mwsize,dimension(Nout)  :: mvec_out,nvec_out,size_vec_out
      mwsize                  :: nmesh

      !! outputs from main subroutine (float)
      !! - need to be double to pass to matlab (plhs)
      !! - use dble
      real*8,dimension(ii*jj*ny3d)  :: y3d
      real*8,dimension(ii*jj*nxy4d) :: y4d
      real*8,allocatable            :: ymesh(:)
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
      logical           :: critter
      logical,parameter :: PRINT_INFO = .false.
      character(len=12) :: cj,cm,cn,cm0,cn0,cdum
      character(len=80) :: msg


      ! =================================================
C     Check for proper number of arguments.
      if(nrhs.ne.Nin) then

         write(cdum,'(i2.2)') Nin
         msg  = trim(cdum)//' inputs required.\n'
         k     = mexPrintf(trim(msg))
         !!
         write(cdum,'(i2.2)') nrhs
         msg   = trim(cdum)//' inputs given.\n'
         k     = mexPrintf(trim(msg))
         !!
         call mexErrMsgTxt('Number of inputs')

      elseif(nlhs.ne.Nout) then

         write(cdum,'(i2.2)') Nout
         msg  = trim(cdum)//' outputs required.\n'
         k     = mexPrintf(trim(msg))
         !!
         write(cdum,'(i2.2)') nlhs
         msg   = trim(cdum)//' outputs requested.\n'
         k     = mexPrintf(trim(msg))
         !!
         call mexErrMsgTxt('Number of outputs')

      endif
      ! =================================================

      !call write_status()

      
      !don't know size of xmesh in advance
      nmesh = mxGetM(prhs(3)) !mwsize (integer*8)
      allocate(xmesh(nmesh))
      nmeshi   = int(nmesh)   !integer*4

      !call write_status()

      ! =================================================
      ! expected sizes of inputs:
      mvec_in_0      = (/ii*jj*nxy4d   !Sdir      
     &                 ,ii*jj*nx3d     !ice_arrays
     &                 ,nmeshi         !xmesh
     &                 ,1              !int_prams
     &                 ,1              !real_prams
     &                 ,1              !Tp_io
     &                 ,1              !mwd_io
     &                 ,1/)            !nmesh_e

      nvec_in_0      = (/1    !Sdir      
     &                 ,1     !ice_arrays
     &                 ,1     !xmesh
     &                 ,ni    !int_prams
     &                 ,nr    !real_prams
     &                 ,1     !Tp_io
     &                 ,1     !mwd_io
     &                 ,1/)   !nmesh_e

      size_vec_in_0  = mvec_in_0*nvec_in_0
      do j=1,Nin

         !test if input is numeric
         if(mxIsNumeric(prhs(j)) .eq. 0) then
            k  = mexPrintf('Input '//cj//' not numeric.')
            call mexErrMsgTxt('Input '//cj//' not numeric.')
         endif

         !! get actual input sizes
         mvec_in(j)     = mxGetM(prhs(j))
         nvec_in(j)     = mxGetN(prhs(j))
         size_vec_in(j) = mvec_in(j)*nvec_in(j)

         !! for error printout (if there is one)
         write(cj,'(i1.1)')      j
         write(cm0,'(i12.12)')   mvec_in_0(j) !expected M
         write(cn0,'(i4.4)')     nvec_in_0(j) !expected N
         write(cm,'(i12.12)')    mvec_in(j)   !given M
         write(cn,'(i4.4)')      nvec_in(j)   !given N

         !! test sizes
         critter  = (    (int(mvec_in(j)).ne.mvec_in_0(j))              &
     &               .or.(int(nvec_in(j)).ne.nvec_in_0(j))    )

         if (critter) then
            k  = mexPrintf('Input '//cj//' dimensions:\n')
            k  = mexPrintf(cm//' x '//cn//'\n')
            k  = mexPrintf('Should be:'//'\n')
            k  = mexPrintf(cm0//' x '//cn0//'\n')
            call mexErrMsgTxt('Wrong dimensions for input '//cj)
         endif

      end do

      if (PRINT_INFO) then
         k  = mexPrintf('All inputs correct sizes...\n')
         k  = mexPrintf('All inputs numeric...\n')
      endif
      ! =================================================

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!pointers to input memory
      x4d_pr      = mxGetPr(prhs(1))
      x3d_pr      = mxGetPr(prhs(2))
      xmesh_pr    = mxGetPr(prhs(3))
      ip_pr       = mxGetPr(prhs(4))
      rp_pr       = mxGetPr(prhs(5))
      Tp_pr       = mxGetPr(prhs(6))
      MD_pr       = mxGetPr(prhs(7))
      nmesh_e_pr  = mxGetPr(prhs(8))

      !!copy input memory to fortran arrays
      call mxCopyPtrToReal8(x4d_pr    ,x4d    ,size_vec_in(1))
      call mxCopyPtrToReal8(x3d_pr    ,x3d    ,size_vec_in(2))
      call mxCopyPtrToReal8(xmesh_pr  ,xmesh  ,size_vec_in(3))
      call mxCopyPtrToReal8(ip_pr     ,ip     ,size_vec_in(4))
      call mxCopyPtrToReal8(rp_pr     ,rp     ,size_vec_in(5))
      call mxCopyPtrToReal8(Tp_pr     ,Tp_io  ,size_vec_in(6))
      call mxCopyPtrToReal8(MD_pr     ,mwd_io ,size_vec_in(7))
      call mxCopyPtrToReal8(nmesh_e_pr,nmesh_e,size_vec_in(8))

      if (PRINT_INFO) then
         k  = mexPrintf('All input pointers copied to arrays...\n')
      endif
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     Create matrices for the return arguments.
      mvec_out       = (/ii*jj*nxy4d   !y4d
     &                  ,ii*jj*ny3d    !y3d
     &                  ,nmeshi/)      !ymesh

      nvec_out       = (/1    !y4d
     &                  ,1    !y3d
     &                  ,1/)  !ymesh

      size_vec_out   = mvec_out*nvec_out
      do j=1,Nout
         plhs(j)  = mxCreateDoubleMatrix(mvec_out(j),nvec_out(j),0)
      end do

      !!pointers to output memory
      y4d_pr   = mxGetPr(plhs(1))
      y3d_pr   = mxGetPr(plhs(2))
      ymesh_pr = mxGetPr(plhs(3))

      !! initialise output arrays
      y3d   = 0.0

      !! in/out arrays
      allocate(ymesh(nmesh))
      ymesh = xmesh
      y4d   = x4d

      if (PRINT_INFO) then
         k  = mexPrintf('Initialised output arrays...\n')
         k  = mexPrintf('Have pointers, about to do computations...\n')
         k  = mexPrintf('Entering mex_io_gate_vSdir_mesh...\n')
      end if


      call mex_io_gate_vSdir_mesh(y4d,y3d,ymesh
     &                    ,x3d,nint(ip),rp
     &                    ,Tp_io,mwd_io
     &                    ,ni,nr,nx3d,ny3d,nxy4d,nint(nmesh_e))

      if (PRINT_INFO) then
         k  = mexPrintf('Finished computations...\n')
      endif

C     Load the data into pointers, which are the outputs to MATLAB
      call mxCopyReal8ToPtr(y4d,y4d_pr,size_vec_out(1))
      call mxCopyReal8ToPtr(y3d,y3d_pr,size_vec_out(2))
      call mxCopyReal8ToPtr(ymesh,ymesh_pr,size_vec_out(3))

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      deallocate(xmesh)
      deallocate(ymesh)
      return

      end subroutine mexFunction


      subroutine write_status()
      implicit none

      open(unit=20,file='status.txt',status='replace')
      write(20,'(a)') 'Mex reached this point'
      close(20)

      end subroutine write_status
