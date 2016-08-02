#include "fintrf.h"
C
#if 0
C     
C     CALL:
C     out_arrays=...
C        WIM2d_run_io_mex_vSdir(in_arrays,param_vec);
C     
#endif
C     
C     This is a MEX-file for MATLAB.

      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
C-----------------------------------------------------------------------
C
      use mod_WIM2d_run, only: mex_io_gate_v2,n_param_vec
      implicit none

      mwpointer plhs(*), prhs(*)
      mwpointer mxCreateDoubleMatrix, mxGetPr
      integer*4   :: mexPrintf,k

      !!pointers to inputs
      mwpointer x3d_pr!collect all 2d input arrays into 1 3d array (input as vector)
      mwpointer ip_pr
      mwpointer rp_pr
      integer,parameter :: Nin   = 2   ! no of inputs
      integer,parameter :: nx3d  = 6   ! size(x3d) = ii*jj*nx3d
      integer,parameter :: nr    = n_param_vec  ! length of real_prams

      !!pointers to outputs
      mwpointer y3d_pr!collect all 2d output arrays into 1 3d array (input as vector)
      integer,parameter :: Nout  = 1   ! no of outputs
      integer,parameter :: ny3d  = 6   ! size(y3d) = ii*jj*ny3d

C-----------------------------------------------------------------------
C

      include "grid_info.h"

      integer nlhs, nrhs, mxIsNumeric
      mwsize mxGetM, mxGetN
      mwsize m,n,size      !dimensions of x3d         (in)
      mwsize m3,n3,size3   !dimensions of real_prams  (in)
      mwsize m4,n4,size4   !dimensions of y3d         (out)
      integer,parameter :: DO_TEST  = 0   !!test i-o for reals,integers

      !!inputs from matlab - need to be double precision
      real*8,dimension(ii*jj*nx3d)  :: x3d
      real*8,dimension(nr)          :: rp
      integer  :: j,n3d
      logical  :: critter

      !! outputs from main subroutine (float)
      !! - need to be double to pass to matlab (plhs)
      real*8,dimension(ii*jj*ny3d) :: y3d
      character(len=12)          :: cii,cjj,cj,cm,cn,cdum
      character(len=80)          :: msg

C     Check for proper number of arguments. 
      if(nrhs .ne. Nin) then
         write(cdum,'(i2.2)') Nin
         call mexErrMsgTxt(trim(cdum)//' inputs required.')
      elseif(nlhs .ne. Nout) then
         write(cdum,'(i2.2)') Nout
         call mexErrMsgTxt(trim(cdum)//' outputs required.')
      endif

      m     = mxGetM(prhs(1))
      n     = mxGetN(prhs(1))
      size  = m*n
      m3    = mxGetM(prhs(2))
      n3    = mxGetN(prhs(2))
      size3 = m3*n3

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     Get the size of the input arrays.

      !call write_status()
      !! x3d
      do j=1,1

         n3d   = ii*jj*nx3d
         write(cdum,'(i12.12)') n3d
         write(cii,'(i4.4)') ii
         write(cii,'(i4.4)') ii
         write(cjj,'(i4.4)') jj
         write(cj,'(i1.1)')  j
         write(cm,'(i4.4)')  mxGetM(prhs(j))
         write(cn,'(i4.4)')  mxGetN(prhs(j))

C        matrices should be ii*jj
         critter  = ((mxGetM(prhs(j)).ne.n3d)                           &
     &               .or.(mxGetN(prhs(j)).ne.1))
         if (critter) then
            k  = mexPrintf('Input '//cj//' dimensions:\n')
            k  = mexPrintf(cm//' x '//cn//'\n')
            k  = mexPrintf('Should be:'//'\n')
            k  = mexPrintf(cdum//' x 1\n')
            call mexErrMsgTxt('Input '//cj//' dimensions should be '    &
     &                        //cii//' x '//cjj)
         endif
      end do

      !! real_prams
      do j=2,2
         write(cii,'(i4.4)') nr
         write(cjj,'(i4.4)') 1
         write(cj,'(i1.1)')  j
         write(cm,'(i4.4)')  m3
         write(cn,'(i4.4)')  n3

C        matrices should be 1*nr
         if((m3.ne.nr).or.(n3.ne.1)) then
            k  = mexPrintf('Input '//cj//' dimensions:\n')
            k  = mexPrintf(cm//' x '//cn//'\n')
            k  = mexPrintf('Should be:'//'\n')
            k  = mexPrintf(cii//' x '//cjj//'\n')
            call mexErrMsgTxt('Input '//cj//' wrong size')
         endif
      end do
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
C     Check to ensure the array is numeric (not strings).
      do j=1,Nin
         write(cj,'(i1.1)')  j
         if(mxIsNumeric(prhs(j)) .eq. 0) then
            k  = mexPrintf('Input '//cj//' must be a numeric array.\n')
            call mexErrMsgTxt('Input '//cj//' must be a numeric array.')
         endif
      end do

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!pointers to input memory
      x3d_pr   = mxGetPr(prhs(1))
      rp_pr    = mxGetPr(prhs(2))

      !!copy input memory to fortran arrays
      call mxCopyPtrToReal8(x3d_pr,x3d,size)
      call mxCopyPtrToReal8(rp_pr,rp,size3)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     Create matrix for the return arguments.
      m4       = ii*jj*ny3d
      n4       = 1
      size4    = m4*n4
      plhs(1)  = mxCreateDoubleMatrix(m4,n4,0)

      !!pointers to output memory
      y3d_pr   = mxGetPr(plhs(1))

      !! initialise output array
      y3d = 0.0
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      k  = mexPrintf('Have pointers, about to do computations...\n')
      k  = mexPrintf('Entering mex_io_gate_v2...\n')

      call mex_io_gate_v2(y3d,x3d,rp
     &                    ,ii,jj,nx3d,ny3d)

      k  = mexPrintf('Finished computations...\n')

C     Load the data into pointers, which are the outputs to MATLAB
      call mxCopyReal8ToPtr(y3d,y3d_pr,size4)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      return
      end subroutine mexFunction

      subroutine write_status()
      implicit none

      open(unit=20,file='status.txt',status='replace')
      write(20,'(a)') 'Mex reached this point'
      close(20)

      end subroutine write_status
