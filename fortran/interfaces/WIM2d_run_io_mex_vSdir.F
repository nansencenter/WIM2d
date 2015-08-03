#include "fintrf.h"
C
#if 0
C
C     CALL:
C     [Sdir_io,out_arrays]=...
C        WIM2d_run_io_mex_vSdir(Sdir_io,ice_arrays,int_prams,...
C                               real_prams,Tp_io,mwd_io);
C
#endif
C
C     This is a MEX-file for MATLAB.

      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
C-----------------------------------------------------------------------
C
      use mod_WIM2d_run, only: mex_io_gate_vSdir
      implicit none

      include "grid_info.h"
      include "wave_info.h"

      mwpointer plhs(*), prhs(*)
      mwpointer mxCreateDoubleMatrix, mxGetPr
      integer*4   :: mexPrintf,k

      !!pointers to inputs
      mwpointer x3d_pr!collect all 2d input arrays into 1 3d array (input as vector)
      mwpointer x4d_pr!pointer to input of Sdir
      mwpointer ip_pr
      mwpointer rp_pr
      mwpointer Tp_pr
      mwpointer MD_pr
      integer,parameter :: nxy4d  = n_wavdir*n_wave_freq  ! size(xy4d) = ii*jj*nxy4d
      integer,parameter :: Nin   = 6   ! no of inputs
      integer,parameter :: nx3d  = 3   ! size(x3d) = ii*jj*nx3d -- ice_arrays rearranged to a vector
      integer,parameter :: ni    = 8   ! length of int_prams
      integer,parameter :: nr    = 3   ! length of real_prams

      !!pointers to outputs
      integer,parameter :: Nout  = 2   ! no of outputs

      mwpointer y3d_pr!collect all 2d output arrays into 1 3d array (output as vector)
      integer,parameter :: ny3d  = 5   ! size(y3d) = ii*jj*ny3d

      mwpointer y4d_pr!pointer to output of Sdir


C-----------------------------------------------------------------------
C

      integer nlhs, nrhs, mxIsNumeric
      mwsize mxGetM, mxGetN
      mwsize m0,n0,size0   !dimensions of x4d,y4d     (in,out)
      mwsize m,n,size      !dimensions of x3d         (in)
      mwsize m2,n2,size2   !dimensions of int_prams   (in)
      mwsize m3,n3,size3   !dimensions of real_prams  (in)
      mwsize m4,n4,size4   !dimensions of y3d         (out)
      mwsize m5,n5,size5   !dimensions of Tp_io       (in)
      mwsize m6,n6,size6   !dimensions of mwd_io      (in)
      integer,parameter :: DO_TEST  = 0   !!test i-o for reals,integers

      !!inputs from matlab - need to be double precision
      real*8,dimension(ii*jj*nx3d)  :: x3d
      real*8,dimension(ii*jj*nxy4d) :: x4d
      real*8,dimension(ni)          :: ip
      real*8,dimension(nr)          :: rp
      !!
      real*8   :: Tp_io,mwd_io
      integer  :: j,n3d,n4d
      logical  :: critter

      !! outputs from main subroutine (float)
      !! - need to be converted to double before passing to matlab (plhs)
      !! - use dble
      real,dimension(ii*jj*ny3d)    :: y3d
      real,dimension(ii*jj*nxy4d)   :: y4d
      character(len=12)             :: cii,cjj,cj,cm,cn,cdum
      character(len=80)             :: msg

C     Check for proper number of arguments.
      if(nrhs .ne. Nin) then
         write(cdum,'(i2.2)') Nin
         call mexErrMsgTxt(trim(cdum)//' inputs required.')
      elseif(nlhs .ne. Nout) then
         write(cdum,'(i2.2)') Nout
         call mexErrMsgTxt(trim(cdum)//' outputs required.')
      endif

      m0    = mxGetM(prhs(1))
      n0    = mxGetN(prhs(1))
      size0 = m0*n0
      m     = mxGetM(prhs(2))
      n     = mxGetN(prhs(2))
      size  = m*n
      m2    = mxGetM(prhs(3))
      n2    = mxGetN(prhs(3))
      size2 = m2*n2
      m3    = mxGetM(prhs(4))
      n3    = mxGetN(prhs(4))
      size3 = m3*n3
      m5    = mxGetM(prhs(5))
      n5    = mxGetN(prhs(5))
      size5 = m5*n5
      m6    = mxGetM(prhs(6))
      n6    = mxGetN(prhs(6))
      size6 = m6*n6

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     Get the size of the input arrays.

      !! x4d (Sdir_io)
      do j=1,1

         n4d   = ii*jj*nxy4d
         write(cdum,'(i12.12)') n4d
         write(cii,'(i4.4)') ii
         write(cii,'(i4.4)') ii
         write(cjj,'(i4.4)') jj
         write(cj,'(i1.1)')  j
         write(cm,'(i12.12)')  mxGetM(prhs(j))
         write(cn,'(i4.4)')  mxGetN(prhs(j))

C        matrices should be ii*jj
         critter  = ((mxGetM(prhs(j)).ne.n4d)                           &
     &               .or.(mxGetN(prhs(j)).ne.1))
         if (critter) then
            k  = mexPrintf('Input '//cj//' dimensions:\n')
            k  = mexPrintf(cm//' x '//cn//'\n')
            k  = mexPrintf('Should be:'//'\n')
            k  = mexPrintf(cdum//' x 1\n')
            call mexErrMsgTxt('Wrong dimensions for input '//cj)
         endif
      end do

      !! x3d (ice_arrays)
      do j=2,2

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
            call mexErrMsgTxt('Wrong dimensions for input '//cj)
         endif

      end do

      !! int_prams
      do j=3,3
         write(cii,'(i4.4)') 1
         write(cjj,'(i4.4)') ni
         write(cj,'(i1.1)')  j
         write(cm,'(i1.1)')  m2
         write(cn,'(i1.1)')  n2

C        matrices should be 1*ni
         if((m2.ne.1).or.(n2.ne.ni)) then
            k  = mexPrintf('Input '//cj//' dimensions:\n')
            k  = mexPrintf(cm//' x '//cn//'\n')
            k  = mexPrintf('Should be:'//'\n')
            k  = mexPrintf(cii//' x '//cjj//'\n')
            call mexErrMsgTxt('Wrong dimensions for input '//cj)
         endif
      end do

      !! real_prams
      do j=4,4
         write(cii,'(i4.4)') 1
         write(cjj,'(i4.4)') nr
         write(cj,'(i1.1)')  j
         write(cm,'(i1.1)')  m3
         write(cn,'(i1.1)')  n3

C        matrices should be 1*nr
         if((m3.ne.1).or.(n3.ne.nr)) then
            k  = mexPrintf('Input '//cj//' dimensions:\n')
            k  = mexPrintf(cm//' x '//cn//'\n')
            k  = mexPrintf('Should be:'//'\n')
            k  = mexPrintf(cii//' x '//cjj//'\n')
            call mexErrMsgTxt('Wrong dimensions for input '//cj)
         endif
      end do

      !! Tp_io
      do j=5,5
         write(cii,'(i4.4)') 1
         write(cjj,'(i4.4)') 1
         write(cj,'(i1.1)')  j
         write(cm,'(i1.1)')  m5
         write(cn,'(i1.1)')  n5

C        matrices should be 1*1
         if((m5.ne.1).or.(n5.ne.1)) then
            k  = mexPrintf('Input '//cj//' dimensions:\n')
            k  = mexPrintf(cm//' x '//cn//'\n')
            k  = mexPrintf('Should be:'//'\n')
            k  = mexPrintf(cii//' x '//cjj//'\n')
            call mexErrMsgTxt('Wrong dimensions for input '//cj)
         endif
      end do

      !! mwd_io
      do j=6,6
         write(cii,'(i4.4)') 1
         write(cjj,'(i4.4)') 1
         write(cj,'(i1.1)')  j
         write(cm,'(i1.1)')  m6
         write(cn,'(i1.1)')  n6

C        matrices should be 1*1
         if((m6.ne.1).or.(n6.ne.1)) then
            k  = mexPrintf('Input '//cj//' dimensions:\n')
            k  = mexPrintf(cm//' x '//cn//'\n')
            k  = mexPrintf('Should be:'//'\n')
            k  = mexPrintf(cii//' x '//cjj//'\n')
            call mexErrMsgTxt('Wrong dimensions for input '//cj)
         endif
      end do
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !k  = mexPrintf('All inputs correct sizes...\n')

C     Check to ensure the array is numeric (not strings).
      do j=1,Nin
         write(cj,'(i1.1)')  j
         if(mxIsNumeric(prhs(j)) .eq. 0) then
            call mexErrMsgTxt('Input '//cj//' must be a numeric array.')
         endif
      end do
      !k  = mexPrintf('All inputs numeric...\n')

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!pointers to input memory
      x4d_pr   = mxGetPr(prhs(1))
      x3d_pr   = mxGetPr(prhs(2))
      ip_pr    = mxGetPr(prhs(3))
      rp_pr    = mxGetPr(prhs(4))
      Tp_pr    = mxGetPr(prhs(5))
      MD_pr    = mxGetPr(prhs(6))

      !!copy input memory to fortran arrays
      call mxCopyPtrToReal8(x4d_pr,x4d,size0)
      call mxCopyPtrToReal8(x3d_pr,x3d,size)
      call mxCopyPtrToReal8(ip_pr,ip,size2)
      call mxCopyPtrToReal8(rp_pr,rp,size3)
      call mxCopyPtrToReal8(Tp_pr,Tp_io,size5)
      call mxCopyPtrToReal8(MD_pr,mwd_io,size6)
      y4d   = real(x4d)!cast to real*4: this goes into the main function with intent(inout)
      !!
      !k  = mexPrintf('All pointers copied to arrays...\n')
!     write(msg,'(f12.12)') Tp_io
!     k  = mexPrintf('T for single freq (s): '//trim(msg)//'\n')
!     write(msg,'(e7.3)') mwd_io
!     k  = mexPrintf('mwd for single dir (degrees-from): '//            &
!    &               trim(msg)//'\n')
      

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     Create matrix for the return arguments.
      plhs(1)  = mxCreateDoubleMatrix(m0,n0,0)
      !!
      m4       = ii*jj*ny3d
      n4       = 1
      size4    = m4*n4
      plhs(2)  = mxCreateDoubleMatrix(m4,n4,0)

      !!pointers to output memory
      y4d_pr   = mxGetPr(plhs(1))
      y3d_pr   = mxGetPr(plhs(2))

      !! initialise output arrays (NB not in/out arrays!)
      y3d = 0.0

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !k  = mexPrintf('Initialised output arrays...\n')

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !k  = mexPrintf('Have pointers, about to do computations...\n')
      !k  = mexPrintf('Entering mex_io_gate_vSdir...\n')

      call mex_io_gate_vSdir(y4d,y3d,real(x3d),nint(ip),real(rp)        &
     &                    ,real(Tp_io),real(mwd_io)                     &
     &                    ,ni,nr,nx3d,ny3d,nxy4d)

      k  = mexPrintf('Finished computations...\n')
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     Load the data into pointers, which are the outputs to MATLAB
      call mxCopyReal8ToPtr(dble(y3d),y3d_pr,size4)
      call mxCopyReal8ToPtr(dble(y4d),y4d_pr,size0)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      return
      end subroutine mexFunction
