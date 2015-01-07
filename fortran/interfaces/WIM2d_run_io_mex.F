#include "fintrf.h"
C
#if 0
C     
C     CALL:
C     [Dmax,Hs,Tp,taux,tauy]=...
C        WIM2d_run_io_mex(c,h,Dmax,Hs,Tp,mwd,int_prams);
C        int_prams   = [SOLVER,ADV_DIM]
C     
#endif
C     
C     This is a MEX-file for MATLAB.

      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
C-----------------------------------------------------------------------
C
      use mod_WIM2d_run, only: mex_io_gate,mex_io_gate_test
      implicit none

      mwpointer plhs(*), prhs(*)
      mwpointer mxCreateDoubleMatrix, mxGetPr
      integer*4   :: mexPrintf,k

      !!pointers to inputs
      integer,parameter :: Nin   = 7
      mwpointer x1_pr
      mwpointer x2_pr
      mwpointer x3_pr
      mwpointer x4_pr
      mwpointer x5_pr
      mwpointer x6_pr
      mwpointer ip_pr
      mwpointer rp_pr

      !!pointers to outputs
      integer,parameter :: Nout  = 5
      mwpointer y1_pr
      mwpointer y2_pr
      mwpointer y3_pr
      mwpointer y4_pr
      mwpointer y5_pr
      
C-----------------------------------------------------------------------
C

      include "grid_info.h"

      integer nlhs, nrhs, mxIsNumeric
      mwsize mxGetM, mxGetN
      mwsize m,n,size
      mwsize m2,n2,size2
      mwsize m3,n3,size3
      integer,parameter :: ni       = 2
      integer,parameter :: nr       = 3
      integer,parameter :: DO_TEST  = 0   !!test i-o for reals,integers

      !!inputs from matlab - need to be double precision
      real*8,dimension(ii*jj) :: x1,x2,x3,x4,x5,x6
      real*8,dimension(ni)    :: ip
      real*8,dimension(nr)    :: rp
      integer  :: j
      logical  :: critter

      !! outputs from main subroutine (float)
      !! - need to be converted to double before passing to matlab (plhs)
      !! - use dble
      real,dimension(ii*jj)   :: y1,y2,y3,y4,y5
      integer,dimension(ni)   :: y6             !DO_TEST = 1
      character(len=4)        :: cii,cjj,cj,cm,cn
      character(len=80)       :: msg

C     Check for proper number of arguments. 
      if(nrhs .ne. Nin) then
         call mexErrMsgTxt('7 inputs required.')
      elseif(nlhs .ne. Nout) then
         call mexErrMsgTxt('5 outputs required.')
      endif

      m     = mxGetM(prhs(1))
      n     = mxGetN(prhs(1))
      size  = m*n
      m2    = mxGetM(prhs(7))
      n2    = mxGetN(prhs(7))
      size2 = m2*n2
      m3    = mxGetM(prhs(8))
      n3    = mxGetN(prhs(8))
      size3 = m3*n3

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     Get the size of the input arrays.
      do j=8,8
         write(cii,'(i4.4)') 1
         write(cjj,'(i4.4)') nr
         write(cj,'(i1.1)')  j
         write(cm,'(i1.1)')  m3
         write(cn,'(i1.1)')  n3

C        matrices should be 1*nr
            k  = mexPrintf('Input '//cj//' dimensions:\n')
            k  = mexPrintf(cm//' x '//cn//'\n')
            k  = mexPrintf('Should be:'//'\n')
            k  = mexPrintf(cii//' x '//cjj//'\n')
         if((m3.ne.1).or.(n3.ne.nr)) then
            k  = mexPrintf('Input '//cj//' dimensions:\n')
            k  = mexPrintf(cm//' x '//cn//'\n')
            k  = mexPrintf('Should be:'//'\n')
            k  = mexPrintf(cii//' x '//cjj//'\n')
            call mexErrMsgTxt('Input '//cj//' wrong size')
         endif
      end do
      do j=7,7
         write(cii,'(i4.4)') 1
         write(cjj,'(i4.4)') ni
         write(cj,'(i1.1)')  j
         write(cm,'(i1.1)')  m2
         write(cn,'(i1.1)')  n2

C        matrices should be 1*ni
            k  = mexPrintf('Input '//cj//' dimensions:\n')
            k  = mexPrintf(cm//' x '//cn//'\n')
            k  = mexPrintf('Should be:'//'\n')
            k  = mexPrintf(cii//' x '//cjj//'\n')
         if((m2.ne.1).or.(n2.ne.ni)) then
            k  = mexPrintf('Input '//cj//' dimensions:\n')
            k  = mexPrintf(cm//' x '//cn//'\n')
            k  = mexPrintf('Should be:'//'\n')
            k  = mexPrintf(cii//' x '//cjj//'\n')
            call mexErrMsgTxt('Input '//cj//' wrong size')
         endif
      end do
      do j=1,6
         write(cii,'(i4.4)') ii
         write(cjj,'(i4.4)') jj
         write(cj,'(i1.1)')  j
         write(cm,'(i4.4)')  mxGetM(prhs(j))
         write(cn,'(i4.4)')  mxGetN(prhs(j))

C        matrices should be ii*jj
         critter  = ((mxGetM(prhs(j)).ne.ii)                            &
     &               .or.(mxGetN(prhs(j)).ne.jj))
         if (critter) then
            k  = mexPrintf('Input '//cj//' dimensions:\n')
            k  = mexPrintf(cm//' x '//cn//'\n')
            k  = mexPrintf('Should be:'//'\n')
            k  = mexPrintf(cii//' x '//cjj//'\n')
            call mexErrMsgTxt('Input '//cj//' dimensions should be '    &
     &                        //cii//' x '//cjj)
         endif
      end do
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
C     Check to ensure the array is numeric (not strings).
      do j=1,8
         write(cj,'(i1.1)')  j
         if(mxIsNumeric(prhs(j)) .eq. 0) then
            call mexErrMsgTxt('Input '//cj//' must be a numeric array.')
         endif
      end do

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!pointers to input memory
      x1_pr    = mxGetPr(prhs(1))
      x2_pr    = mxGetPr(prhs(2))
      x3_pr    = mxGetPr(prhs(3))
      x4_pr    = mxGetPr(prhs(4))
      x5_pr    = mxGetPr(prhs(5))
      x6_pr    = mxGetPr(prhs(6))
      ip_pr    = mxGetPr(prhs(7))
      rp_pr    = mxGetPr(prhs(8))

      !!copy input memory to fortran arrays
      call mxCopyPtrToReal8(x1_pr,x1,size)
      call mxCopyPtrToReal8(x2_pr,x2,size)
      call mxCopyPtrToReal8(x3_pr,x3,size)
      call mxCopyPtrToReal8(x4_pr,x4,size)
      call mxCopyPtrToReal8(x5_pr,x5,size)
      call mxCopyPtrToReal8(x6_pr,x6,size)
      call mxCopyPtrToReal8(ip_pr,ip,size2)
      call mxCopyPtrToReal8(rp_pr,rp,size2)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     Create matrix for the return arguments.
      do j=1,4
         plhs(j)  = mxCreateDoubleMatrix(m,n,0)
      end do
      if (DO_TEST.eq.0) then
         plhs(5)  = mxCreateDoubleMatrix(m,n,0)
      else
         plhs(5)  = mxCreateDoubleMatrix(m2,n2,0)
      end if

      !!pointers to output memory
      y1_pr    = mxGetPr(plhs(1))
      y2_pr    = mxGetPr(plhs(2))
      y3_pr    = mxGetPr(plhs(3))
      y4_pr    = mxGetPr(plhs(4))
      y5_pr    = mxGetPr(plhs(5))
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      k  = mexPrintf('Have pointers, about to do computations...\n')
C     Call the computational subroutine.
      y1 = 0.0
      y2 = 0.0
      y3 = 0.0
      y4 = 0.0
      y5 = 0.0
      y6 = 0

      if (.false.) then
         !!simple test
         k  = mexPrintf('Doing simple test...\n')
         y1 = real(x1)**2
         y2 = real(x2)**2
         y3 = real(x3)**2
         y4 = real(x4)**2+real(x5)**2+real(x6)**2
         y5 = real(x5)**2
         !!
         write(msg,'(a,i2.2,a,i2.2,a)') 'int_prams = (',                &
     &                  nint(ip(1)),',',nint(ip(2)),')\n'
         k  = mexPrintf(trim(msg))
         k  = mexPrintf('\n')
         !!
         y6 = nint(ip)**2
      elseif (.false.) then
         k  = mexPrintf('Entering mex_io_gate_test...\n')
         call mex_io_gate_test(y1,y2,y3,y4,y5                           &
     &                         ,real(x1),real(x2),real(x3)              &
     &                         ,real(x4),real(x5),real(x6)              &
     &                         ,nint(ip),ii,jj,ni)
      else
         k  = mexPrintf('Entering mex_io_gate...\n')
         call mex_io_gate(y1,y2,y3,y4,y5                                &
     &                    ,real(x1),real(x2),real(x3)                   &
     &                    ,real(x4),real(x5),real(x6)                   &
     &                    ,nint(ip),real(rp),ii,jj,ni,nr)
      end if
!     if (DO_TEST.eq.0) then
!        k  = mexPrintf('Entering mex_io_gate...')
!        call mex_io_gate(y1,y2,y3,y4,y5                                &
!    &                    ,real(x1),real(x2),real(x3)                   &
!    &                    ,real(x4),real(x5),real(x6)                   &
!    &                    ,nint(ip),ii,jj,ni)
!     elseif (.false.) then
!        k  = mexPrintf('Entering mex_io_gate_test...')
!        call mex_io_gate_test(y1,y2,y3,y4,y6                           &
!    &                         ,real(x1),real(x2),real(x3)              &
!    &                         ,real(x4),real(x5),real(x6)              &
!    &                         ,nint(ip),ii,jj,ni)
!     else
!        k  = mexPrintf('Doing simple test...')
!        y1 = real(x1)**2
!        y2 = real(x2)**2
!        y3 = real(x3)**2
!        y4 = real(x4)**2+real(x5)**2+real(x6)**2
!        y6 = nint(ip)**2
!     end if
      k  = mexPrintf('Finished computations...\n')
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     Load the data into pointers, which are the outputs to MATLAB
      call mxCopyReal8ToPtr(dble(y1),y1_pr,size)
      call mxCopyReal8ToPtr(dble(y2),y2_pr,size)
      call mxCopyReal8ToPtr(dble(y3),y3_pr,size)
      call mxCopyReal8ToPtr(dble(y4),y4_pr,size)
      if (DO_TEST.eq.0) then
         call mxCopyReal8ToPtr(dble(y5),y5_pr,size)
      else
         call mxCopyReal8ToPtr(dble(y6),y5_pr,size2)
      end if
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      return
      end subroutine mexFunction
