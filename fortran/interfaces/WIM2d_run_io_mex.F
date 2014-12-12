#include "fintrf.h"
C
#if 0
C     
C     matsq.F
C     .F file needs to be preprocessed to generate .for equivalent
C     
#endif
C     
C     matsq.f
C
C     squares the input matrix
      
C     This is a MEX-file for MATLAB.
C     Copyright 1984-2006 The MathWorks, Inc.
C     $Revision: 1.13.2.3 $

      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
C-----------------------------------------------------------------------
C
      use mod_WIM2d_run, only: mex_io_gate
      implicit none

      mwpointer plhs(*), prhs(*)
      mwpointer mxCreateDoubleMatrix, mxGetPr

      !!pointers to inputs
      mwpointer x1_pr, x2_pr, x3_pr, x4_pr, x5_pr, x6_pr,ip_pr

      !!pointers to outputs
      mwpointer y1_pr, y2_pr, y3_pr, y4_pr, y5_pr, y6_pr
C-----------------------------------------------------------------------
C

      include "grid_info.h"

      integer nlhs, nrhs, mxIsNumeric
      mwsize mxGetM, mxGetN
      mwsize m, n, size
      integer,parameter :: ni = 2

      !!inputs from matlab - need to be double precision
      real*8   :: x1(ii*jj),x2(ii*jj),x3(ii*jj)                         &
     &            ,x4(ii*jj),x5(ii*jj),x6(ii*jj),ip(ni)
      integer  :: j,int_prams(ni)

      !! outputs from main subroutine (float)
      !! - need to be converted to double before passing to matlab (plhs)
      !! - use dble
      real   :: y1(ii*jj),y2(ii*jj),y3(ii*jj)                           &
     &         ,y4(ii*jj),y5(ii*jj),y6(ii*jj)
      character(len=4)   :: cii,cjj

C     Check for proper number of arguments. 
      if(nrhs .ne. 6) then
         call mexErrMsgTxt('6 inputs required.')
      elseif(nlhs .ne. 6) then
         call mexErrMsgTxt('6 outputs required.')
      endif

C     Get the size of the input arrays.
      m     = mxGetM(prhs(1))
      n     = mxGetN(prhs(1))
      size  = m*n
      write(cii,'(i4.4)') ii
      write(cjj,'(i4.4)') jj

C     matrices shoulb be ii*jj
      if((m.ne.ii).or.(n.ne.jj)) then
         call mexErrMsgTxt('Input dimensions should be '                &
     &                     //cii//' x '//cjj)
      endif
      
C     Check to insure the array is numeric (not strings).
      if(mxIsNumeric(prhs(1)) .eq. 0) then
         call mexErrMsgTxt('Input must be a numeric array.')
      endif

C     Create matrix for the return arguments.
      do j=1,5
         plhs(j)  = mxCreateDoubleMatrix(m,n,0)
      end do

      !!pointers to output memory
      y1_pr    = mxGetPr(plhs(1))
      y2_pr    = mxGetPr(plhs(2))
      y3_pr    = mxGetPr(plhs(3))
      y4_pr    = mxGetPr(plhs(4))
      y5_pr    = mxGetPr(plhs(5))
      y6_pr    = mxGetPr(plhs(6))

      !!pointers to input memory
      x1_pr    = mxGetPr(prhs(1))
      x2_pr    = mxGetPr(prhs(2))
      x3_pr    = mxGetPr(prhs(3))
      x4_pr    = mxGetPr(prhs(4))
      x5_pr    = mxGetPr(prhs(5))
      x6_pr    = mxGetPr(prhs(6))
      ip_pr    = mxGetPr(prhs(7))

      !!copy input memory to fortran arrays
      call mxCopyPtrToReal8(x1_pr,x1,size)
      call mxCopyPtrToReal8(x2_pr,x2,size)
      call mxCopyPtrToReal8(x3_pr,x3,size)
      call mxCopyPtrToReal8(x4_pr,x4,size)
      call mxCopyPtrToReal8(x5_pr,x5,size)
      call mxCopyPtrToReal8(x6_pr,x6,size)
      call mxCopyPtrToReal8(ip_pr,ip,ni)

C     Call the computational subroutine.
      int_prams   = nint(ip)
      call mex_io_gate(y1,y2,y3,y4,y5                                   &
     &                 ,real(x1),real(x2),real(x3)                      &
     &                 ,real(x4),real(x5),real(x6)                      &
     &                 ,int_prams,ii,jj,ni)

C     Load the data into pointers, which are the outputs to MATLAB
      call mxCopyReal8ToPtr(dble(y1),y1_pr,size)     
      call mxCopyReal8ToPtr(dble(y2),y2_pr,size)     
      call mxCopyReal8ToPtr(dble(y3),y3_pr,size)     
      call mxCopyReal8ToPtr(dble(y4),y4_pr,size)     
      call mxCopyReal8ToPtr(dble(y5),y5_pr,size)     

      return
      end subroutine mexFunction
