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
      integer nlhs, nrhs, mxIsNumeric
      mwsize mxGetM, mxGetN
      integer*4   :: mexPrintf,k,j

C-----INPUTS------------------------------------------------------------
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
      integer,parameter :: ni    = 9   ! length of int_prams
      integer,parameter :: nr    = 4   ! length of real_prams

      mwsize,dimension (Nin)  :: mvec_in,nvec_in,size_vec_in
      integer,dimension(Nin)  :: mvec_in_0,nvec_in_0,size_vec_in_0

      !!inputs from matlab - need to be double precision
      real*8,dimension(ii*jj*nx3d)  :: x3d
      real*8,dimension(ii*jj*nxy4d) :: x4d
      real*8,dimension(ni)          :: ip
      real*8,dimension(nr)          :: rp
      real*8                        :: Tp_io,mwd_io
C-----------------------------------------------------------------------


C-----OUTPUTS-----------------------------------------------------------
      integer,parameter :: Nout  = 2   ! no of outputs

      !!pointers to outputs
      mwpointer y3d_pr!collect all 2d output arrays into 1 3d array (output as vector)
      mwpointer y4d_pr!pointer to output of Sdir

      integer,parameter       :: ny3d  = 5   ! size(y3d) = ii*jj*ny3d
      mwsize,dimension(Nout)  :: mvec_out,nvec_out,size_vec_out

      !! outputs from main subroutine (float)
      !! - need to be converted to double before passing to matlab (plhs)
      !! - use dble
      real,dimension(ii*jj*ny3d)    :: y3d
      real,dimension(ii*jj*nxy4d)   :: y4d
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
      logical           :: critter
      logical,parameter :: PRINT_INFO = .false.
      character(len=12) :: cj,cm,cn,cm0,cn0,cdum
      character(len=80) :: msg


      ! =================================================
C     Check for proper number of arguments.
      if(nrhs .ne. Nin) then
         write(cdum,'(i2.2)') Nin
         call mexErrMsgTxt(trim(cdum)//' inputs required.')
      elseif(nlhs .ne. Nout) then
         write(cdum,'(i2.2)') Nout
         call mexErrMsgTxt(trim(cdum)//' outputs required.')
      endif
      ! =================================================


      ! =================================================
      ! expected sizes of inputs:
      ! Sdir,ice_arrays,int_prams,real_prams,Tp_io,mwd_io
      mvec_in_0      = (/ii*jj*nxy4d,ii*jj*nx3d,1 ,1 ,1,1/)
      nvec_in_0      = (/1    ,1   ,ni,nr,1,1/)
      size_vec_in_0  = mvec_in_0*nvec_in_0

      do j=1,Nin
         !! get actual input sizes
         mvec_in(j)     = mxGetM(prhs(j))
         nvec_in(j)     = mxGetN(prhs(j))
         size_vec_in(j) = mvec_in(j)*nvec_in(j)

         !! for error printout (if there is one)
         write(cj,'(i1.1)')         j
         write(cm0,'(i4.4)')        mvec_in_0(j)
         write(cn0,'(i4.4)')        nvec_in_0(j)
         write(cm,'(i12.12)')       mvec_in(j)
         write(cn,'(i4.4)')         nvec_in(j)

         !! test sizes
         critter  = (    (mvec_in(j).ne.mvec_in_0(j))                   &
     &               .or.(nvec_in(j).ne.nvec_in_0(j))    )

         if (critter) then
            k  = mexPrintf('Input '//cj//' dimensions:\n')
            k  = mexPrintf(cm//' x '//cn//'\n')
            k  = mexPrintf('Should be:'//'\n')
            k  = mexPrintf(cm0//' x '//cn0//'\n')
            call mexErrMsgTxt('Wrong dimensions for input '//cj)
         endif

         !test if input is numeric
         if(mxIsNumeric(prhs(j)) .eq. 0) then
            call mexErrMsgTxt('Input '//cj//' must be a numeric array.')
         endif
      end do

      if (PRINT_INFO) then
         k  = mexPrintf('All inputs correct sizes...\n')
         k  = mexPrintf('All inputs numeric...\n')
      endif
      ! =================================================

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!pointers to input memory
      x4d_pr   = mxGetPr(prhs(1))
      x3d_pr   = mxGetPr(prhs(2))
      ip_pr    = mxGetPr(prhs(3))
      rp_pr    = mxGetPr(prhs(4))
      Tp_pr    = mxGetPr(prhs(5))
      MD_pr    = mxGetPr(prhs(6))

      !!copy input memory to fortran arrays
      call mxCopyPtrToReal8(x4d_pr,x4d    ,size_vec_in(1))
      call mxCopyPtrToReal8(x3d_pr,x3d    ,size_vec_in(2))
      call mxCopyPtrToReal8(ip_pr,ip      ,size_vec_in(3))
      call mxCopyPtrToReal8(rp_pr,rp      ,size_vec_in(4))
      call mxCopyPtrToReal8(Tp_pr,Tp_io   ,size_vec_in(5))
      call mxCopyPtrToReal8(MD_pr,mwd_io  ,size_vec_in(6))
      y4d   = real(x4d)!cast to real*4: this goes into the main function with intent(inout)

      if (PRINT_INFO) then
         k  = mexPrintf('All input pointers copied to arrays...\n')
      endif
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     Create matrices for the return arguments.
      mvec_out       = (/ii*jj*nxy4d,ii*jj*ny3d/)
      nvec_out       = (/1,1/)
      size_vec_out   = mvec_out*nvec_out

      do j=1,Nout
         plhs(j)  = mxCreateDoubleMatrix(mvec_out(j),nvec_out(j),0)
      end do

      !!pointers to output memory
      y4d_pr   = mxGetPr(plhs(1))
      y3d_pr   = mxGetPr(plhs(2))

      !! initialise output arrays (NB not in/out arrays!)
      y3d = 0.0

      if (PRINT_INFO) then
         k  = mexPrintf('Initialised output arrays...\n')
         k  = mexPrintf('Have pointers, about to do computations...\n')
         k  = mexPrintf('Entering mex_io_gate_vSdir...\n')
      end if

      call mex_io_gate_vSdir(y4d,y3d,real(x3d),nint(ip),real(rp)        &
     &                    ,real(Tp_io),real(mwd_io)                     &
     &                    ,ni,nr,nx3d,ny3d,nxy4d)
c     subroutine mex_io_gate_vSdir(xy4d,y3d,x3d,int_prams,real_prams    &
c    &                  ,Tp_io,mwd_io,ni,nr,nx3d,ny3d,nxy4d)

      if (PRINT_INFO) then
         k  = mexPrintf('Finished computations...\n')
      endif

C     Load the data into pointers, which are the outputs to MATLAB
      call mxCopyReal8ToPtr(dble(y4d),y4d_pr,size_vec_out(1))
      call mxCopyReal8ToPtr(dble(y3d),y3d_pr,size_vec_out(2))
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      return
      end subroutine mexFunction
