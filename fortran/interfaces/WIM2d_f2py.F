      !WIM2d_f2py.F
      !interface to python
      !Author: Timothy Williams
      !Date:   20141118, 17:20:34 CET

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine py_wim2d_run

      use mod_wim2d_run, only: wim2d_run

      call wim2d_run()

      end subroutine py_wim2d_run
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine py_wim2d_run_io(in_arrays,out_arrays                   &
     &                  ,int_prams,real_prams,nx,ny)

      use mod_wim2d_run, only: wim2d_run_io

      implicit none

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!variables
      integer,intent(in)                  :: nx,ny,int_prams(2)
      real,intent(in)                     :: real_prams(3)
      real,dimension(nx,ny,6),intent(in)  :: in_arrays!!c,h,Dmax,Hs,Tp,mwd
      real,dimension(nx,ny,5),intent(out) :: out_arrays!!Dmax,Hs,Tp,taux,tauy
!f2py intent(hide)   :: nx,ny

      call wim2d_run_io(in_arrays,out_arrays,int_prams,real_prams,nx,ny)

      end subroutine py_wim2d_run_io
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
