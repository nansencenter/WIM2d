      subroutine atten_young(outputs,inputs,kguess)
      !wrapper for f2py for RTparam_outer (make sure doubles are changed to reals if necessary)

      use mod_atten_young, only: RTparam_outer
      implicit none
      include "common_wim.h"

      real*8,intent(in)   :: inputs(5)
         !(/h,om,young,drag_rp,visc_ws/)
      real*8,intent(in),optional :: kguess
      real*8,intent(out)  :: outputs(8)
         !(/damping,kice,kwtr,int_adm,ac,modT,argR,argT/)

      real :: damping,kice,kwtr,int_adm,ac,modT,argR,argT
      real :: om,guess

      !!set some params manually
      om                = real(inputs(2))
      wim_info%young    = real(inputs(3))
      wim_info%drag_rp  = real(inputs(4))
      wim_info%visc_ws  = real(inputs(5))
      if (present(kguess)) then
         guess = kguess
      else
         guess = om**2/wim_info%gravity
      end if

      call RTparam_outer(damping,kice,kwtr,int_adm,
     &                  ac,modT,argR,argT,
     &                  real(inputs(1:2)),guess)

      outputs  = dble((/damping,kice,kwtr,int_adm
     &                 ,ac,modT,argR,argT/))

      end subroutine atten_young
