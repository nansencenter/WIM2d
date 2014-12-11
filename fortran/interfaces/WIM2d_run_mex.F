#include "fintrf.h"
C
#if 0
C     
C     matsq.F
C     .F file needs to be preprocessed to generate .for equivalent
C     
#endif
C     

      subroutine mexFunction(nlhs, plhs, nrhs, prhs)

      use mod_wim2d_run, only: wim2d_run
C-----------------------------------------------------------------------
C
      mwpointer plhs(*), prhs(*)
C-----------------------------------------------------------------------
C

      integer nlhs, nrhs

C     Check for proper number of arguments. 
      if(nrhs .gt. 0) then
         call mexErrMsgTxt('No inputs.')
      elseif(nlhs .gt. 0) then
         call mexErrMsgTxt('No outputs.')
      endif

      !call main subroutine
      call wim2d_run()

      return
      end subroutine mexFunction
