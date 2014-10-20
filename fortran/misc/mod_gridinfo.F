      !! test_advection_weno.F
      !! Author: Timothy Williams
      !! Date:   20140821, 12:22:17 CEST
      module mod_gridinfo

      implicit none

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!global variables
      !integer, parameter   :: nbdy     = 3
      !integer, parameter   :: margin   = nbdy
      integer, parameter   :: ii = 49
      integer, parameter   :: jj = 51
      real*8, parameter    :: dx = 4e3!m
      real*8, parameter    :: dy = 4e3!m
      !!
      integer, parameter   :: idm   = ii
      integer, parameter   :: jdm   = jj
      !!
      real*8,dimension(idm,jdm)   ::                                    &
     &   X,Y,scuy,scvx,scp2,scp2i,LANDMASK
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module mod_gridinfo
