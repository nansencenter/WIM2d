      module mod_common_wim
         !!global variables: ice conc/thickness/floe size

      use mod_gridinfo

      implicit none

      real,dimension(ii,jj)   :: icec,iceh,dfloe,ICE_MASK
      real,dimension(ii,jj)   :: WTR_MASK
      real,dimension(idm,jdm) :: Hs,Tp,mwd,WAVE_MASK  !set in waves_init
      !!
      real,parameter          :: unif_c            = .7
      real,parameter          :: unif_h            = 1.5
      real,parameter          :: dfloe_pack_init   = 300.0

      end module mod_common_wim
