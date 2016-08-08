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
     &                  ,real_prams,nx,ny,nr)

      use mod_wim2d_run, only: wim2d_run_io

      implicit none

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!variables
      integer,intent(in)                  :: nx,ny,nr
      real*8,intent(in)                   :: real_prams(nr)
      real*8,dimension(nx,ny,6),intent(in)  :: in_arrays!!c,h,Dmax,Hs,Tp,mwd
      real*8,dimension(nx,ny,6),intent(out) :: out_arrays!!Dmax,Hs,Tp,taux,tauy,mwd
!f2py intent(hide)   :: nx,ny,nr

      call wim2d_run_io(in_arrays,out_arrays,real_prams,nx,ny)

      end subroutine py_wim2d_run_io
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine py_wim2d_run_vSdir(Sdir_out,out_arrays                 &
     &                 ,Sdir_in,ice_arrays                              &
     &                 ,real_prams,ndir,nw,nx,ny,nr)

      use mod_common_wim,  only: icec,iceh,dfloe                        &
     &                           ,sdf_dir,sdf_inc                       &
     &                           ,freq_vec,wavdir
      use mod_wim2d_run,   only: wim2d_run_vSdir,read_param_vec,STEADY
      use mod_wavesice,    only: WAVE_MASK2

      implicit none

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! array dimensions (detected by python but need to declare them)
      integer,intent(in)   :: nx,ny,ndir,nw,nr
!f2py intent(hide)   :: nx,ny,nr

      !!arrays
      real*8,dimension(nx,ny,3),intent(in)        :: ice_arrays  !! pass in  icec,iceh,dfloe
      real*8,dimension(nx,ny,6),intent(out)       :: out_arrays  !! pass out Dmax,Hs,Tp,taux,tauy,mwd
      real*8,dimension(nx*ny*ndir*nw),intent(in)  :: Sdir_in     !! pass in/out directional spec
      real*8,dimension(nx*ny*ndir*nw),intent(out) :: Sdir_out    !! pass in/out directional spec

      !! parameters
      real*8,intent(in) :: real_prams(nr)
      !!
      integer  :: i,j,w,wth
      real     :: PI,adv_dir
      real     :: Tp_io,mwd_io,init_cons(6)

      PI       = 4.*atan(1.)
      call read_param_vec(init_cons,real_prams)
      Tp_io    = init_cons(2)
      mwd_io   = init_cons(3)

      !! initial ice conditions
      icec     = real(ice_arrays(:,:,1))
      iceh     = real(ice_arrays(:,:,2))
      dfloe    = real(ice_arrays(:,:,3))

      !! get input directional wave spectrum
      sdf_dir  = reshape(real(Sdir_in),(/nx,ny,ndir,nw/))
      sdf_inc  = sdf_dir

      print*,nw,(nw.eq.1),Tp_io
      if (nw.eq.1) then
         freq_vec(1) = 1./Tp_io
         print*,Tp_io,freq_vec
      end if
      if (ndir.eq.1) then
         wavdir(1)   = mwd_io
         print*,mwd_io
      end if

      !wim2d_run_vSdir(out_arrays,int_prams,real_prams,nmesh_e,mesh_e)
      ! mesh_e is optional, but need a place holder for nmesh_e
      call wim2d_run_vSdir(out_arrays,1)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! get updated directional wave spectrum

      !! correction for STEADY (copied from mod_waves_ice.F)
      if (STEADY.eq.1) then
         do w=1,nw
         do wth=1,ndir
            adv_dir  = -PI/180.0*(wavdir(wth)+90.0)
            if (cos(adv_dir).ge.0.0) then
               where (WAVE_MASK2.gt.0.0)
                   sdf_dir(:,:,wth,w)   = sdf_inc(:,:,wth,w)
               end where
            end if
         end do
         end do
      end if

      !!TODO update Hs,Tp,... also?
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      Sdir_out = reshape(dble(sdf_dir),(/nx*ny*ndir*nw/))

      end subroutine py_wim2d_run_vSdir
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine py_wim2d_run_vSdir_mesh(Sdir_out,out_arrays            &
     &                  ,Sdir_in,ice_arrays,mesh_in,mesh_out            &
     &                  ,real_prams,ndir,nw,nmesh_e,nx,ny,nr)

      use mod_common_wim,  only: icec,iceh,dfloe                        &
     &                           ,sdf_dir,sdf_inc                       &
     &                           ,freq_vec,wavdir
      use mod_wim2d_run,   only: wim2d_run_vSdir,read_param_vec,STEADY
      use mod_wavesice,    only: WAVE_MASK2

      implicit none

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! array dimensions (detected by python but need to declare them)
      integer,intent(in)   :: nx,ny,ndir,nw,nmesh_e,nr
!f2py intent(hide)   :: nx,ny,nr

      !!arrays
      real*8,dimension(nx,ny,3),intent(in)        :: ice_arrays   !! pass in  icec,iceh,dfloe
      real*8,dimension(nx,ny,6),intent(out)       :: out_arrays   !! pass out Dmax,Hs,Tp,taux,tauy,mwd
      real*8,dimension(nx*ny*ndir*nw),intent(in)  :: Sdir_in      !! pass in directional spec
      real*8,dimension(nx*ny*ndir*nw),intent(out) :: Sdir_out     !! pass out directional spec

      integer,parameter    :: nmesh_vars = 6
      real*8,dimension(nmesh_e*nmesh_vars),intent(in)
     &   :: mesh_in                                        !! pass in mesh_in
      real*8,dimension(nmesh_e*nmesh_vars),intent(out)
     &   :: mesh_out                                       !! pass out mesh_out
      real*8,dimension(nmesh_e,nmesh_vars)  :: mesh_e_dummy  !! reshaped mesh_e

      !! parameters
      real*8,intent(in) :: real_prams(nr)
      !!
      integer  :: i,j,w,wth
      real     :: PI,adv_dir
      real     :: Tp_io,mwd_io,init_cons(6)

      PI       = 4.*atan(1.)
      call read_param_vec(init_cons,real_prams)
      Tp_io    = init_cons(2)
      mwd_io   = init_cons(3)

      !! initial ice conditions
      icec     = ice_arrays(:,:,1)
      iceh     = ice_arrays(:,:,2)
      dfloe    = ice_arrays(:,:,3)

      !! get input directional wave spectrum
      sdf_dir  = reshape(real(Sdir_in),(/nx,ny,ndir,nw/))
      sdf_inc  = sdf_dir

      print*,nw,(nw.eq.1),Tp_io
      if (nw.eq.1) then
         freq_vec(1) = 1./Tp_io
         print*,Tp_io,freq_vec
      end if
      if (ndir.eq.1) then
         wavdir(1)   = mwd_io
         print*,mwd_io
      end if

      !wim2d_run_vSdir(out_arrays,int_prams,real_prams,nmesh_e,mesh_e)
      ! mesh_e is optional, but need a place holder for nmesh_e
      mesh_e_dummy   = reshape(mesh_in,(/nmesh_e,nmesh_vars/))
      call wim2d_run_vSdir(out_arrays
     &                     ,nmesh_e,mesh_e_dummy)
      mesh_out       = reshape(mesh_e_dummy,(/nmesh_e*nmesh_vars/))

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! get updated directional wave spectrum

      !! correction for STEADY (copied from mod_waves_ice.F)
      if (STEADY.eq.1) then
         do w=1,nw
         do wth=1,ndir
            adv_dir  = -PI/180.0*(wavdir(wth)+90.0)
            if (cos(adv_dir).ge.0.0) then
               where (WAVE_MASK2.gt.0.0)
                   sdf_dir(:,:,wth,w)   = sdf_inc(:,:,wth,w)
               end where
            end if
         end do
         end do
      end if

      !!TODO update Hs,Tp,... also?
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      Sdir_out = reshape(dble(sdf_dir),(/nx*ny*ndir*nw/))


      end subroutine py_wim2d_run_vSdir_mesh
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine wim2d_dimensions(nx,ny,nw,ndir)

      implicit none

      include "grid_info.h"
      include "wave_info.h"

      integer, intent(out) :: nx,ny,nw,ndir

      nx    = ii
      ny    = jj
      nw    = n_wave_freq
      ndir  = n_wavdir

      end subroutine wim2d_dimensions
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine wim2d_freq_vec(nw,freq_vec)

      implicit none

      include "wave_info.h"

      integer, intent(in)  :: nw
      real,    intent(out) :: freq_vec(nw)

      real     :: freq,fmin,fmax,df
      integer  :: w

      freq_vec = -1.
      if (nw.gt.1) then

         fmin     = 1./Tmax
         fmax     = 1./Tmin
         df       = (fmax-fmin)/real(nw-1)
         freq     = fmin

         do w=1,nw
            freq_vec(w) = freq
            freq        = freq+df
         end do
                  
      end if

      end subroutine wim2d_freq_vec
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine atten_youngs(outputs,inputs)
      use mod_common_wim, only: wim_info

      real,intent(in)   :: inputs(4)
      real,intent(out)  :: outputs(8)

      real  :: damping,kice,kwtr,int_adm,ac,modT,argR,argT

      real  :: guess,om

      ! h        = inputs(1)
      ! om       = inputs(2)
      wim_info%young    = inputs(3)
      wim_info%viscrp   = inputs(4)
      om    = inputs(2)
      guess = om**2/wim_info%gravity

      call RTparam_outer(damping,kice,kwtr,int_adm,
     &                  ac,modT,argR,argT,
     &                  inputs(1:2),guess)

      outputs  = (/damping,kice,kwtr,int_adm,
     &                  ac,modT,argR,argT/)

      end subroutine atten_youngs
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
