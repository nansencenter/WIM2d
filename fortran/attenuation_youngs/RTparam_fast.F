      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!Main routine here;
      subroutine RTparam_fast(ac,modT,argR,argT,
     &    alp_nd,hnd,int_adm)

      real, intent(out)  :: ac,modT,argR,argT
      real, intent(out)  :: alp_nd,hnd,int_adm

      real   :: t_a,t_h,modR            !!inputs to RTparam_inner;
      real   :: alp_lin3,alp_lin4
      integer  :: OPT,LOW
      !!
      real   :: alp_nd_lims(6) =
     &   (/1.0e-6,0.005e0,0.3e0,1.5e0,2.5e0,3.5e0/)
      real   :: mc_alplin(2)   =
     &   (/-3.323529252398524e0,3.119943407349375e0/)
      real   :: y0_ll          = 40.0e0
      real   :: dy_ll          = 120.0e0
      integer  :: n_ll           = 3
      real   :: h1_ll          = .4e0
      
      real   :: hnd_lims(3)    = (/1.0e-2,0.2e0,0.4e0/)
      integer  :: LOG_A_vec(5)   = (/1,1,1,0,1/)
      integer  :: LOG_A
      real   :: a0,a1,h0,h1,l0,l1,dtmp,PI
      
      !!start of approximately linear regime;
      alp_lin3       = mc_alplin(2)+mc_alplin(1)*log(hnd)
      alp_nd_lims(5) = alp_lin3
      
      !!end of calculatable results;
      PI       = 4.0*atan(1.0)
      dtmp     = cos(hnd/h1_ll*PI/2.0)
      alp_lin4 = y0_ll+dy_ll*exp(n_ll*log(dtmp))
      
      alp_nd_lims(6) = alp_lin4
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!determine if low/high thickness regime;
      if (hnd<=hnd_lims(1)) then
         hnd    = hnd_lims(1) 
         LOW    = 1
      
      else if (hnd>=hnd_lims(3)) then
         hnd    = hnd_lims(3)
         LOW    = 0
      
      !!now between hnd_lims[0] and hnd_lims[1];
      else if (hnd<hnd_lims(2)) then
         LOW = 1
      else
         LOW = 0
      end if
      
      !!Interpolation value (-1<=t<=1) to use for hnd;
      h0  = hnd_lims(2-LOW)
      h1  = hnd_lims(3-LOW)
      
      if (LOW==0) then
         !!linear interpolation wrt hnd;
         t_h   = -1.0+2*(hnd-h0)/(h1-h0)
      else
         !!interpolation wrt log(hnd);
         !printf("\n(h0,h1)=(%f,%f)\n",h0,h1);
         l0  = log(h0)
         l1  = log(h1)
         t_h = -1.0+2.0*(log(hnd)-l0)/(l1-l0)
      end if
      !printf("\n(t_a,t_h)=(%f,%f)\n",t_a,t_h);
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!determine frequency regime;
      if (alp_nd<=alp_nd_lims(1)) then
         alp_nd = alp_nd_lims(1)
         OPT    = 1
      else if (alp_nd>=alp_nd_lims(6)) then
         alp_nd = alp_nd_lims(6)
         OPT    = 5
      !!
      else if (alp_nd<alp_nd_lims(2)) then
         OPT = 1
      else if (alp_nd<alp_nd_lims(3)) then
         OPT = 2
      else if (alp_nd<alp_nd_lims(4)) then
         OPT = 3
      else if (alp_nd<alp_nd_lims(5)) then
         OPT = 4
      else
         OPT = 5
      end if
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      a0     = alp_nd_lims(OPT)
      a1     = alp_nd_lims(OPT+1)
      LOG_A  = LOG_A_vec(OPT)
      if (LOG_A==1) then
         !!interpolation wrt log(alp_nd);
         l0  = log(a0)
         l1  = log(a1)
         t_a = -1.0+2.0*(log(alp_nd)-l0)/(l1-l0)
      else
         !!interpolation wrt alp_nd;
         t_a = -1.0+2.0*(alp_nd-a0)/(a1-a0)
      end if
        
      !!call main interpolation routine;
      call RTparam_inner(ac,modT,argR,argT,
     &              t_a,t_h,OPT,LOW,int_adm)

      end subroutine RTparam_fast
      
      subroutine RTparam_inner(ac,modT,argR,argT,
     &                  tx,ty,OPT,LOW,int_adm)
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!declare "input/output variables";
      real, intent(out)   :: ac,modT,argR,argT
      real, intent(in)    :: tx,ty,int_adm
      integer, intent(in)   :: OPT,LOW
      !!
      
      !!auxiliary variables;
      real, allocatable  :: chebys(:)!!use allocate/deallocate later
      !real, allocatable  :: Zthings(:)!!use deallocate/deallocate later
      real               :: Zthings(4)!!use deallocate/deallocate later
      !!
      real    :: z,Rr,Ri,Tr,Ti,modR
      integer   :: Ncx,Ncy,ncol
      integer   :: i,j,r
      integer   :: Nthings
      integer   :: INTERP_MODE(5) = (/1,1,3,2,1/)
      integer   :: IM
      
      !!Order of chebyshev expansion for lower thicknesses;
      integer  :: Ncx_L(5) = (/10,10,10,10,3/) !!for alp_nd;
      integer  :: Ncy_L(5) = (/10,10,10,10,10/)!!for hnd;
      
      !!Order of chebyshev expansion for higher thicknesses;
      integer  :: Ncx_H(5) = (/10,10,10,10,4/) !!for alp_nd;
      integer  :: Ncy_H(5) = (/10,10,10,10,10/)!!for hnd;
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      !!Mode of interpolation;
      IM  =INTERP_MODE(OPT)
      if (IM==1) then
         !!Interpolate log(ac), Arg[R,T];
         Nthings   = 3
      else if (IM==2) then
         !!Interpolate (ac), Arg[R,T];
         Nthings   = 3
      else if (IM==3) then
         !!Interpolate Re[R,T], Im[R,T];
         Nthings   = 4
      end if
      
      !print*,'LOW,OPT=',LOW,OPT
      
      if (LOW==1) then
         !!Lower non-dimensional thicknesses;
         Ncx = Ncx_L(OPT)
         Ncy = Ncy_L(OPT)
         call get_Zthings(Zthings,tx,ty,Nthings,Ncx,Ncy,LOW,OPT)
      else
         !!Higher non-dimensional thicknesses;
         Ncx = Ncx_H(OPT)
         Ncy = Ncy_H(OPT)
         call get_Zthings(Zthings,tx,ty,Nthings,Ncx,Ncy,LOW,OPT)
      end if
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!interpret interpolated quantities;
      call RTparam_get_ac(ac,modT,argR,argT,Zthings,IM,int_adm)
      
      end subroutine
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine RTparam_get_ac(ac,modT,argR,argT,
     &      Zthings,IM,int_adm)
      
      real ,intent(out)  :: ac,modT,argR,argT
      real ,intent(in)   :: Zthings(4),int_adm
      integer,intent(in)   :: IM
      !!
      real   ::Rr,Ri,Tr,Ti,modR
      
      if (IM==1) then
         ac     = exp(Zthings(1))
         argR   = Zthings(2)
         argT   = Zthings(3)
         modT   = sqrt(exp(-ac/2.0)/int_adm)
      else if (IM==2) then
         ac     = Zthings(1)
         argR   = Zthings(2)
         argT   = Zthings(3)
         modT   = sqrt(exp(-ac/2.0)/int_adm)
      else if (IM==3) then
         Rr     = Zthings(1)
         Ri     = Zthings(2)
         Tr     = Zthings(3)
         Ti     = Zthings(4)
         !!
         argR   = atan2(Ri,Rr)
         argT   = atan2(Ti,Tr)
         modR   = sqrt(Rr*Rr+Ri*Ri)
         modT   = sqrt(Tr*Tr+Ti*Ti)
         ac     = -2*log(1-modR*modR)
      end if
      
      end subroutine RTparam_get_ac
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
