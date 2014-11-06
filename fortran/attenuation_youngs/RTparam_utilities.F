      !! RTparam_utilities.F
      !! Author: Timothy Williams
      !! Date:   20141009, 10:28:05 CEST
      !! Contains utility functions like chebyshev interpolation (1-3d)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine OP_chebinterp1d(z,tx,chebys,Ncx)

      real,  intent(out) :: z
      real,  intent(in)  :: tx
      integer, intent(in)  :: Ncx
      !!
      real, intent(in) :: chebys(Ncx+1)
      !!
      real   :: an
      real   :: Tn0, Tn1, Tn
      integer  :: s,nx
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      z     = 0.0
      Tn0   = 1.0
      Tn1   = tx
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !! Need to calculate z=\sum_{n=0}^Ncx a_n*T_n(tx),

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! n=0 term;
      nx    = 0
      s     = nx+1
      an    = chebys(s)
      Tn    = Tn0
      z     = z +an*Tn
      !!
      nx    = 1
      s     = nx+1
      an    = chebys(s)
      Tn    = Tn1
      z     = z +an*Tn
      !!
      do nx=2,Ncx
         Tn    = 2.0*tx*Tn1-Tn0
         Tn0   = Tn1
         Tn1   = Tn
         !!
         s  = nx+1
         an = chebys(s)
         z  = z+an*Tn
      end do

      end subroutine OP_chebinterp1d
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine OP_chebinterp2d(z,tx,ty,chebys,Ncx,Ncy)

      real,  intent(out) :: z
      real,  intent(in)  :: tx,ty
      integer, intent(in)  :: Ncx,Ncy
      !!
      real, intent(in)   :: chebys((Ncx+1)*(Ncy+1))
      !!
      real   :: Amn(Ncx+1,Ncy+1)
      real   :: cheb1d(Ncx+1)
      real   :: Tn0, Tn1, Tn
      integer  :: s,n
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Amn   = reshape(chebys,(/ Ncx+1,Ncy+1 /))
      Tn0   = 1.0
      Tn1   = ty
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! make cheb1d, the coefficients for the 1d interpolation
      n        = 0
      Tn       = Tn0
      s        = n+1
      cheb1d   = Amn(:,s)*Tn
      !!
      n        = 1
      Tn       = Tn1
      s        = n+1
      cheb1d   = cheb1d+Amn(:,s)*Tn
      !!
      do n=2,Ncy
         Tn    = 2.0*ty*Tn1-Tn0
         Tn0   = Tn1
         Tn1   = Tn
         !!
         s        = n+1
         cheb1d   = cheb1d+Amn(:,s)*Tn
      end do

      !! now use cheb1d to do the 1d interpolation
      call OP_chebinterp1d(z,tx,cheb1d,Ncx)

      end subroutine OP_chebinterp2d
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine OP_chebinterp3d(z,tx,ty,tz,chebys,Ncx,Ncy,Ncz)

      real,  intent(out) :: z
      real,  intent(in)  :: tx,ty,tz
      integer, intent(in)  :: Ncx,Ncy,Ncz
      !!
      real, intent(in)   :: chebys((Ncx+1)*(Ncy+1)*(Ncz+1))
      !!
      real   :: Amn((Ncx+1)*(Ncy+1),Ncz+1)
      real   :: cheb2d((Ncx+1)*(Ncy+1))
      real   :: Tn0, Tn1, Tn
      integer  :: s,n
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Amn   = reshape(chebys,(/ (Ncx+1)*(Ncy+1),Ncz+1 /))
      Tn0   = 1.0
      Tn1   = tz
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! make cheb2d, the coefficients for the 2d interpolation
      n        = 0
      Tn       = Tn0
      s        = n+1
      cheb2d   = Amn(:,s)*Tn
      !!
      n        = 1
      Tn       = Tn1
      s        = n+1
      cheb2d   = cheb2d+Amn(:,s)*Tn
      !!
      do n=2,Ncz
         Tn    = 2.0*tz*Tn1-Tn0
         Tn0   = Tn1
         Tn1   = Tn
         !!
         s        = n+1
         cheb2d   = cheb2d+Amn(:,s)*Tn
      end do
      
      !! now use cheb2d to do the 2d interpolation
      call OP_chebinterp2d(z,tx,ty,cheb2d,Ncx,Ncy)

      end subroutine OP_chebinterp3d
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
