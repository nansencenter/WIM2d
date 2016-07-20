      module mod_interp
      
      implicit none

      private
      
      real ( kind = 8 ), parameter :: r8_huge = 1.79769313486231571D+308 !from r8lib.f90
      !integer ( kind = 4 )  :: r8vec_bracket5 !from r8lib.f90
      
      public   :: pwl_interp_2d,pwl_interp_2d_real
      
      contains
      
      subroutine pwl_interp_2d_real( nxd, nyd, xd, yd, zd, ni, xi, yi, zi )
      
      integer,intent(in)         :: nxd,nyd,ni
      real(kind=4),intent(in)    :: xd(nxd),yd(nyd),zd(nxd,nyd),xi(ni),yi(ni)
      real(kind=4),intent(out)   :: zi(ni)
      real(kind=8)   :: zi8(ni)
      real(kind=8)   :: xi8(ni),yi8(ni),xd8(nxd),yd8(nyd),zd8(nxd,nyd)
      
      xd8   = xd
      yd8   = yd
      zd8   = zd
      xi8   = xi
      yi8   = yi
      call pwl_interp_2d ( nxd, nyd, xd8, yd8, zd8,      &
     &         ni, xi8, yi8, zi8 )
!     call pwl_interp_2d ( nxd, nyd, double(xd), double(yd), double(zd),      &
!    &         ni, double(xi), double(yi), zi8 )
      zi = zi8

      end subroutine pwl_interp_2d_real
      
      
      subroutine pwl_interp_2d ( nxd, nyd, xd, yd, zd, ni, xi, yi, zi )
      
      !*****************************************************************************80
      !
      !! PWL_INTERP_2D: piecewise linear interpolant to data defined on a 2D grid.
      !
      !  Licensing:
      !
      !    This code is distributed under the GNU LGPL license.
      !
      !  Modified:
      !
      !    14 October 2012
      !
      !  Author:
      !
      !    John Burkardt
      !
      !  Parameters:
      !
      !    Input, integer ( kind = 4 ) NXD, NYD, the number of X and Y data values.
      !
      !    Input, real ( kind = 8 ) XD(NXD), YD(NYD), the sorted X and Y data.
      !
      !    Input, real ( kind = 8 ) ZD(NXD,NYD), the Z data.
      !
      !    Input, integer ( kind = 4 ) NI, the number of interpolation points.
      !
      !    Input, real ( kind = 8 ) XI(NI), YI(NI), the coordinates of the
      !    interpolation points.
      !
      !    Output, real ( kind = 8 ) ZI(NI), the value of the interpolant.
      !
        implicit none
      
        integer ( kind = 4 ) ni
        integer ( kind = 4 ) nxd
        integer ( kind = 4 ) nyd
      
        real ( kind = 8 ) alpha
        real ( kind = 8 ) beta
        real ( kind = 8 ) det
        real ( kind = 8 ) dxa
        real ( kind = 8 ) dxb
        real ( kind = 8 ) dxi
        real ( kind = 8 ) dya
        real ( kind = 8 ) dyb
        real ( kind = 8 ) dyi
        real ( kind = 8 ) gamma
        integer ( kind = 4 ) i
        integer ( kind = 4 ) j
        integer ( kind = 4 ) k
        !real ( kind = 8 ) r8_huge
        !integer ( kind = 4 ) r8vec_bracket5
        real ( kind = 8 ) xd(nxd)
        real ( kind = 8 ) xi(ni)
        real ( kind = 8 ) yd(nyd)
        real ( kind = 8 ) yi(ni)
        real ( kind = 8 ) zd(nxd,nyd)
        real ( kind = 8 ) zi(ni)
      
        do k = 1, ni
      
          i = r8vec_bracket5 ( nxd, xd, xi(k) )
          if ( i == -1 ) then
            zi(k) = r8_huge !( )
            cycle
          end if
      
          j = r8vec_bracket5 ( nyd, yd, yi(k) )
          if ( j == -1 ) then
            zi(k) = r8_huge !( )
            cycle
          end if
      
          if ( yi(k) < yd(j+1) &
            + ( yd(j) - yd(j+1) ) * ( xi(i) - xd(i) ) / ( xd(i+1) - xd(i) ) ) then
      
            dxa = xd(i+1) - xd(i)
            dya = yd(j)   - yd(j)
      
            dxb = xd(i)   - xd(i)
            dyb = yd(j+1) - yd(j)
      
            dxi = xi(k)   - xd(i)
            dyi = yi(k)   - yd(j)
      
            det = dxa * dyb - dya * dxb
      
            alpha = ( dxi * dyb - dyi * dxb ) / det
            beta =  ( dxa * dyi - dya * dxi ) / det
            gamma = 1.0D+00 - alpha - beta
      
            zi(k) = alpha * zd(i+1,j) + beta * zd(i,j+1) + gamma * zd(i,j)
      
          else
      
            dxa = xd(i)   - xd(i+1)
            dya = yd(j+1) - yd(j+1)
      
            dxb = xd(i+1) - xd(i+1)
            dyb = yd(j)   - yd(j+1)
      
            dxi = xi(k)   - xd(i+1)
            dyi = yi(k)   - yd(j+1)
      
            det = dxa * dyb - dya * dxb
      
            alpha = ( dxi * dyb - dyi * dxb ) / det
            beta =  ( dxa * dyi - dya * dxi ) / det
            gamma = 1.0D+00 - alpha - beta
      
            zi(k) = alpha * zd(i,j+1) + beta * zd(i+1,j) + gamma * zd(i+1,j+1)
      
          end if
      
        end do
      
        return
      end
      
      function r8vec_bracket5 ( nd, xd, xi )

      !*****************************************************************************80
      !
      !! R8VEC_BRACKET5 brackets data between successive entries of a sorted R8VEC.
      !
      !  Discussion:
      !
      !    We assume XD is sorted.
      !
      !    If XI is contained in the interval [XD(1),XD(N)], then the returned 
      !    value B indicates that XI is contained in [ XD(B), XD(B+1) ].
      !
      !    If XI is not contained in the interval [XD(1),XD(N)], then B = -1.
      !
      !    This code implements a version of binary search which is perhaps more
      !    understandable than the usual ones.
      !
      !  Licensing:
      !
      !    This code is distributed under the GNU LGPL license.
      !
      !  Modified:
      !
      !    14 October 2012
      !
      !  Author:
      !
      !    John Burkardt
      !
      !  Parameters:
      !
      !    Input, integer ( kind = 4 ) ND, the number of data values.
      !
      !    Input, real ( kind = 8 ) XD(N), the sorted data.
      !
      !    Input, real ( kind = 8 ) XD, the query value.
      !
      !    Output, integer ( kind = 4 ) R8VEC_BRACKET5, the bracket information.
      !
        implicit none

        integer ( kind = 4 ) nd

        integer ( kind = 4 ) b
        integer ( kind = 4 ) l
        integer ( kind = 4 ) m
        integer ( kind = 4 ) r
        integer ( kind = 4 ) r8vec_bracket5
        real ( kind = 8 ) xd(nd)
        real ( kind = 8 ) xi

        if ( xi < xd(1) .or. xd(nd) < xi ) then

          b = -1

        else

          l = 1
          r = nd

          do while ( l + 1 < r )
            m = ( l + r ) / 2
            if ( xi < xd(m) ) then
              r = m
            else
              l = m
            end if
          end do

          b = l

        end if

        r8vec_bracket5 = b

        return
      end
      end module mod_interp
