module gamma_distribution

! Based on John Burkhardt's Fortran 90 code with minimal changes to achieve C/C++ interoperability
! Used here under the terms of the GNU Lesser General Public License


  use,intrinsic :: iso_c_binding
  implicit none


  contains


    real (c_double) function gammad ( x, p, ifault )  bind(c)

    !*****************************************************************************80
    !
    !! GAMMAD computes the Incomplete Gamma Integral
    !
    !  Auxiliary functions:
    !
    !    ALNGAM = logarithm of the gamma function, 
    !    ALNORM = algorithm AS66
    !
    !  Modified:
    !
    !    20 January 2008
    !
    !  Author:
    !
    !    Original FORTRAN77 version by B Shea.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    B Shea,
    !    Algorithm AS 239:
    !    Chi-squared and Incomplete Gamma Integral,
    !    Applied Statistics,
    !    Volume 37, Number 3, 1988, pages 466-473.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) X, P, the parameters of the incomplete 
    !    gamma ratio.  0 <= X, and 0 < P.
    !
    !    Output, integer ( kind = 4 ) IFAULT, error flag.
    !    0, no error.
    !    1, X < 0 or P <= 0.
    !
    !    Output, real ( kind = 8 ) GAMMAD, the value of the incomplete 
    !    Gamma integral.
    !
      use,intrinsic :: iso_c_binding
      implicit none

      real  (c_double), intent(in)  :: x
      real  (c_double), intent(in)  :: p      
      integer (c_int), intent(out) :: ifault
      
      real ( c_double ) a
!      real ( kind = 8 ) alnorm
!      real ( kind = 8 ) alngam
      real ( c_double ) an
      real ( c_double ) arg
      real ( c_double ) b
      real ( c_double ) c
      real ( c_double ), parameter :: elimit = - 88.0D+00
!rlk      real ( kind = 8 ) gammad
!rlk      integer ( kind = 4 ) ifault
      real ( c_double ), parameter :: oflo = 1.0D+37
!rlk      real ( kind = 8 ) p
      real ( c_double ), parameter :: plimit = 1000.0D+00
      real ( c_double ) pn1
      real ( c_double ) pn2
      real ( c_double ) pn3
      real ( c_double ) pn4
      real ( c_double ) pn5
      real ( c_double ) pn6
      real ( c_double ) rn
      real ( c_double ), parameter :: tol = 1.0D-14
      logical (c_bool) upper
!rlk      real ( kind = 8 ) x
      real ( c_double ), parameter :: xbig = 1.0D+08

      gammad = 0.0D+00
    !
    !  Check the input.
    !
      if ( x < 0.0D+00 ) then
        ifault = 1
        return
      end if

      if ( p <= 0.0D+00 ) then
        ifault = 1
        return
      end if

      ifault = 0

      if ( x == 0.0D+00 ) then
        gammad = 0.0D+00
        return
      end if
    !
    !  If P is large, use a normal approximation.
    !
      if ( plimit < p ) then

        pn1 = 3.0D+00 * sqrt ( p ) * ( ( x / p )**( 1.0D+00 / 3.0D+00 ) &
        + 1.0D+00 / ( 9.0D+00 * p ) - 1.0D+00 )

        upper = .false.
        gammad = alnorm ( pn1, upper )
        return

      end if
    !
    !  If X is large set GAMMAD = 1.
    !
      if ( xbig < x ) then
        gammad = 1.0D+00
        return
      end if
    !
    !  Use Pearson's series expansion.
    !  (Note that P is not large enough to force overflow in ALOGAM).
    !  No need to test IFAULT on exit since P > 0.
    !
      if ( x <= 1.0D+00 .or. x < p ) then

        arg = p * log ( x ) - x - alngam ( p + 1.0D+00, ifault )
        c = 1.0D+00
        gammad = 1.0D+00
        a = p

        do

          a = a + 1.0D+00
          c = c * x / a
          gammad = gammad + c

          if ( c <= tol ) then
            exit
          end if

        end do

        arg = arg + log ( gammad )

        if ( elimit <= arg ) then
          gammad = exp ( arg )
        else
          gammad = 0.0D+00
        end if
    !
    !  Use a continued fraction expansion.
    !
      else 

        arg = p * log ( x ) - x - alngam ( p, ifault )
        a = 1.0D+00 - p
        b = a + x + 1.0D+00
        c = 0.0D+00
        pn1 = 1.0D+00
        pn2 = x
        pn3 = x + 1.0D+00
        pn4 = x * b
        gammad = pn3 / pn4

        do

          a = a + 1.0D+00
          b = b + 2.0D+00
          c = c + 1.0D+00
          an = a * c
          pn5 = b * pn3 - an * pn1
          pn6 = b * pn4 - an * pn2

          if ( pn6 /= 0.0D+00 ) then

            rn = pn5 / pn6

            if ( abs ( gammad - rn ) <= min ( tol, tol * rn ) ) then
              exit
            end if

            gammad = rn

          end if

          pn1 = pn3
          pn2 = pn4
          pn3 = pn5
          pn4 = pn6
    !
    !  Re-scale terms in continued fraction if terms are large.
    !
          if ( oflo <= abs ( pn5 ) ) then
            pn1 = pn1 / oflo
            pn2 = pn2 / oflo
            pn3 = pn3 / oflo
            pn4 = pn4 / oflo
          end if

        end do

        arg = arg + log ( gammad )

        if ( elimit <= arg ) then
          gammad = 1.0D+00 - exp ( arg )
        else
          gammad = 1.0D+00
        end if

      end if

      return
    end function gammad
    
    real (c_double) function alngam ( xvalue, ifault ) bind(c)

    !*****************************************************************************80
    !
    !! ALNGAM computes the logarithm of the gamma function.
    !
    !  Modified:
    !
    !    13 January 2008
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Allan Macleod.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Allan Macleod,
    !    Algorithm AS 245,
    !    A Robust and Reliable Algorithm for the Logarithm of the Gamma Function,
    !    Applied Statistics,
    !    Volume 38, Number 2, 1989, pages 397-402.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) XVALUE, the argument of the Gamma function.
    !
    !    Output, integer ( kind = 4 ) IFAULT, error flag.
    !    0, no error occurred.
    !    1, XVALUE is less than or equal to 0.
    !    2, XVALUE is too big.
    !
    !    Output, real ( kind = 8 ) ALNGAM, the logarithm of the gamma function of X.
    !
    
      use,intrinsic :: iso_c_binding
      implicit none

      real  (c_double), intent(in)  :: xvalue  
      integer  (c_int), intent(out) :: ifault

!rlk      real ( kind = 8 ) alngam
      real ( c_double ), parameter :: alr2pi = 0.918938533204673D+00
!rlk      integer ( kind = 4 ) ifault
      real ( c_double ), dimension ( 9 ) :: r1 = (/ &
        -2.66685511495D+00, &
        -24.4387534237D+00, &
        -21.9698958928D+00, &
         11.1667541262D+00, &
         3.13060547623D+00, &
         0.607771387771D+00, &
         11.9400905721D+00, &
         31.4690115749D+00, &
         15.2346874070D+00 /)
      real ( c_double ), dimension ( 9 ) :: r2 = (/ &
        -78.3359299449D+00, &
        -142.046296688D+00, &
         137.519416416D+00, &
         78.6994924154D+00, &
         4.16438922228D+00, &
         47.0668766060D+00, &
         313.399215894D+00, &
         263.505074721D+00, &
         43.3400022514D+00 /)
      real ( c_double ), dimension ( 9 ) :: r3 = (/ &
        -2.12159572323D+05, &
         2.30661510616D+05, &
         2.74647644705D+04, &
        -4.02621119975D+04, &
        -2.29660729780D+03, &
        -1.16328495004D+05, &
        -1.46025937511D+05, &
        -2.42357409629D+04, &
        -5.70691009324D+02 /)
      real ( c_double ), dimension ( 5 ) :: r4 = (/ &
         0.279195317918525D+00, &
         0.4917317610505968D+00, &
         0.0692910599291889D+00, &
         3.350343815022304D+00, &
         6.012459259764103D+00 /)
      real ( c_double ) x
      real ( c_double ) x1
      real ( c_double ) x2
      real ( c_double ), parameter :: xlge = 5.10D+05
      real ( c_double ), parameter :: xlgst = 1.0D+30
!rlk      real ( kind = 8 ) xvalue
      real ( c_double ) y

      x = xvalue
      alngam = 0.0D+00
    !
    !  Check the input.
    !
      if ( xlgst <= x ) then
        ifault = 2
        return
      end if

      if ( x <= 0.0D+00 ) then
        ifault = 1
        return
      end if

      ifault = 0
    !
    !  Calculation for 0 < X < 0.5 and 0.5 <= X < 1.5 combined.
    !
      if ( x < 1.5D+00 ) then

        if ( x < 0.5D+00 ) then

          alngam = - log ( x )
          y = x + 1.0D+00
    !
    !  Test whether X < machine epsilon.
    !
          if ( y == 1.0D+00 ) then
            return
          end if

        else

          alngam = 0.0D+00
          y = x
          x = ( x - 0.5D+00 ) - 0.5D+00

        end if

        alngam = alngam + x * (((( &
            r1(5)   * y &
          + r1(4) ) * y &
          + r1(3) ) * y &
          + r1(2) ) * y &
          + r1(1) ) / (((( &
                      y &
          + r1(9) ) * y &
          + r1(8) ) * y &
          + r1(7) ) * y &
          + r1(6) )

        return

      end if
    !
    !  Calculation for 1.5 <= X < 4.0.
    !
      if ( x < 4.0D+00 ) then

        y = ( x - 1.0D+00 ) - 1.0D+00

        alngam = y * (((( &
            r2(5)   * x &
          + r2(4) ) * x &
          + r2(3) ) * x &
          + r2(2) ) * x &
          + r2(1) ) / (((( &
                      x &
          + r2(9) ) * x &
          + r2(8) ) * x &
          + r2(7) ) * x &
          + r2(6) )
    !
    !  Calculation for 4.0 <= X < 12.0.
    !
      else if ( x < 12.0D+00 ) then

        alngam = (((( &
            r3(5)   * x &
          + r3(4) ) * x &
          + r3(3) ) * x &
          + r3(2) ) * x &
          + r3(1) ) / (((( &
                      x &
          + r3(9) ) * x &
          + r3(8) ) * x &
          + r3(7) ) * x &
          + r3(6) )
    !
    !  Calculation for 12.0 <= X.
    !
      else

        y = log ( x )
        alngam = x * ( y - 1.0D+00 ) - 0.5D+00 * y + alr2pi

        if ( x <= xlge ) then

          x1 = 1.0D+00 / x
          x2 = x1 * x1

          alngam = alngam + x1 * ( ( &
                 r4(3)   * &
            x2 + r4(2) ) * &
            x2 + r4(1) ) / ( ( &
            x2 + r4(5) ) * &
            x2 + r4(4) )

        end if

      end if

      return
    end function alngam
    
       real (c_double) function alnorm ( x, upper ) bind(c)

    !*****************************************************************************80
    !
    !! ALNORM computes the cumulative density of the standard normal distribution.
    !
    !  Modified:
    !
    !    13 January 2008
    !
    !  Author:
    !
    !    Original FORTRAN77 version by David Hill.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    David Hill,
    !    Algorithm AS 66:
    !    The Normal Integral,
    !    Applied Statistics,
    !    Volume 22, Number 3, 1973, pages 424-427.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) X, is one endpoint of the semi-infinite interval
    !    over which the integration takes place.
    !
    !    Input, logical UPPER, determines whether the upper or lower
    !    interval is to be integrated:
    !    .TRUE.  => integrate from X to + Infinity;
    !    .FALSE. => integrate from - Infinity to X.
    !
    !    Output, real ( kind = 8 ) ALNORM, the integral of the standard normal
    !    distribution over the desired interval.
    !
    
      use,intrinsic :: iso_c_binding
      implicit none
      
      real    (c_double), intent(in)  :: x
      logical (c_bool),   intent(in)  :: upper
       

      
      real ( c_double ), parameter :: a1 = 5.75885480458D+00
      real ( c_double ), parameter :: a2 = 2.62433121679D+00
      real ( c_double ), parameter :: a3 = 5.92885724438D+00
!rlk      real ( kind = 8 ) alnorm
      real ( c_double ), parameter :: b1 = -29.8213557807D+00
      real ( c_double ), parameter :: b2 = 48.6959930692D+00
      real ( c_double ), parameter :: c1 = -0.000000038052D+00
      real ( c_double ), parameter :: c2 = 0.000398064794D+00
      real ( c_double ), parameter :: c3 = -0.151679116635D+00
      real ( c_double ), parameter :: c4 = 4.8385912808D+00
      real ( c_double ), parameter :: c5 = 0.742380924027D+00
      real ( c_double ), parameter :: c6 = 3.99019417011D+00
      real ( c_double ), parameter :: con = 1.28D+00
      real ( c_double ), parameter :: d1 = 1.00000615302D+00
      real ( c_double ), parameter :: d2 = 1.98615381364D+00
      real ( c_double ), parameter :: d3 = 5.29330324926D+00
      real ( c_double ), parameter :: d4 = -15.1508972451D+00
      real ( c_double ), parameter :: d5 = 30.789933034D+00
      real ( c_double ), parameter :: ltone = 7.0D+00
      real ( c_double ), parameter :: p = 0.398942280444D+00
      real ( c_double ), parameter :: q = 0.39990348504D+00
      real ( c_double ), parameter :: r = 0.398942280385D+00
      logical up
!rlk      logical upper
      real ( c_double ), parameter :: utzero = 18.66D+00
!rlk      real ( kind = 8 ) x
      real ( c_double ) y
      real ( c_double ) z

      up = upper
      z = x

      if ( z < 0.0D+00 ) then
        up = .not. up
        z = - z
      end if

      if ( ltone < z .and. ( ( .not. up ) .or. utzero < z ) ) then

        if ( up ) then
          alnorm = 0.0D+00
        else
          alnorm = 1.0D+00
        end if

        return

      end if

      y = 0.5D+00 * z * z

      if ( z <= con ) then

        alnorm = 0.5D+00 - z * ( p - q * y &
          / ( y + a1 + b1 &
          / ( y + a2 + b2 & 
          / ( y + a3 ))))

      else

        alnorm = r * exp ( - y ) &
          / ( z + c1 + d1 &
          / ( z + c2 + d2 &
          / ( z + c3 + d3 &
          / ( z + c4 + d4 &
          / ( z + c5 + d5 &
          / ( z + c6 ))))))

      end if

      if ( .not. up ) then
        alnorm = 1.0D+00 - alnorm
      end if

      return
    end function alnorm
    
end module gamma_distribution
