# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "circle_rule.hpp"

//****************************************************************************80

void half_circle_rule ( int nt, long double w[], long double t[] )

//****************************************************************************80
//
//    The integral I(f) is approximated by
//
//      Q(f) = pi * sum ( 1 <= i <= NT ) W(i) * F ( cos(T(i)), sin(T(i)) ).
//
//  Parameters:
//
//    Input, int NT, the number of angles to use.
//
//    Output, long double W[NT], the weights for the rule.
//
//    Output, long double T[NT], the angles for the rule.
//
{
  int it;
  long double r8_pi = 3.141592653589793238462643383279502884;

  for ( it = 0; it < nt; it++ )
  {
      if (it == 0 || it == nt-1) {
          w[it] = 0.5 / ( long double ) ( nt - 1 );
      } else {
          w[it] = 1.0 / ( long double ) ( nt - 1 );
      }
      t[it] = r8_pi * ( long double ) ( it ) / ( long double ) ( nt - 1 );
    
      // w[it] = 1.0 / ( long double ) ( nt );
    // t[it] = r8_pi * ( long double ) ( 2.0 * it + 1.0 ) / ( long double ) ( 2.0 * nt );
  }
  return;
}
//****************************************************************************80

//****************************************************************************80

void circle_rule ( int nt, long double x[], long double y[], long double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    CIRCLE_RULE computes a quadrature rule for the unit circle.
//
//  Discussion:
//
//    The unit circle is the region:
//
//      x * x + y * y = 1.
//
//    The integral I(f) is then approximated by
//
//      Q(f) = 2 * pi * sum ( 1 <= i <= NT ) W(i) * F ( cos(T(i)), sin(T(i)) ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 April 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NT, the number of angles to use.
//
//    Output, long double W[NT], the weights for the rule.
//
//    Output, long double T[NT], the angles for the rule.
//
{
  int it;
  long double t;
  long double r8_pi = 3.141592653589793238462643383279502884;

  for ( it = 0; it < nt; it++ )
  {
    w[it] = 1.0 / ( long double ) ( nt );
    t = 2.0 * r8_pi * ( long double ) ( it ) / ( long double ) ( nt );
    x[it] = cos(t);
    y[it] = sin(t);
  }
  return;
}

//****************************************************************************80

void circle_rule ( int nt, long double w[], long double t[] )

//****************************************************************************80
//
//  Purpose:
//
//    CIRCLE_RULE computes a quadrature rule for the unit circle.
//
//  Discussion:
//
//    The unit circle is the region:
//
//      x * x + y * y = 1.
//
//    The integral I(f) is then approximated by
//
//      Q(f) = 2 * pi * sum ( 1 <= i <= NT ) W(i) * F ( cos(T(i)), sin(T(i)) ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 April 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NT, the number of angles to use.
//
//    Output, long double W[NT], the weights for the rule.
//
//    Output, long double T[NT], the angles for the rule.
//
{
  int it;
  long double r8_pi = 3.141592653589793238462643383279502884;

  for ( it = 0; it < nt; it++ )
  {
    w[it] = 1.0 / ( long double ) ( nt );
    t[it] = 2.0 * r8_pi * ( long double ) ( it ) / ( long double ) ( nt );
  }
  return;
}
//****************************************************************************80

long double circle01_monomial_integral ( int e[2] )

//****************************************************************************80
//
//  Purpose:
//
//    CIRCLE01_MONOMIAL_INTEGRAL returns monomial integrals on the unit circle.
//
//  Discussion:
//
//    The integration region is 
//
//      X^2 + Y^2 = 1.
//
//    The monomial is F(X,Y) = X^E(1) * Y^E(2).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Philip Davis, Philip Rabinowitz,
//    Methods of Numerical Integration,
//    Second Edition,
//    Academic Press, 1984, page 263.
//
//  Parameters:
//
//    Input, int E[2], the exponents of X and Y in the 
//    monomial.  Each exponent must be nonnegative.
//
//    Output, long double CIRCLE01_MONOMIAL_INTEGRAL, the integral.
//
{
  long double arg;
  int i;
  long double integral;

  if ( e[0] < 0 || e[1] < 0 )
  {
    cout << "\n";
    cout << "CIRCLE01_MONOMIAL_INTEGRAL - Fatal error!\n";
    cout << "  All exponents must be nonnegative.\n";
    cout << "  E[0] = " << e[0] << "\n";
    cout << "  E[1] = " << e[1] << "\n";
    exit ( 1 );
  }

  if ( ( e[0] % 2 ) == 1 || ( e[1] % 2 ) == 1 )
  {
    integral = 0.0;
  }
  else
  {
    integral = 2.0;

    for ( i = 0; i < 2; i++ )
    {
      arg = 0.5 * ( long double ) ( e[i] + 1 );
      integral = integral * r8_gamma_circle ( arg );
    }

    arg = 0.5 * ( long double ) ( e[0] + e[1] + 2 );
    integral = integral / r8_gamma_circle ( arg );
  }
  return integral;
}
//****************************************************************************80

long double r8_gamma_circle ( long double x )

//****************************************************************************80
//
//  Purpose:
//
//    r8_gamma_circle evaluates Gamma(X) for a real argument.
//
//  Discussion:
//
//    The C MATH library includes a function GAMMA ( X ) which should be
//    invoked instead of this function.
//
//    This routine calculates the gamma function for a real argument X.
//
//    Computation is based on an algorithm outlined in reference 1.
//    The program uses rational functions that approximate the gamma
//    function to at least 20 significant decimal digits.  Coefficients
//    for the approximation over the interval (1,2) are unpublished.
//    Those for the approximation for 12 <= X are from reference 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 January 2008
//
//  Author:
//
//    Original FORTRAN77 version by William Cody, Laura Stoltz.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    William Cody,
//    An Overview of Software Development for Special Functions,
//    in Numerical Analysis Dundee, 1975,
//    edited by GA Watson,
//    Lecture Notes in Mathematics 506,
//    Springer, 1976.
//
//    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
//    Charles Mesztenyi, John Rice, Henry Thatcher,
//    Christoph Witzgall,
//    Computer Approximations,
//    Wiley, 1968,
//    LC: QA297.C64.
//
//  Parameters:
//
//    Input, long double X, the argument of the function.
//
//    Output, long double r8_gamma_circle, the value of the function.
//
{
  long double c[7] = {
   -1.910444077728E-03, 
    8.4171387781295E-04, 
   -5.952379913043012E-04, 
    7.93650793500350248E-04, 
   -2.777777777777681622553E-03, 
    8.333333333333333331554247E-02, 
    5.7083835261E-03 };
  long double eps = 2.22E-16;
  long double fact;
  int i;
  int n;
  long double p[8] = {
  -1.71618513886549492533811E+00,
   2.47656508055759199108314E+01, 
  -3.79804256470945635097577E+02,
   6.29331155312818442661052E+02, 
   8.66966202790413211295064E+02,
  -3.14512729688483675254357E+04, 
  -3.61444134186911729807069E+04,
   6.64561438202405440627855E+04 };
  bool parity;
  long double pi = 3.141592653589793238462643383279502884;
  long double q[8] = {
  -3.08402300119738975254353E+01,
   3.15350626979604161529144E+02, 
  -1.01515636749021914166146E+03,
  -3.10777167157231109440444E+03, 
   2.25381184209801510330112E+04,
   4.75584627752788110767815E+03, 
  -1.34659959864969306392456E+05,
  -1.15132259675553483497211E+05 };
  long double res;
  long double sqrtpi = 0.9189385332046727417803297;
  long double sum;
  long double value;
  long double xbig = 171.624;
  long double xden;
  long double xinf = 1.79E+308;
  long double xminin = 2.23E-308;
  long double xnum;
  long double y;
  long double y1;
  long double ysq;
  long double z;

  parity = false;
  fact = 1.0;
  n = 0;
  y = x;
//
//  Argument is negative.
//
  if ( y <= 0.0 )
  {
    y = - x;
    y1 = ( long double ) ( int ) ( y );
    res = y - y1;

    if ( res != 0.0 )
    {
      if ( y1 != ( long double ) ( int ) ( y1 * 0.5 ) * 2.0 )
      {
        parity = true;
      }

      fact = - pi / sin ( pi * res );
      y = y + 1.0;
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
//
//  Argument is positive.
//
  if ( y < eps )
  {
//
//  Argument < EPS.
//
    if ( xminin <= y )
    {
      res = 1.0 / y;
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
  else if ( y < 12.0 )
  {
    y1 = y;
//
//  0.0 < argument < 1.0.
//
    if ( y < 1.0 )
    {
      z = y;
      y = y + 1.0;
    }
//
//  1.0 < argument < 12.0.
//  Reduce argument if necessary.
//
    else
    {
      n = ( int ) ( y ) - 1;
      y = y - ( long double ) ( n );
      z = y - 1.0;
    }
//
//  Evaluate approximation for 1.0 < argument < 2.0.
//
    xnum = 0.0;
    xden = 1.0;
    for ( i = 0; i < 8; i++ )
    {
      xnum = ( xnum + p[i] ) * z;
      xden = xden * z + q[i];
    }
    res = xnum / xden + 1.0;
//
//  Adjust result for case  0.0 < argument < 1.0.
//
    if ( y1 < y )
    {
      res = res / y1;
    }
//
//  Adjust result for case 2.0 < argument < 12.0.
//
    else if ( y < y1 )
    {
      for ( i = 1; i <= n; i++ )
      {
        res = res * y;
        y = y + 1.0;
      }
    }
  }
  else
  {
//
//  Evaluate for 12.0 <= argument.
//
    if ( y <= xbig )
    {
      ysq = y * y;
      sum = c[6];
      for ( i = 0; i < 6; i++ )
      {
        sum = sum / ysq + c[i];
      }
      sum = sum / y - y + sqrtpi;
      sum = sum + ( y - 0.5 ) * log ( y );
      res = exp ( sum );
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
//
//  Final adjustments and return.
//
  if ( parity )
  {
    res = - res;
  }

  if ( fact != 1.0 )
  {
    res = fact / res;
  }

  value = res;

  return value;
}
