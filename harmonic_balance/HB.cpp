#include "HB.hpp"
#include <vector>
#include <iostream>
#include <cassert>
#include <cmath>  // for std::atan
#include <utility>
#include <complex>

// This is what I thought that Liu used, but the damping ratio is lower.

HB get_duffing_low() {
    HB duff;
    duff.amp = 0.025;// f(t) = amp * sin( omega t )
    duff.circular_frequency = 1.0;
    duff.linear = false;
    duff.xi = 0.01;// u'' + 2 xi u' + u + u^3 = f(t)
    return duff;
  }

// This is what Liu uses.  Everyone seems to define things
// differently.
// u'' + 2 xi u' + u + u^3 = f(t),  f(t) = a sin(wt)
HB get_duffing() {
    HB duff;
    duff.amp = 0.025;// 1/4, 5/4
    duff.circular_frequency = 1.0;
    duff.linear = false;
    duff.xi = 0.1;// u'' + 2 xi u' + u + u^3 = f(t)
    return duff;
  }

  HB get_duffing_linear()
  {
    HB duff;
    duff.amp = 0.025; // 1/4, 5/4
    duff.circular_frequency = 0.25;
    duff.linear = true;
    duff.xi = 0.1; // u'' + 2 xi u' + u + u^3 = f(t)
    return duff;
  }

std::vector<double>
interpolate( double vm, double vp, double dvm, double dvp, double h) {
 double d2 = (dvp - dvm)/h; //v''(0) = (v'(h/2) - v'(-h/2))/h
 // v(0) = ( v(h/2) + v(-h/2) - v''(0) h^2/4 )/2
 double v = 0.5*( vp + vm - 0.25 * d2 * h * h);
  /*
 //v(x) - v(-x) = v'(0) 2x + v'''(0) x^3/3
 // v'(x) + v'(-x) = 2 v'(0) + v'''(0) x^2
 // d = v'(0),   c = v'''(0)
 */
 double diff = vp - vm;  // a = v(x) - v(-x)
 double sum =  dvp + dvm; // b = v'(x) + v'(-x)
 /*
 2x    x^3/3  d  = a
 2      x^2   c    b,     y = 1/x,
 [2x   x^3/3]  [3y/4        -1/4]
 [2      x^2]  [-3y^3/2   3y^2/2] = I

 v'(0) = 3ay/4 -b/4
 v'''(0) = -(3a/2)y^3  + (3b/2) y^2
 x = h/2,  y = 2/h,  3y/4 = 3/(2h)
 (3/2)y^2 = 3 * 4/(2 h^2) = 6/h^2
 (3/2)y^3 = 3 * 8/(2 h^3) = 12/h^3
 */
 double d = 1.5*diff/h - 0.25*sum;  // v'(0) = 3a/(2h) -b/4
 double d3 = (sum - 2.0*diff/h)*6.0/(h*h);
 // v'''(0) = -12 a/h^3  + 6b/h^2= (b-2a/h)6/(h^2)
 std::vector<double> coef = {v,d,d2,d3};
 return coef;
}


std::vector<double>
interpolateHermite( double vm, double vp, double dvm, double dvp, double h) {

double x = 0.5*h;
double y = 1/x;
double y2 = y*y;
double y3 = y*y2;
std::vector<double> vec = {0.25 * vm* y3, 0.25*vp*y3, 4.0*dvm*y2, 4.0*dvp*y2};

double x2 = 1.0/y2;
double x3 = 1.0/y3;
std::vector<double> mat_col = {
1.0,      0,    -3*x2,2*x3,
-1.0,0,3*x2,2*x3,
1.0,-x,-x2,x3,
1.0,x,-x2,x3};
(void) mat_col;  // really return mat_col * vec.
return vec;
}

std::pair<double,double>
init_disp( double um, double up, double vm, double vp, double dvm, double dvp, double h) {
  std::vector<double> coef = interpolate(vm,vp,dvm,dvp,h);
  double v = coef[0];
  double d = coef[1];
  double d2 = coef[2];
  double d3 = coef[3];
  assert( vm*vp <= 0.0);
  double vmin = v;
  double amin = d;
  double t1 = -v/d;// v + d x = 0,   x = -v/d
  double u1 = 0;
  bool verbose = false;
  bool normal =  ( vm > 0.0 && vp < 0.0 );
  if (!normal) {
	    std::cout << " left bracket  v("<< -h/2 <<")= " << vm << "\n";
	    std::cout << " rightt bracket  v("<< h/2 <<")= " << vp << "\n";
  }
  assert( normal );  //  v(-h/2) = vm < 0,   v(h/2) = vp > 0
  std::vector<double> cu = interpolate(um,up,vm,vp,h);
  if ( std::abs(t1) < 0.5*h )  {
    double a1 = d + t1*( d2 + 0.5*d3*t1); // v'(x) =  d + d2 x + d3 x2/2
    double v1 = v + t1*( d + 0.5*t1*( d2 + d3*t1/3.0));
    if ( std::abs(v1) >= std::abs(vmin)) {
      std::cout << "v(0)=" << v << " smaller than v("<<t1<<")="<< v1 << "\n";
    }
    assert ( std::abs(v1) < std::abs(vmin));
    double argmin = t1;
    vmin = v1;
    amin = a1;
    // The quadratic coefficient is independent of the
    // values at the endpoints.  In other words, the quadratic part
    // of the cubit interpolant is not helpful here.
    t1 = argmin;
    for (int iter = 0; iter < 2; iter++)
    {
      double dt = -vmin / amin;
      t1 += dt;
      if ( t1 < -0.5*h || t1 > 0.5*h ) {
        std::cout <<  "iteration " << iter << "  interval " << -0.5*h << "," << 0.5*h << " init Newton step " << t1 << "\n";
        verbose = true;
      }
      u1 = cu[0] + t1 * (cu[1] + 0.5 * t1 * (cu[2] + cu[3] * t1 / 3.0));
      vmin = v + t1 * (d + 0.5 * t1 * (d2 + d3 * t1 / 3.0));
      amin = d + t1 * (d2 + 0.5 * d3 * t1); // v'(x) =  d + d2 x + d3 x2/2
      if (verbose) {
        std::cout << "iteration" << iter << "  t " << t1 << "  u " << u1 << "  v " << vmin << "  a " << amin << "\n";
      }
    }
  } else {
	  // only needed once, when there was a negative local maximum of u(t),  u(t1
    std::cout <<  " init_disp: first Newton step " << t1 <<
    " not in " << -0.5*h << "," << 0.5*h  << "\n";
    verbose = true;
    double mid_val = v;  // left_val = vm, right_val  = vp,
    double left = -0.5*h, right  = 0.5*h, mid = 0.0;
    if (verbose) {
      std::cout << "bisection\n";
    }
    for (int iter = 0; iter < 10; iter++) {
      if ( mid_val < 0.0 ) {
        right = mid;
      } else {
        left = mid;
      }
      mid = 0.5*(left+right);
      mid_val = v + mid * (d + 0.5 * mid * (d2 + d3 * mid / 3.0));
      if (verbose) {
        std::cout << " iter " << iter << "  t  " << mid << "   v  " << mid_val << "\n";
      }
    }
    t1 = mid;
    u1 = cu[0] + t1 * (cu[1] + 0.5 * t1 * (cu[2] + cu[3] * t1 / 3.0));
  }
  return {t1,u1};
}
// The idea here is to have an atan function
// with range  -pi , 0.  Yes, this could be replaced
// by atan2.
double my_atan(double reH, double imH)
{
  assert(imH >= 0.0);
  constexpr double tol = 1.e-14;
  if (reH > tol)
  {
    return atan(-imH / reH);
  }
  if (reH < -tol)
  {
    return atan(-imH / reH) - M_PI;
  }
  return -0.5 * M_PI;
}

// linear case: drop u^3 term.  duffing_linear returns
// helpful information about the periodic solution
// 0 < xi << 1,    w ~ 1  (relative angular frequency, near resonance), a > 0.
// H = -w^2 + 2 sqrt(-1) xi w + 1, b=a/| H(w) |,
// quantity of interest: sup_t u(t) ~ a/xi, and 'quality' ~ 1/qoi .
// u(t) = b sin(w t + p),  
//  u'  =  b w cos(wt+p)
//  u''  = -w^2 u,   H = X + i Y =  H co + i H so
//  H b sin(wt+p) co + H b cos(wt+p) so = a sin(wt),   b = a/H
//  co = cos(-p)= (1-w^2)/H   so = sin(-p) = 2 xi w/H  > 0
void duffing_linear(HB xx) {
assert( xx.amp > 0.0 );
assert( 0 <  xx.xi  && xx.xi < 0.25 );
auto w = xx.circular_frequency;
auto reH = 1.0-w*w;
auto imH = 2.0*xx.xi*w;
constexpr double tol = 1.e-14;
auto phase = my_atan(reH, imH);
std::complex<double> H = {reH, imH};
auto ah = std::abs(H);
assert( std::abs( cos(-phase) - reH/ah ) < tol );
assert( std::abs( sin(-phase) - imH/ah ) < tol );

// modulus:  m(w) = |H(w)|,
// auto w_worst = std::sqrt(1 - 2.0* xx.xi * xx.xi );
//   m(w) => m(v) =  x sqrt(1- x^2)

auto to = (0.5 * M_PI - phase )/ w;
// b = sup u(t) = a/m(w),  m(w) = |H(w)|
auto b = xx.amp/ah;
double uo = b*sin(phase);
double vo = w*b*cos(phase);
assert( std::abs( sin(phase)*reH + cos(phase)*imH ) < tol );
std::cout << "(u,v,w)(0) "  << uo << ", "  << vo << " , " << -w*w*uo << "\n";
std::cout << " u( " << to << ") =" << b <<  " v( " << to << ") =0\n";
}