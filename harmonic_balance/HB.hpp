#pragma once

#include <vector>

struct HB {
  double circular_frequency = 1.1e1;
  int order = 1;  // low
  bool linear = true;
  double amp = 1;   // large
  double xi = 0.0;  // damping ratio
};

HB get_duffing();      // Duffing of Liu
HB get_duffing_low();  // Duffing similar to Liu, but lightly damped
HB get_duffing_linear();

// given v(-h/2), v(h/2),  v'(-h/2), v'(h/2)
// return v(0), v'(0), v''(0), v'''(0) such that
// v(x) = v(0) + v'(0) x + v''(0) x^2/2 + v'''(0) x^3/6
// interpolates the data.
std::vector<double> interpolate(double vm, double vp, double dvm, double dvp, double h);

// 2nd equivalent formulation, in terms of Hermite polynomials,
// maybe use this to check the original?
std::vector<double> interpolateHermite(double vm, double vp, double dvm, double dvp, double h);

std::pair<double, double> init_disp(
    double um, double up, double vm, double vp, double dvm, double dvp, double h);
// given  v(kh) < 0 < v((k+1)h),  find t : v(t)=0.
// what is u(t)?  1. [kh, (k+1)h] -> [-h/2,h/2].
// interpolate v(t), a(t) = v'(t)

void duffing_linear(HB xx);