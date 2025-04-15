
#include <cmath>
#include <cassert>
#include <algorithm> // for fill, max, copy
#include <iostream>
#include <complex>    // for complex

void ThrowRequireMsg(bool worked, std::string description) {
  if (!worked) {
    std::cout << description << std::endl;
  }
}

void EXPECT_NEAR(double a, double b, double tol) { assert(std::abs(b - a) < tol); }

std::complex<double> conjugate(std::complex<double> z) { return {real(z), -imag(z)}; }

double modulus(std::complex<double> z) { return sqrt(real(z * conjugate(z))); }
