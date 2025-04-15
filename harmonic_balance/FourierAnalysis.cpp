#include <cassert>
#include <algorithm>  // for fill, max, copy
#include <valarray>   // for valarray
#include <vector>     // for vector
#include <cstdlib>
#include <complex>
#include "FourierCoefficients.hpp"
#include "Tests.hpp"


int TestOddDFT() {
  int n = 7;
  int nmod2 = n % 2;
  double half = 0.5 * static_cast<double>(n);
  double tol = 1.e-14;
  std::valarray<double> constant(1.0, n);
  std::valarray<std::complex<double>> c0_hat(0.0, n);
  FourierCoefficients::Dft(constant, c0_hat);
  EXPECT_NEAR(1.0, real(c0_hat[0])/static_cast<double>(n), tol);

  for (int i=1;i<n;i++) {
    EXPECT_NEAR(0, modulus(c0_hat[i]), tol);
  }

  for (int i=1; i<n; i++ ) {
    EXPECT_NEAR(0.0, modulus(c0_hat[i]), tol);
  }
  int degree = (n - nmod2) / 2;

  std::valarray<double> c(0.0, n);
  std::valarray<double> s(0.0, n);
  for (int L = 1; L <= degree; L++)
  {
    for (int j = 0; j < n; j++)
    {
      int Lj = (L*j)%n;// skirt floating point arithmetic issue
      double a = M_PI * static_cast<double>(2 * Lj) / static_cast<double>(n);
      c[j] = cos(a);
      s[j] = sin(a);
    }
    std::valarray<std::complex<double>> c_hat(0.0, n);
    std::valarray<std::complex<double>> s_hat(0.0, n);
    FourierCoefficients::Dft(c, c_hat);
    EXPECT_NEAR(half, real(c_hat[L]), tol);
    EXPECT_NEAR(half, real(c_hat[n - L]), tol);
    FourierCoefficients::Dft(s, s_hat);
    EXPECT_NEAR(-half, imag(s_hat[L]), tol);
    EXPECT_NEAR(half, imag(s_hat[n - L]), tol);
    std::valarray<std::complex<double>> c2(0.0, n);
    std::valarray<std::complex<double>> s2(0.0, n);
    FourierCoefficients::Idft(c_hat, c2);
    FourierCoefficients::Idft(s_hat, s2);
    for (int j = 1; j < n; j++)
    {
      EXPECT_NEAR(c[j], real(c2[j]), tol);
      EXPECT_NEAR(s[j], real(s2[j]), tol);
    }
  } // L
  return EXIT_SUCCESS;
}

int TestEvenDFT()
{
  int n = 6;
  int nmod2 = n % 2;
  double half = .5 * static_cast<double>(n);
  double tau = 1.e-14;
  std::valarray<double> c(1.0, n);
  std::valarray<double> s(0.0, n);
  std::valarray<std::complex<double>> c0_hat(0.0, n);
  FourierCoefficients::Dft(c, c0_hat);
  auto whole = static_cast<double>(n);
  EXPECT_NEAR(whole, real(c0_hat[0]), tau);
  int degree = (nmod2 == 0 ? n / 2 - 1 : (n - 1) / 2);

  for (int L = 1; L <= degree; L++) {
    for (int j = 0; j < n; j++) {
      int Lj = (L*j)%n;
      double a = M_PI * static_cast<double>(2 * Lj) / static_cast<double>(n);
      c[j] = cos(a);
      s[j] = sin(a);
    }
    std::valarray<std::complex<double>> c_hat(0.0, n);
    std::valarray<std::complex<double>> s_hat(0.0, n);
    FourierCoefficients::Dft(c, c_hat);
    EXPECT_NEAR(half, real(c_hat[L]), tau);
    EXPECT_NEAR(half, real(c_hat[n - L]), tau);
    FourierCoefficients::Dft(s, s_hat);
    EXPECT_NEAR(-half, imag(s_hat[L]), tau);
    EXPECT_NEAR(half, imag(s_hat[n - L]), tau);
    std::valarray<std::complex<double>> c2(0.0, n);
    std::valarray<std::complex<double>> s2(0.0, n);
    FourierCoefficients::Idft(c_hat, c2);
    FourierCoefficients::Idft(s_hat, s2);
    for (int j = 1; j < n; j++) {
      EXPECT_NEAR(c[j], real(c2[j]), tau);
      EXPECT_NEAR(s[j], real(s2[j]), tau);
    }
  } // L
  int L = degree + 1;
  for (int j = 0; j < n; j++)
  {
    double a = M_PI * static_cast<double>(2 * L * j) / static_cast<double>(n);
    c[j] = cos(a);
  }
  c0_hat = 0.0;
  FourierCoefficients::Dft(c, c0_hat);
  EXPECT_NEAR(whole, real(c0_hat[L]), tau);
  return EXIT_SUCCESS;
}

int TestHelp() {
  int neq(1);
  int nstep(1);
  double zero(0.0);
  std::valarray<double> x(zero, static_cast<size_t>(neq * nstep));
  FourierCoefficients c = FourierCoefficients(x, neq, nstep);
  assert(0 == c.NumTerms());
  return EXIT_SUCCESS;
}
