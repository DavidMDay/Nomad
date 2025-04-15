#pragma once

#include <complex>  // for complex

void ThrowRequireMsg(bool worked, std::string description);

void EXPECT_NEAR(double a, double b, double tol);

double modulus(std::complex<double> z);
std::complex<double> conjugate(std::complex<double> z);
struct material {
  double E = 0;
  double density = 0;
  double poissons_ratio = 0;
};

int TestOddDFT();
int TestEvenDFT();
int TestLiuOde();
int TestLightlyDampedOde();
int TestSpaceTimeTransOde();
int TestIncreasedOmega();
int TestIncreasedOmegaLiu();
int TestLinearized();
int TestOperators();
