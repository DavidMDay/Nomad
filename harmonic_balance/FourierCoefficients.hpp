#pragma once
#include <algorithm>  // for fill, max, copy
#include <complex>    // for complex
#include <iosfwd>     // for ostream
#include <valarray>   // for valarray
#include <vector>     // for vector

// Motivation, intent: a systems of neq odes, the solution
// is known at numstep equally spaced points:  0,1,...., numstep-1
// points.  We'll need to compute a Fourier series of each component.
// The period of the Fourier series corresponds to time step numstep.
// The actual time step size does not matter.

class FourierCoefficients {
 public:
  FourierCoefficients()=default;
  FourierCoefficients(const std::valarray<double>& inputArray, int neq, int numstep);
  double GetMaxBeta0() const { return m_maxBeta0; }
  FourierCoefficients norm1() const;
  void Print(int idof, std::ostream& out) const;
  int NumTerms() const;
  double ConstantTerm(int idof) const;
  double CosineTerm(int idof, int istep) const;
  double SineTerm(int idof, int istep) const;
  std::vector<double> CosineField(int) const;
  std::vector<double> SineField(int) const;
  static void Dft(const std::valarray<double>& x, std::valarray<std::complex<double>>& y);
  static void Idft(const std::valarray<std::complex<double>>& y, std::valarray<std::complex<double>>& x);

 private:
  int Degree() const;
  int m_neq = 0;
  int m_numStep = 0;
  std::vector<double> m_data;
  double m_maxBeta0 = 0.0;
};
