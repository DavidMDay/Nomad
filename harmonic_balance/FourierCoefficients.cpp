#include "FourierCoefficients.hpp"
#include "Tests.hpp"
#include <cassert>
#include <cmath>                            // for abs
#include <complex>                          // for complex, operator*, imag
#include <numeric>                          // for accumulate
#include <ostream>                          // for operator<<, basic_ostream...
#include <iostream>
#include <valarray>                         // for valarray, _SClos, slice
#include <vector>                           // for vector

// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

// The idea here is to handle an array of real Fourier series. 


FourierCoefficients::FourierCoefficients(const std::valarray<double>& inputArray, int neq, int numstep)
    : m_neq(neq), m_numStep(numstep) {
  m_data.resize(numstep * neq, 0.0);
  std::valarray<std::complex<double>> uhat(numstep);

  double normalize0 = 1.0 / static_cast<double>(numstep);
  double normalize = 2.0 * normalize0;
  m_maxBeta0 = 0.0;
  int maxhat = numstep * neq;

  int nmod2 = numstep % 2;

  int jump = Degree() + (1 - nmod2);

  for (int dof = 0; dof < neq; dof++) {
    std::valarray<double> u = inputArray[std::slice(dof, numstep, neq)];
    ThrowRequireMsg(static_cast<int>(u.size()) == numstep, " wrong slice \n");
    Dft(u, uhat);
    m_maxBeta0 = (std::abs(imag(uhat[0])) > m_maxBeta0 ? std::abs(imag(uhat[0])) : m_maxBeta0);
    m_data[dof] = real(uhat[0]) * normalize0;

    for (int order = 1; order <= Degree(); order++) {
      ThrowRequireMsg(dof + (order + jump) * neq < maxhat, " indexing error\n");
      m_data[dof + order * neq] = real(uhat[order]) * normalize;
      m_data[dof + (order + jump) * neq] = -imag(uhat[order]) * normalize;
    }
    if (nmod2 == 0) {
      m_data[dof + jump * neq] = real(uhat[jump]) * normalize0;
    }
  }
}


// simple static inefficient Discrete Fourier transform
void FourierCoefficients::Dft(const std::valarray<double>& x, std::valarray<std::complex<double>>& y) {
  using namespace std::complex_literals;
  int n = static_cast<int>(x.size());
  double r = -(2.0 * M_PI) / static_cast<double>(n);  // W_n =exp(-2pi i/n)
  std::complex<double> zero = 0.0;
  y = zero;
  for (int k = 0; k < n; k++) {
    for (int j = 0; j < n; j++) {
      double angle = r * static_cast<double>(j * k);
      y[k] += x[j] * std::exp(1i * angle);  // Y_k= sum(0<=j<n) X_j    W_n^(jk)
    }
  }
}

// simple static inefficient Discrete Inverse Fourier transform
void FourierCoefficients::Idft(const std::valarray<std::complex<double>>& y,
                               std::valarray<std::complex<double>>& x) {
  using namespace std::complex_literals;
  int n = static_cast<int>(y.size());
  double r = (2.0 * M_PI) / static_cast<double>(n);  // W_n =exp(-2pi i/n)
  std::complex<double> zero = 0.0;
  x = zero;
  for (int j = 0; j < n; j++) {
    for (int k = 0; k < n; k++) {
      double angle = r * static_cast<double>(j * k);
      x[j] += y[k] * std::exp(1i * angle);  // X_j = (1/n) sum(0<=k<n) Y_j    W_n^(-jk)
    }
  }
  std::complex<double> reciprocal = static_cast<std::complex<double>>(1.0 / static_cast<double>(n));
  for (int k = 0; k < n; k++) {
    x[k] *= reciprocal;
  }
}



FourierCoefficients FourierCoefficients::norm1() const {
  FourierCoefficients accumulatedValues;
  accumulatedValues.m_neq = 1;
  accumulatedValues.m_numStep = m_numStep;
  accumulatedValues.m_maxBeta0 = m_maxBeta0;
  accumulatedValues.m_data.resize(m_numStep, 0.0);

  auto norm1 = [](auto sum, auto scalar) { return sum + std::abs(scalar); };
  for (int k = 0; k < m_numStep; k++) {
    accumulatedValues.m_data[k] =
        std::accumulate(m_data.begin() + k * m_neq, m_data.begin() + (k + 1) * m_neq, 0.0, norm1);
  }
  return accumulatedValues;
}

int FourierCoefficients::Degree() const {
  int nmod2 = m_numStep % 2;
  return (nmod2 == 0 ? m_numStep / 2 - 1 : (m_numStep - 1) / 2);
}

double FourierCoefficients::ConstantTerm(int idof) const {
  assert(idof >= 0 && idof < m_neq);
  return m_data[idof];
}

int FourierCoefficients::NumTerms() const {
  int nmod2 = m_numStep % 2;
  if (nmod2 == 0) {
    return Degree() + 1;
  } else {
    return Degree();
  }
}

double FourierCoefficients::CosineTerm(int idof, int istep) const {
  assert(idof >= 0 && idof < m_neq);
  assert(istep >= 1);
  int nmod2 = m_numStep % 2;
  if (istep <= Degree()) {
    return m_data[istep * m_neq + idof];
  } else if (nmod2 == 0 && istep == Degree() + 1) {
    int jump = Degree() + (1 - nmod2);
    ThrowRequireMsg( jump == (m_numStep-nmod2)/2, " confused developer\n");
    return m_data[jump * m_neq + idof];
  } else {
           std::cout << "Invalid arguemnt to FourierCoefficients::CosineTerm" << std::endl;
    return 0.0;
  }
}

std::vector<double> FourierCoefficients::CosineField(int j) const {
  ThrowRequireMsg(j < 0 || j >= m_numStep, "CosineField illegal argument");
  std::vector<double> f( m_data.begin() + j*m_neq, m_data.begin() + (j+1) * m_neq);
  return f;
}

std::vector<double> FourierCoefficients::SineField(int k) const {
  ThrowRequireMsg(k == 0 , "SineField 0 is not defined");
  int nmod2 = m_numStep % 2;
  int jump = (m_numStep-nmod2)/2;
  int j = k + jump;
  return CosineField(j);
}


double FourierCoefficients::SineTerm(int idof, int istep) const {
  assert(idof >= 0 && idof < m_neq);
  assert(istep >= 1);
  int nmod2 = m_numStep % 2;
  int jump = Degree() + (1 - nmod2);
  ThrowRequireMsg( jump == (m_numStep-nmod2)/2, " another confused developer\n");
  if (istep <= Degree()) {
    return m_data[(istep + jump) * m_neq + idof];
  } else if (nmod2 == 0 && istep == Degree() + 1) {
    return 0.0; // ?
  } else {
          std::cout << "Invalid arguemnt to FourierCoefficients::SineTerm" << std::endl;
    return 0.0;
  }
}

void FourierCoefficients::Print(int idof, std::ostream& out) const {
  out << "         cosine             sine \n";
  out << 0 << "       " << ConstantTerm(idof) << "\n";
  for (int k = 1; k <= NumTerms(); ++k) {
    out << k << ": " << CosineTerm(idof, k) << "      " << SineTerm(idof, k) << "\n";
  }
}
