#include "Operators.hpp"
#include <vector>
#include <iostream>
#include <cmath>  // for std::sqrt
#include <cassert>
#include "util.hpp"
#include "smp.hpp"
namespace {

  double sep(const std::vector<double> &A, std::vector<double> B);
  double diffIdenity(const std::vector<double> &A);

  double sep(const std::vector<double> &A, std::vector<double> B)
  {
    axpy(A, -1.0, B);
    return normInf(B);
  }
  double diffIdenity(const std::vector<double> &A)
  {
    std::vector<double> Id(A);
    auto n = std::sqrt(Id.size());
    for (int i = 0; i < n; i++)
    {
      Id[i + n * i] -= 1.0;
    }
    return normInf(Id);
  }

  std::vector<double> E1()
  {
    constexpr double a = 1.0 / 3.0;
    const double b = std::sqrt(1.0 / 3.0);
    std::vector<double> e = {a, 2.0 * a, 0, a, -a, b, a, -a, -b};
    return e;
  }

  double e1sep()
  {
    std::vector<double> E = getE(1);
    std::vector<double> e1 = E1();
    return sep(E, e1);
  }

  std::vector<double> invE1()
  {
    constexpr double c = 0.5;
    const double s = 0.5 * std::sqrt(3.0);
    std::vector<double> ei = {1.0, 1.0, 1.0, 1.0, -c, -c, 0, s, -s};
    return ei;
  }

  double e1invsep()
  {
    std::vector<double> Ei = getInvE(1);
    std::vector<double> ei1 = invE1();
    return sep(Ei, ei1);
  }

  std::vector<double> D1()
  {
    const double s = 1.0/std::sqrt(3.0);
    std::vector<double> d = {0, s, -s, -s, 0, s, s, -s, 0};
    return d;
  }

  double d1sep()
  {
    std::vector<double> D = getD(1);
    std::vector<double> d1 = D1();
    return sep(D, d1);
  }

double RealDFSTest(const std::vector<double> &E, const std::vector<double> &invE)
{
  int n = std::sqrt(static_cast<int>(E.size()));
  std::vector<double> X(E);
  for (int col = 0; col < n; col++)
  {
    for (int row = 0; row < n; row++)
    {
      if (row == 0)
      {
        X[row + n * col] = E[row + n * col] * static_cast<double>(n);
      }
      else
      {
        X[row + n * col] = E[row + n * col] * static_cast<double>(n) * 0.5;
      }
    }
  }
  auto XT = square_matrix_transpose(X);
  return sep(invE,XT);
}

int TestLiuOperators(int order)
{
  constexpr double tol = 1.e-12;
  std::vector<double> A = getA(order);
  std::vector<double> A2 = getAsquared(order);
  std::vector<double> PA = square_matrix_product_col(A,A);
  assert( sep(A2,PA) < tol) ;
  std::vector<double> E = getE(order);
  std::vector<double> invE = getInvE(order);
  assert( RealDFSTest(E, invE ) < tol );
  std::vector<double> Id = square_matrix_product_col(E,invE);
  assert(  diffIdenity(Id)  < tol );
  std::vector<double> AE = square_matrix_product_col(A,E);
  return EXIT_SUCCESS;
}
} // namespace

int TestOperators()
{
  constexpr double tol = 1.e-12;
  assert( e1sep() < tol);
  assert( e1invsep() < tol);
  assert( d1sep() < tol);
  for (int order = 0; order < 5; order++)
  {
    if (TestLiuOperators(order) == EXIT_FAILURE)
    {
      return EXIT_FAILURE;
    }
  }
  return EXIT_SUCCESS;
}
