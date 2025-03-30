
#include "matrix2x2.hpp"
#include <vector>
#include <stddef.h>  // for size_t
#include <cassert>
#include <cmath>  // for std::abs
#include "util.hpp" // for scale_vector


double determinant2x2( const std::vector<double>& matrix) {
    return matrix[0] * matrix[3] - matrix[1] * matrix[2];
}

std::vector<double> mat_vec_2x2(const std::vector<double> &x, const std::vector<double> &matrix)
{
  std::vector<double> y = {
      matrix[0] * x[0] + matrix[2] * x[1], matrix[1] * x[0] + matrix[3] * x[1]};
  return y;
}

std::vector<double>
solve_2by2(const std::vector<double> &rhs, const std::vector<double> &matrix)
{
  std::vector<double> inv = {matrix[3], -matrix[1], -matrix[2], matrix[0]};
  double det = determinant2x2( matrix);
  assert(std::abs(det) > 0.0);
  scale_vector(inv, 1.0/det);
  std::vector<double> solution = mat_vec_2x2(rhs, inv);
  return solution;
}
std::vector<double>
residual_2by2(const std::vector<double> &rhs, const std::vector<double> &matrix,
              const std::vector<double> &sol)
{
  std::vector<double> y = mat_vec_2x2(sol, matrix);
  std::vector<double> resid = {rhs[0] - y[0], rhs[1] - y[1]};
  return resid;
}