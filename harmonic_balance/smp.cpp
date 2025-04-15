#include "smp.hpp"

#include <cmath>
#include <iostream>
#include <utility>
#include <vector>

#include "util.hpp"
// column major order,  (i,j) i + nr*j

std::vector<double> square_matrix_product_col(const std::vector<double> &Lvec,
                                              const std::vector<double> &Rvec) {
  int n = std::sqrt(static_cast<int>(Lvec.size()));
  std::pair<int, int> szLeft = {n, n};
  std::vector<std::vector<double>> Left = matrixOfVector(true, szLeft, Lvec);
  std::pair<int, int> szRight = szLeft;
  std::vector<std::vector<double>> Right = matrixOfVector(true, szRight, Rvec);
  std::vector<std::vector<double>> matrix(Right);
  gemm_nn(1.0, szLeft, Left, szRight, Right, 0.0, matrix);
  std::vector<double> P = vectorOfMatrix(true, matrix);
  return P;
}
