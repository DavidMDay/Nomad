#include "getSquareUpperTri.hpp"
#include <vector>
#include <utility>
#include <cmath>
#include <cassert>
#include <iostream>
//#include "triuSol.hpp"  // for apply_triu, apply_inv_triu

std::vector<std::vector<double>> getSquareUpperTri(std::pair<int, int> sz,
                                                   double diag,
                                                   double upper_diag) {
  assert(sz.first == sz.second);
  int n = sz.first;
  auto R = std::vector<std::vector<double>>(n, std::vector<double>(n, 0.0));
  for (int row = 0; row < n; row++) {
    R[row][row] = diag;
    for (int col = row + 1; col < n; col++) {
      R[row][col] = upper_diag;
    }
  }
  return R;
}
// NOLINTBEGIN
// x  u   u   u      y   -v1   v2 -v3
// 0  x   u   u      0    y   -v1  v2
// 0  0   x   u      0    0    y  -v1
// 0  0   0   x      0    0    0    y
// v0 = 1/x
// v1 = v0 * u * v0;
// v2 = v1 *(u*v0-1)
// v3 = v2 *(u*v0-1)
// NOLINTEND
std::vector<std::vector<double>> getSquareUpperTriInv(std::pair<int, int> sz,
                                                      double diag,
                                                      double upper_diag) {
  assert(sz.first == sz.second);
  int n = sz.first;
  if (n <= 0) {
    std::vector<std::vector<double>> R;
    return R;
  }
  assert(diag != 0.0);
  std::vector<double> v(n, 1.0 / diag);
  if (n == 1) {
    auto R = std::vector<std::vector<double>>(n, std::vector<double>(n, 0.0));
    R[0][0] = 1.0 / diag;
    return R;
  }
  v[1] = v[0] * (upper_diag * v[0]);
  for (int i = 2; i < n; i++) {
    v[i] = v[i - 1] * (upper_diag * v[0] - 1.0);
  }
  auto R = std::vector<std::vector<double>>(n, std::vector<double>(n, 0.0));
  for (int row = 0; row < n; row++) {
    for (int col = row; col < n; col++) {
      R[row][col] = v[col - row];
    }
  }
  return R;
}
