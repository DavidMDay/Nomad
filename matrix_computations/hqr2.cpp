#include <cassert>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "util.hpp"
void hqr2(std::vector<std::vector<double>>& A, std::vector<double>& tau) {
  auto sz = get_size(A);
  int nn = sz.second;
  tau.resize(nn, 0.0);
  int mm = sz.first;
  for (int j = 0; j < nn; j++) {
    double normx = norm_lower_tri_col(sz, A, j);
    double s = -std::copysign(1.0, A[j][j]);
    double u1 = A[j][j] - s * normx;
    assert(-s * u1 > 0.0);

    std::vector<double> w(mm - j, 0.0);
    w[0] = 1.0;
    for (int row = j + 1; row < mm; row++) {
      w[row - j] = A[row][j] / u1;
    }
    for (int row = j + 1; row < mm; row++) {
      A[row][j] = w[row - j];
    }
    A[j][j] = s * normx;
    tau[j] = -s * u1 / normx;
    if (j + 1 == nn) continue;

    std::vector<double> y(nn - j - 1, 0.0);
    for (int row = j; row < mm; row++) {
      for (int col = j + 1; col < nn; col++) {
        y[col - j - 1] += w[row - j] * A[row][col];
      }
    }

    for (int k = 0; k < nn - j - 1; k++) {
      y[k] *= tau[j];
    }

    for (int row = j; row < mm; row++) {
      for (int col = j + 1; col < nn; col++) {
        A[row][col] -= w[row - j] * y[col - j - 1];
      }
    }
  }
}
