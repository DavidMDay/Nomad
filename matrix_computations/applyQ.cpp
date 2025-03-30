#include <cassert>
#include <iostream>
#include <vector>
#include "util.hpp"

void applyQ(const std::vector<std::vector<double>>& QR,
            const std::vector<double>& tau,
            std::vector<std::vector<double>>& QX) {
  auto sz = get_size(QR);
  auto mm = sz.first;
  auto nn = sz.second;
  auto mm2 = static_cast<int>(QX.size());
  assert(mm2 == mm);
  std::vector<double> y(nn, 0.0);
  for (int j = nn - 1; j >= 0; j--) {
    std::vector<double> w(mm - j, 0.0);
    w[0] = 1.0;
    for (int row = j + 1; row < mm; row++) {
      w[row - j] = QR[row][j];
    }
    std::fill(y.begin(), y.end(), 0.0);
    for (int row = j; row < mm; row++) {
      for (int col = 0; col < nn; col++) {
        y[col] += w[row - j] * QX[row][col];
      }
    }
    for (int k = 0; k < nn; k++) {
      y[k] *= tau[j];
    }
    for (int row = j; row < mm; row++) {
      for (int col = 0; col < nn; col++) {
        QX[row][col] -= w[row - j] * y[col];
      }
    }
  }
}
