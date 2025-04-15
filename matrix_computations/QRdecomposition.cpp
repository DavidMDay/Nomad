#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include "util.hpp"

void QRdecomposition(const std::vector<std::vector<double> >& A,
                     std::vector<std::vector<double> >& Q,
                     std::vector<std::vector<double> >& R) {
  auto sz = get_size(A);
  int m = sz.first;
  Q = std::vector<std::vector<double> >(m, std::vector<double>(m, 0.0));
  for (int i = 0; i < m; i++) {
    Q[i][i] = 1.0;
  }
  std::pair<int, int> szQ = {m, m};
  int inc = 1;
  R = A;
  int n = sz.second;
  for (int j = 0; j < n; j++) {
    double normx = norm_lower_tri_col(sz, R, j);
    double s = -std::copysign(1.0, R[j][j]);
    double u1 = R[j][j] - s * normx;
    assert(-s * u1 > 0.0);
    double scale = 1 / u1;
    std::vector<double> w = extract_column_and_scale(sz, R, j, scale);
    w[j] = 1.0;
    double tau = -s * u1 / normx;
    auto y = gemvt(sz, -tau, R, w, j);
    gerSpecial(sz, w, inc, y, inc, R, j);
    auto x = gemv(sz, -tau, Q, w, j);
    ger(szQ, x, inc, w, inc, Q, j);
  }
}
