#include "util.hpp"
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <iomanip>
#include <string>
#include <utility> // for pair
#include <vector>
// matrix R, m by n, return R(col:end, col) * scale,  over-write w[row]
std::vector<double> extract_column_and_scale(const std::pair<int, int> sz,
                                             const std::vector<std::vector<double>> &R,
                                             const int col,
                                             double scale) {
  assert(col < sz.second);
  auto m = sz.first;
  std::vector<double> w(m, 0.0);
  for (int row = col; row < sz.first; row++) {
    w[row] = R[row][col] * scale;
  }
  return w;
}

std::pair<int, int> get_size(const std::vector<std::vector<double>> &A) {
  int m = static_cast<int>(A.size());
  int n = static_cast<int>(A[0].size());
  return {m, n};
}

// <matrix(row,start:end), vector(start:end)>
double ip_mat_row_vec(const std::vector<double> &w,
                      const std::pair<int, int> sz,
                      const std::vector<std::vector<double>> &Q,
                      const int row,
                      const int first_col) {
  assert(row < sz.first);
  double innerProduct = 0.0;
  for (int col = first_col; col < sz.first; col++) {
    innerProduct += Q[row][col] * w[col];
  }
  return innerProduct;
}

// m-by-n R, return w(first_row:end)  R( first_row:end, col)
double ip_vec_mat_column(const std::vector<double> &w,
                         const std::pair<int, int> sz,
                         const std::vector<std::vector<double>> &R,
                         const int first_row,
                         const int column) {
  double innerProduct = 0.0;
  // std::vector wR(sz.second -column,0.0);

  for (int row = first_row; row < sz.first; row++) {
    innerProduct += w[row] * R[row][column];
  }
  return innerProduct;
}

// norm_lower_tri_col, m-by-n R, return norm( R(col:end, col)
double norm_lower_tri_col(const std::pair<int, int> sz,
                          const std::vector<std::vector<double>> &R,
                          const int col) {
  assert(col < sz.second);
  double normx = 0.0;
  for (int row = col; row < sz.first; row++) {
    normx += R[row][col] * R[row][col];
  }
  normx = std::sqrt(normx);
  return normx;
}

double norm1(const std::pair<int, int> sz, const std::vector<std::vector<double>> &B) {
  double x = 0.0;
  for (int row = 0; row < sz.first; row++) {
    double xrow = 0.0;
    for (int col = 0; col < sz.second; col++) {
      xrow = std::max(xrow, std::abs(B[row][col]));
    }
    x = std::max(x, xrow);
  }
  return x;
}

double normInf(const std::pair<int, int> sz, const std::vector<std::vector<double>> &B) {
  double x = 0.0;
  for (int row = 0; row < sz.first; row++) {
    double xrow = 0.0;
    for (int col = 0; col < sz.second; col++) {
      xrow += std::abs(B[row][col]);
    }
    x = std::max(x, xrow);
  }
  return x;
}
// tril is a MATLAB function that returns the lower triangular part
// of a matrix.  This implementation of tril sets the upper triangular
// part of a matrix to zero, leaving the strictly lower triangular part.
void tril(std::pair<int, int> sz, std::vector<std::vector<double>> &X) {
  for (int row = 0; row < sz.second; row++) {
    for (int col = row; col < sz.second; col++) {
      X[row][col] = 0.0;
    }
  }
}
// m-by-n X,  R = triu(X(1:n,:));
std::vector<std::vector<double>> triu(std::pair<int, int> sz, std::vector<std::vector<double>> &X) {
  auto n = sz.second;
  assert(sz.first >= n);
  auto R = std::vector<std::vector<double>>(n, std::vector<double>(n, 0.0));
  for (int row = 0; row < n; row++) {
    for (int col = row; col < n; col++) {
      R[row][col] = X[row][col];
    }
  }
  return R;
}

std::vector<double> gemv(std::pair<int, int> sz,
                         const double alpha,
                         const std::vector<std::vector<double>> &B,
                         const std::vector<double> &w,
                         const int start) {
  std::vector<double> y(sz.first, 0.0);
  for (int row = 0; row < sz.first; row++) {
    for (int col = start; col < sz.first; col++) {
      y[row] += B[row][col] * w[col];
    }
    y[row] *= alpha;
  }
  return y;
}

// m-by-1 w, m-by-n B, return w(start:end)  B( start:end,start:end)
std::vector<double> gemvt(std::pair<int, int> sz,
                          const double alpha,
                          const std::vector<std::vector<double>> &B,
                          const std::vector<double> &w,
                          const int start) {
  std::vector<double> y(sz.second - start, 0.0);
  for (int row = start; row < sz.first; row++) {
    for (int column = start; column < sz.second; column++) {
      y[column - start] += w[row] * B[row][column];
    }
  }
  for (int column = start; column < sz.second; column++) {
    y[column - start] *= alpha;
  }
  return y;
}

void trtrs(int /*unused*/, const std::vector<std::vector<double>> &R, std::vector<double> &rhs) {
  auto sz = get_size(R);
  for (int row = sz.second - 1; row >= 0; row--) {
    for (int col = row + 1; col < sz.second; col++) {
      rhs[row] -= R[row][col] * rhs[col];
    }
    rhs[row] /= R[row][row];
  }
}

void axpy(double alpha,
          std::pair<int, int> sz,
          const std::vector<std::vector<double>> &X,
          std::vector<std::vector<double>> &Y) {
  for (int row = 0; row < sz.first; row++) {
    for (int col = 0; col < sz.second; col++) {
      Y[row][col] = Y[row][col] + X[row][col] * alpha;
    }
  }
}

// assert(szLeft.second == szRight.first);
// P = std::vector<std::vector<double>>(szLeft.first,
// std::vector<double>(szRight.second,0.0));
void gemm_nn(double alpha,
             std::pair<int, int> szLeft,
             const std::vector<std::vector<double>> &Left,
             std::pair<int, int> szRight,
             const std::vector<std::vector<double>> &Right,
             double beta,
             std::vector<std::vector<double>> &P) {
  for (int row = 0; row < szLeft.first; row++) {
    for (int col = 0; col < szRight.second; col++) {
      double s = 0;
      for (int k = 0; k < szRight.first; k++) {
        s += Left[row][k] * Right[k][col];
      }
      P[row][col] = s * alpha + P[row][col] * beta;
    }
  }
}

// assert(szLeft.first == szRight.first);
// A = std::vector<std::vector<double>>(szLeft.second,
// std::vector<double>(szRight.second, 0.0));
void gemm_tn(double alpha,
             std::pair<int, int> szLeft,
             const std::vector<std::vector<double>> &Left,
             std::pair<int, int> szRight,
             const std::vector<std::vector<double>> &Right,
             double beta,
             std::vector<std::vector<double>> &P) {
  for (int row = 0; row < szLeft.second; row++) {
    for (int col = 0; col < szRight.second; col++) {
      double s = 0;
      for (int k = 0; k < szRight.first; k++) {
        s += Left[k][row] * Right[k][col];
      }
      P[row][col] = s * alpha + P[row][col] * beta;
    }
  }
}

void gerSpecial(std::pair<int, int> sz,  // m-by-n
                std::vector<double> &x,
                int /*unused*/,  // size_x >= 1 + (m - 1)*incx;
                std::vector<double> &y,
                int /*unused*/,  // size_y >= 1 + (n - 1)*incy;
                std::vector<std::vector<double>> &B,
                int start) {
  for (int row = start; row < sz.first; row++) {
    for (int column = start; column < sz.second; column++) {
      B[row][column] += x[row] * y[column - start];  // weird
    }
  }
}
// BLAS2 ger + "start".
void ger(std::pair<int, int> sz,  // m-by-n
         std::vector<double> &x,
         int /*unused*/,  // size_x >= 1 + (m - 1)*incx;
         std::vector<double> &y,
         int /*unused*/,  // size_y >= 1 + (n - 1)*incy;
         std::vector<std::vector<double>> &B,
         int start) {
  for (int row = 0; row < sz.first; row++) {
    for (int column = start; column < sz.second; column++) {
      B[row][column] += x[row] * y[column];
    }
  }
}

void print_matrix(const std::string &name, const std::vector<std::vector<double>> &M) {
  auto sz = get_size(M);
  std::cout << name << " is " << sz.first << " by " << sz.second << std::endl;

  for (const auto &row : M) {
    for (const auto &element : row) {
      std::cout << element << " ";
      // std::cout << std::scientific << std::setprecision(10) << element << "   ";
    }
    std::cout << std::endl;
  }
}

void print_square_matrix_col(std::vector<double> m) {
  auto o = static_cast<size_t>(round(std::sqrt(static_cast<double>(m.size())) ));
  assert(o * o == m.size());
  for (size_t row = 0; row < o; row++) {
    for (size_t col = 0; col < o; col++) {
      double x = (std::abs(m[row + o * col]) < 1.e-12 ? 0.0 : m[row + o * col]);
      std::cout << x;
      if (col < o - 1) {
        std::cout << "   ";
      }
    }
    std::cout << "\n";
  }
}

void scale_vector(std::vector<double> &v, double scale) {
  for (double & i : v) {
    i *= scale;
  }
}

std::vector<double> square_matrix_transpose(const std::vector<double> &B) {
  long n = static_cast<long>(round(std::sqrt(static_cast<double>(B.size())) ));
  std::vector<double> T(n * n, 0.0);
  for (int col = 0; col < n; col++) {
    for (int row = 0; row < n; row++) {
      T[col + n * row] = B[row + n * col];
    }
  }
  return T;
}

void axpy(const std::vector<double> &x, double a, std::vector<double> &y) {
  for (size_t i = 0; i < y.size(); i++) {
    y[i] += a * x[i];
  }
}

double normInf(const std::vector<double>& x) {
  double r = 0.0;
  for (double i : x) {
    r = std::max(r, std::abs(i));
  }
  return r;
}

// return B(:,start:end) * w * alpha
std::vector<double> gemv(std::pair<int, int> sz,
                         double alpha,
                         const std::vector<double> &B,
                         const std::vector<double> &w,
                         int start) {
  auto B2 = matrixOfVector(false, sz, B);
  return gemv(sz, alpha, B2, w, start);
}