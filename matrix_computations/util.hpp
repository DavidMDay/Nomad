#pragma once
#include <string> // retab me
#include <vector>
#include <utility> // for pair

// dense matrix utils and other things too
// One representation of a dense matrix is vector<vector<T>>
// The inner vector is a row and the outer vector keeps track
// of the rows.
// Another is as a vector.  The vector could come from a row
// major (op_row) or column major (op_col) ordering.
// The integer pair = number of rows, columns


std::vector<double>
extract_column_and_scale(std::pair<int, int> sz,
                                             const std::vector<std::vector<double>> &R,
                                             int col,
                                             double scale);

std::pair<int, int>
get_size(const std::vector<std::vector<double>> &A);

// row-major ordering
template <typename T>
std::vector<T> vectorOfMatrix(bool transpose, const std::vector<std::vector<T>>& matrix) {
    auto sz = get_size(matrix);
    std::vector<T> lengthy;
    lengthy.reserve(sz.first * sz.second);
    if (transpose) {
      for (int j = 0; j < sz.second; ++j) {
        for (int i = 0; i < sz.first; ++i) {
          lengthy.push_back(matrix[i][j]);
        }
      }
    } else {
      for (int i = 0; i < sz.first; ++i) {
        for (int j = 0; j < sz.second; ++j) {
          lengthy.push_back(matrix[i][j]);
        }
      }
    }
    return lengthy;
}

// row-major ordering
template <typename T>
std::vector<std::vector<T>> matrixOfVector(bool transpose, std::pair<int,int> sz, const std::vector<T>& v) {
    auto X = std::vector<std::vector<double>>(sz.first, std::vector<double>(sz.second, 0.0));
    int k = 0;
    if (transpose) {
      for (int j = 0; j < sz.second; ++j) {
        for (int i = 0; i < sz.first; ++i) {
          X[i][j] = v[k++];
        }
      }
    } else {
      for (int i = 0; i < sz.first; ++i) {
        for (int j = 0; j < sz.second; ++j) {
          X[i][j] = v[k++];
        }
      }
    }
    return X;
}


double ip_mat_row_vec(const std::vector<double> &w,
                      std::pair<int, int> sz,
                      const std::vector<std::vector<double>> &Q,
                      int row,
                      int first_col);

double ip_vec_mat_column(const std::vector<double> &w,
                         std::pair<int, int> sz,
                         const std::vector<std::vector<double>> &R,
                         int first_row,
                         int col);

double norm_lower_tri_col(std::pair<int, int> sz,
                          const std::vector<std::vector<double>> &R,
                          int col);

double norm1(std::pair<int, int> sz, const std::vector<std::vector<double>> &B);
double normInf(std::pair<int, int> sz, const std::vector<std::vector<double>> &B);

void tril(std::pair<int, int> sz, std::vector<std::vector<double>> &X);

std::vector<std::vector<double>> triu(std::pair<int, int> sz, std::vector<std::vector<double>> &X);

std::vector<double> gemv(std::pair<int, int> sz,
                         double alpha,
                         const std::vector<std::vector<double>> &B,
                         const std::vector<double> &w,
                         int start);

std::vector<double> gemvt(std::pair<int, int> sz,
                          double alpha,
                          const std::vector<std::vector<double>> &B,
                          const std::vector<double> &w,
                          int start);

void trtrs( int /*unused*/,  const std::vector<std::vector<double>> &R, std::vector<double> &rhs);


void axpy(double alpha, std::pair<int, int> sz,
             const std::vector<std::vector<double>> &X,
             std::vector<std::vector<double>> &Y);

void axpy( const std::vector<double>& x,  double a, std::vector<double>& y);

void gemm_nn(double alpha, std::pair<int, int> szLeft,
             const std::vector<std::vector<double>> &Left,
             std::pair<int, int> szRight,
             const std::vector<std::vector<double>> &Right,
             double beta,
             std::vector<std::vector<double>> &P);

void gemm_tn(double alpha, std::pair<int, int> szLeft,
             const std::vector<std::vector<double>> &Left,
             std::pair<int, int> szRight,
             const std::vector<std::vector<double>> &Right,
             double beta,
             std::vector<std::vector<double>> &P);

void gerSpecial(std::pair<int, int> sz,  // m-by-n,  a += x y'
         std::vector<double> &x,
         int /*unused*/,  // size_x, >= 1 + (m - 1)*incx;
         std::vector<double> &y,
         int /*unused*/,  // size_y >= 1 + (n - 1)*incy;
         std::vector<std::vector<double>> &A, int start);

void ger(std::pair<int, int> sz,  // m-by-n,  a += x y'
         std::vector<double> &x,
         int /*unused*/,  // size_x, >= 1 + (m - 1)*incx;
         std::vector<double> &y,
         int /*unused*/,  // size_y >= 1 + (n - 1)*incy;
         std::vector<std::vector<double>> &A,
         int start);

void print_matrix(const std::string &name, const std::vector<std::vector<double>> &R);
void print_square_matrix_col( std::vector<double> m);

std::vector<double> square_matrix_transpose(  const std::vector<double>& B);
double normInf( const std::vector<double> x);
void scale_vector( std::vector<double>& v,  double scale);

// return B(:,start:end) * w * alpha
std::vector<double> gemv(std::pair<int, int> sz,
  double alpha,
  const std::vector<double> &B,
  const std::vector<double> &w,
  int start);
