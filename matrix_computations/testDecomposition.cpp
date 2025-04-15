#include <iostream>
#include <memory>
#include <numeric>
#include <array>
#include <cmath>
#include "util.hpp" // for print_matrix
#include "testDecomposition.hpp"

namespace {
double norm2(std::array<double, 3> x) {
  return std::sqrt(std::inner_product(x.begin(), x.end(), x.begin(), 0.0));
}
void print_decomposition(const std::string& name,
                         const std::vector<std::vector<double>>& Q,
                         const std::vector<std::vector<double>>& R,
                         const std::vector<std::vector<double>>& X) {
  std::string nq = "Orthogonal";
  print_matrix(nq, Q);
  std::string ut = "Triangular";
  print_matrix(ut, R);
  print_matrix(name, X);
}

bool checkSize(std::pair<int, int> szQ, std::pair<int, int> szR, std::pair<int, int> szX) {
  if (szX.first != szQ.first || szX.second != szR.second) {
    return false;
  }
  if (szQ.second == szR.second && szR.second == szR.first) {
    return true;
  }
  return (szQ.second == szQ.first && szR.second == szX.second);
}

// norm(X - Q*R); norm(Q'*Q - eye(p)); norm(R(1:p,:)- triu(R))];
std::array<double, 3> checkValues(const std::vector<std::vector<double>>& Q,
                                  const std::vector<std::vector<double>>& R,
                                  const std::vector<std::vector<double>>& X) {
  auto szQ = get_size(Q);
  auto szR = get_size(R);
  auto P = X;
  double alpha = -1.0;
  double beta = 1.0;
  gemm_nn(alpha, szQ, Q, szR, R, beta, P);
  std::array<double, 3> residual{};
  std::pair<int, int> szP = {szQ.first, szR.second};
  residual[0] = norm1(szP, P);

  auto M = std::vector<std::vector<double>>(szQ.second, std::vector<double>(szQ.second, 0.0));
  for (int i = 0; i < szQ.second; i++) {
    M[i][i] = 1.0;
  }
  gemm_tn(alpha, szQ, Q, szQ, Q, beta, M);
  std::pair<int, int> szM = {szQ.second, szQ.second};
  residual[1] = norm1(szM, M);
  auto trilR = R;
  tril(szR, trilR);
  residual[2] = normInf(szR, trilR);
  return residual;
}
}  // namespace
bool testDecomposition(const std::vector<std::vector<double>>& Q1,
                       const std::vector<std::vector<double>>& R1,
                       const std::vector<std::vector<double>>& A,
                       qr_test_decomp_parameters grab_bag) {
  auto szQ1 = get_size(Q1);
  auto szR1 = get_size(R1);
  bool pass = true;
  auto szA = get_size(A);
  if (checkSize(szQ1, szR1, szA)) {
    auto residual = checkValues(Q1, R1, A);
    if (norm2(residual) > grab_bag.tol) {
      if (residual[0] > grab_bag.tol) {
        std::cout << "testDecomposition: norm X - QR " << residual[0] << std::endl;
      }
      if (residual[1] > grab_bag.tol) {
        std::cout << "testDecomposition: norm Qt Q - I " << residual[1] << std::endl;
      }
      if (residual[2] > grab_bag.tol) {
        std::cout << "testDecomposition: norm tril R " << residual[2] << std::endl;
      }
      grab_bag.name += "Values";
      pass = false;
      print_decomposition(grab_bag.name, Q1, R1, A);
    }
  } else {
    grab_bag.name += "Sizes";
    pass = false;
    print_decomposition(grab_bag.name, Q1, R1, A);
  }
  return pass;
}
