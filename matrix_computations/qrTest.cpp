#include "qrTest.hpp"
#include <iostream>
#include <memory>
#include <numeric>
#include <array>
#include <cassert>
#include <cmath>
#include <vector>
#include <utility>
#include "testDecomposition.hpp"  // for qr_test_decomp_parameters
//#include "triuSol.hpp"
#include "QRdecomposition.hpp"  // for QRdecomposition
#include "util.hpp"  // for get_size
#include "matrixCollection.hpp"  // for matrixCollection

bool hqrTestDecomposition() {
  std::string prefix = "Matrix";
  constexpr double threshold = 1.e-12;
  bool pass = true;
  qr_test_decomp_parameters grab_bag;
  for (int tag = 0; tag < 7; tag++) {
    auto A = matrixCollection(tag);
    auto szA = get_size(A);
    auto xA = norm1(szA, A);
    grab_bag.name = prefix + std::to_string(tag);
    grab_bag.tol = threshold * xA;
    std::vector<std::vector<double>> Q1;
    std::vector<std::vector<double>> R1;
    QRdecomposition(A, Q1, R1);
    pass = pass && testDecomposition(Q1, R1, A, grab_bag);
  }
  return pass;
}

bool hqr2Test() {
  std::string prefix = "Compressed";
  constexpr double threshold = 1.e-12;
  bool pass = true;
  qr_test_decomp_parameters grab_bag;
  for (int tag = 0; tag < 7; tag++) {
    auto A = matrixCollection(tag);
    auto sz = get_size(A);
    auto xA = norm1(sz, A);
    grab_bag.name = prefix + std::to_string(tag);
    grab_bag.tol = threshold * xA;
    std::vector<double> tau;
    auto QX = A;
    hqr2(QX, tau);
    auto R1 = triu(sz, QX);
    auto Q1 = std::vector<std::vector<double>>(sz.first, std::vector<double>(sz.second, 0.0));
    for (int i = 0; i < sz.second; i++) {
      Q1[i][i] = 1.0;
    }
    applyQ(QX, tau, Q1);
    auto Id = Q1;
    applyQT(QX, tau, Id);
    for (int i = 0; i < sz.second; i++) {
      Id[i][i] -= 1.0;
    }
    auto residQT = norm1(sz, Id);
    if (residQT > grab_bag.tol) {
      std::cout << "QT test failed" << std::endl;
      pass = false;
    }
    pass = pass && testDecomposition(Q1, R1, A, grab_bag);
  }
  return pass;
}
