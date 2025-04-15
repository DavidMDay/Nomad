#pragma once

#include <string>
struct qr_test_decomp_parameters {
  std::string name;
  double tol = 0.0;
  int tag = 0;
};

bool testDecomposition(const std::vector<std::vector<double>>& Q1,
                       const std::vector<std::vector<double>>& R1,
                       const std::vector<std::vector<double>>& A,
                       qr_test_decomp_parameters grab_bag);
