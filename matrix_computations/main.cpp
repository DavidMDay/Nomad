#include <cstdlib>
#include <iostream>
#include "qrTest.hpp"

int main() {
  if (!triSolTest()) {
    std::cout << "upper triangular matrix solve failure\n";
    return EXIT_FAILURE;
  }
  if (!testRightStridedInnerProduct()) {
    std::cout << "StridedInnerProduct failure: test_right\n";
    return EXIT_FAILURE;
  }

  if (!testLeftStridedInnerProduct()) {
    std::cout << "StridedInnerProduct failure: test_left\n";
    return EXIT_FAILURE;
  }

  if (!hqrTestDecomposition()) {
    std::cout << "full QR failure\n";
    return EXIT_FAILURE;
  }

  if (!hqr2Test()) {
    std::cout << "compressed QR failure\n";
    return EXIT_FAILURE;
  }
  std::cout << "pass\n";
  return EXIT_SUCCESS;
}
