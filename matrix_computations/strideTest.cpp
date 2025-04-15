#include <cstdlib>
#include <iostream>
bool testLeftStridedInnerProduct();
bool testRightStridedInnerProduct();

int main() {
  if (!testRightStridedInnerProduct()) {
    std::cout << "StridedInnerProduct failure: test_right\n";
    return EXIT_FAILURE;
  }

  if (!testLeftStridedInnerProduct()) {
    std::cout << "StridedInnerProduct failure: test_left\n";
    return EXIT_FAILURE;
  }

  std::cout << "pass\n";
  return EXIT_SUCCESS;
}
