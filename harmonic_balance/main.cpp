#include <cstdlib>
#include <iostream>
#include "Tests.hpp"

int main() {
  if (TestOddDFT()==EXIT_FAILURE) {
    std::cout << "DFT Test failure: odd\n";
    return EXIT_FAILURE;
  }
  if (TestEvenDFT()==EXIT_FAILURE) {
    std::cout << "DFT Test failure: even\n";
    return EXIT_FAILURE;
  }
  if (TestLiuOde() == EXIT_FAILURE) {
    std::cout << "Liu Ode failure: test_left\n";
    return EXIT_FAILURE;
  }
  if (TestLightlyDampedOde() == EXIT_FAILURE)
  {
    std::cout << "Space Time Trans Lightly damped Ode failure\n";
    return EXIT_FAILURE;
  }
  if (TestSpaceTimeTransOde() == EXIT_FAILURE)
  {
    std::cout << "Space Time Trans Ode failure\n";
    return EXIT_FAILURE;
  }

  if (TestIncreasedOmega() == EXIT_FAILURE)
  {
    std::cout << "Space Time Trans Increased Omega failure\n";
    return EXIT_FAILURE;
  }
  if (TestIncreasedOmegaLiu() == EXIT_FAILURE) {
    std::cout << "Space Time Trans failure\n";
    return EXIT_FAILURE;
  }
  if (TestLinearized() == EXIT_FAILURE) {
    std::cout << "Space Time Trans Linear ODE failure\n";
    return EXIT_FAILURE;
  }
  if (TestOperators() == EXIT_FAILURE)
  {
    std::cout << "Operators failure\n";
    return EXIT_FAILURE;
  }
  std::cout << "pass\n";
  return EXIT_SUCCESS;
}

