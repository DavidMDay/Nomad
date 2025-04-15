#include <iostream>
#include <vector>
#include <cmath>
namespace {

std::vector<std::vector<double>> order3triangular();
std::vector<std::vector<double>> SO3();
std::vector<std::vector<double>> order3();
std::vector<std::vector<double>> order3rank2();
std::vector<std::vector<double>> fourby3();
std::vector<std::vector<double>> fiveby3();
std::vector<std::vector<double>> sixby3();

std::vector<std::vector<double>> order3triangular() {
  std::vector<std::vector<double>> R = {{1, 4, 6}, {0, 2, 5}, {0, 0, 3}};
  return R;
}

std::vector<std::vector<double>> SO3() {
  auto c = cos(M_PI / 3.0);  // 1/2
  auto s = sin(M_PI / 3.0);  // sqrt(3)/2
  std::vector<std::vector<double>> Q = {{0, c, s}, {1, 0, 0}, {0, s, -c}};
  return Q;
}

std::vector<std::vector<double>> order3() {
  std::vector<std::vector<double>> B = {{2.0, -1.0, 0.0}, {-1.0, 2.0, -1.0}, {0.0, -1.0, 2.0}};
  return B;
}
std::vector<std::vector<double>> order3rank2() {
  std::vector<std::vector<double>> S = {{2.0, -1.0, -1.0}, {-1.0, 2.0, -1.0}, {-1.0, -1.0, 2.0}};
  return S;
}

std::vector<std::vector<double>> fourby3() {
  constexpr double e = 1.e-6;
  std::vector<std::vector<double>> S = {
      {2.0, -1.0, -1.0}, {-1.0, 2.0, -1.0}, {-1.0, -1.0, 2.0}, {e, e, e}};
  return S;
}
std::vector<std::vector<double>> fiveby3() {
  auto c = cos(M_PI / 3.0);  // 1/2
  auto s = sin(M_PI / 3.0);  // sqrt(3)/2
  constexpr double e = 1.e-6;
  std::vector<std::vector<double>> B = {{0, c, s}, {0, s, -c}, {e, 0, 0}, {0, e, 0}, {0, 0, e}};
  return B;
}
std::vector<std::vector<double>> sixby3() {
  constexpr double e = 1.e-6;
  std::vector<std::vector<double>> S = {{2.0, -1.0, -1.0}, {-1.0, 2.0, -1.0}, {-1.0, -1.0, 2.0},
                                        {e, 0, 0},         {0, e, 0},         {0, 0, e}};
  return S;
}
}  // namespace

std::vector<std::vector<double>> matrixCollection(int tag) {
  if (tag < 0) {
    tag = -tag;
  }
  tag = tag % 7;
  switch (tag) {
    case 0:
      return order3triangular();
    case 1:
      return SO3();
    case 2:
      return order3();
    case 3:
      return order3rank2();
    case 4:
      return fourby3();
    case 5:
      return fiveby3();
    case 6:
      return sixby3();
    default:
      std::cout << " some matrix bug\n";
  }
  return order3();
}
