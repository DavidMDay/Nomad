#include <cassert>
#include <vector>
#include <iostream>
#include "stridedInnerProduct.hpp"

std::pair<int, int> test_sip(const std::vector<int> &left,
                             int left_stride,
                             const std::vector<int> &right,
                             int right_stride) {
  auto safe =
      left_stride * static_cast<int>(right.size()) - right_stride * static_cast<int>(left.size());
  if (safe < 0) {
    std::cout << " test_sip : memory test fails by " << safe << std::endl;
    return {1, 0};
  }
  int correct =
      ip_with_stride(left.begin(), left.end(), left_stride, right.begin(), right_stride, 0);
  if (right_stride != 1 && left_stride != 1) {
    std::cout << " test_sip : one stride should be 1\n";
    std::cout << "but the strides are " << left_stride << "  " << right_stride << std::endl;
    return {0, 1};
  }
  int result =
      (left_stride == 1 ? ip_with_strideR(left.begin(), left.end(), right.begin(), right_stride, 0)
                        : ip_with_strideL(left.begin(), left.end(), left_stride, right.begin(), 0));
  return {correct, result};
}

bool testRightStridedInnerProduct() {
  std::vector<int> left = {1, 2, 3, 4, 5};
  int left_stride = 1;
  std::vector<int> right = {-1, -2, -3, -4, -5, -6, -7, -8, -9, -10};
  int right_stride = 2;
  std::pair<int, int> cr = test_sip(left, left_stride, right, right_stride);
  return (cr.second == cr.first);
}

bool testLeftStridedInnerProduct() {
  std::vector<int> left = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  int left_stride = 2;
  std::vector<int> right = {-1, -2, -3, -4, -5};
  int right_stride = 1;
  std::pair<int, int> cr = test_sip(left, left_stride, right, right_stride);
  return (cr.second == cr.first);
}
