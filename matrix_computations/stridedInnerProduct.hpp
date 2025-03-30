#pragma once
#include <iostream>
#include <vector>
#include <numeric>

template <typename LeftIter, typename RightIter, typename T>
T ip_with_stride(LeftIter left_first, LeftIter left_last, int left_stride,
                 RightIter right_first, int right_stride, T init) {
  while (std::distance(left_first, left_last) > 0) {
    init += *left_first * *right_first;
    std::advance(left_first, left_stride);
    std::advance(right_first, right_stride);
  }
  return init;
}

template <typename LeftIter, typename RightIter, typename T>
T ip_with_strideR(LeftIter left_first, LeftIter left_last,
                  RightIter right_first, int right_stride, T init) {
  while (std::distance(left_first, left_last) > 0) {
    init += *left_first * *right_first;
    std::advance(left_first, 1);
    std::advance(right_first, right_stride);
  }
  return init;
}

template <typename LeftIter, typename RightIter, typename T>
T ip_with_strideL(LeftIter left_first, LeftIter left_last, int left_stride,
                  RightIter right_first, T init) {
  while (std::distance(left_first, left_last) > 0) {
    init += *left_first * *right_first;
    std::advance(left_first, left_stride);
    std::advance(right_first, 1);
  }
  return init;
}
