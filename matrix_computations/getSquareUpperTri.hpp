#pragma once
#include <cassert>
#include <cmath>
#include <utility>
#include <vector>
// upper triangular square matrix

// x  u   u   u
// 0  x   u   u
// 0  0   x   u
// 0  0   0   x
std::vector<std::vector<double>> getSquareUpperTri(std::pair<int, int> sz,
                                                   double diag,
                                                   double upper_diag);
// x  u   u   u      y   -v1   v2 -v3
// 0  x   u   u      0    y   -v1  v2
// 0  0   x   u      0    0    y  -v1
// 0  0   0   x      0    0    0    y = I
// v0 = 1/x
// v1 = v0 * u * v0;
// v2 = v1 *(u*v0-1)
// v3 = v2 *(u*v0-1)
std::vector<std::vector<double>> getSquareUpperTriInv(std::pair<int, int> sz,
                                                      double diag,
                                                      double upper_diag);
