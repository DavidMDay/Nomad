#pragma once

#include <vector>
#include <utility>

// an upper triangluar matrix is square
double min_abs_diag(std::pair<int, int> sz, const std::vector<std::vector<double>> &A);
std::vector<double> diag(std::pair<int, int> sz, const std::vector<std::vector<double>> &A);
std::vector<double> apply_triu(std::pair<int, int> sz,
                               const std::vector<std::vector<double>> &R,
                               const std::vector<double> &x);
std::vector<double> apply_inv_triu(std::pair<int, int> sz,
                                   const std::vector<std::vector<double>> &R,
                                   const std::vector<double> &x);
