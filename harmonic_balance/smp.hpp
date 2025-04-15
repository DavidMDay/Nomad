#pragma once
#include <vector>
// column major order,  (i,j) i + nr*j
std::vector<double> square_matrix_product_col(const std::vector<double>& Lvec,
                                              const std::vector<double>& Rvec);