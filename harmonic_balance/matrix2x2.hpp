#pragma once
#include <vector>

std::vector<double> mat_vec_2x2(const std::vector<double> &x, const std::vector<double> &matrix);

std::vector<double> solve_2by2(const std::vector<double> &rhs, const std::vector<double> &matrix);

std::vector<double> residual_2by2(const std::vector<double> &rhs,
                                  const std::vector<double> &matrix,
                                  const std::vector<double> &sol);

double determinant2x2(const std::vector<double> &matrix);
