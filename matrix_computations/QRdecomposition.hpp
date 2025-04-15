#pragma once
#include <vector>

// m-by-n A, m>=n, Q square,  m-by-n R
void QRdecomposition(const std::vector<std::vector<double>> &A,
                     std::vector<std::vector<double>> &Q,
                     std::vector<std::vector<double>> &R);

// m-by-n A, m>=n, reduced in place to a compressed representation
// of both Q and R.
void hqr2(std::vector<std::vector<double>> &A, std::vector<double> &tau);

void applyQ(const std::vector<std::vector<double>> &QR,
            const std::vector<double> &tau,
            std::vector<std::vector<double>> &X);

void applyQT(const std::vector<std::vector<double>> &QR,
             const std::vector<double> &tau,
             std::vector<std::vector<double>> &X);
