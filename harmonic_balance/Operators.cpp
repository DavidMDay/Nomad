
#include "Operators.hpp"
#include <vector>
#include <cmath>
#include <numeric>                          // for accumulate
#include <complex>                          // for complex, operator*, imag
// column major order,  (i,j) i + nr*j
#include "util.hpp"
#include "smp.hpp" // for square_matrix_product

std::vector<double> getA(int order){
    if (order < 0) {
        std::vector<double> X;
        return X;
    }
    int n = 2*order+1;
    std::vector<double> A(n*n,0.0);
    if (order == 0){
        return A;
    }

    for (int i=0;i<order;i++) {
        int col = 2*i + 1;
        int row = col+1;
        int ko = row + n*col;
        A[ko] = static_cast<double>(i+1);
        row = 2*i + 1;
        col = row+1;
        int k = row + n*col;
        A[k] = -A[ko];
    }
    return A;
}
std::vector<double> getE(int order)
{
	if (order < 0) {
        std::vector<double> X;
        return X;
    }
    int n = 2 * order + 1;
    std::vector<double> E(n * n, 0.0);
    using namespace std::complex_literals;
    double r = 2.0 * M_PI / static_cast<double>(n); // t_1 = r
    std::complex<double> zero = 0.0;
    std::complex<double> y = zero;
    const double scale = 2.0/static_cast<double>(n);
    for (int col0=0;col0<n;col0++) {
        int ko = n*col0;
        E[ko] = 0.5*scale;
    }
    for (int col = 0; col < n; col++)
    {
        for (int i = 0; i < order; i++)
        {
            double angle = r * static_cast<double>((i+1) * col);
            y = std::exp(1i * angle);
            int row = 2*i + 1;
            int k = row + n*col;
            E[k]=real(y)*scale;
            E[k+1]=imag(y)*scale; // (row,col+1)
        }
    }
    return E;
}
std::vector<double> getAsquared(int order){
    if (order < 0) {
        std::vector<double> X;
        return X;
    }
    int n = 2*order+1;
    std::vector<double> A2(n*n,0.0);
    for (int i=0;i<order;i++) {
        int row = 2*i + 1;
        int col = row;
        int k = row + n*col;
        A2[k] = - static_cast<double>( (i+1)*(i+1) );
        A2[k+n+1] = A2[k];
    }
    return A2;
}
std::vector<double> getInvE(int order){
    if (order < 0) {
        std::vector<double> X;
        return X;
    }
    int n = 2*order+1;
    std::vector<double> invE(n*n,0.0);
    for (int row=0;row<n;row++) {
        invE[row] = 1.0;
    }
    double r = 2.0 * M_PI / static_cast<double>(n); // t_1 = r
    using namespace std::complex_literals;
    std::complex<double> zero = 0.0;
    std::complex<double> y = zero;
    for (int i = 0; i < order; i++)
    {
        int col = 2*i + 1;
        for (int row = 0; row < n; row++)
        {
            double angle = r * static_cast<double>((i+1) * row);
            y = std::exp(1i * angle);
            int k = row + n*col;
            invE[k]=real(y);
            invE[k+n]=imag(y); // (row,col+1)
        }
    }
    return invE;
}
std::vector<double> getD(int order)
{
    std::vector<double> A = getA(order);
    std::vector<double> E = getE(order);
    std::vector<double> AE = square_matrix_product_col(A, E);
    std::vector<double> invE = getInvE(order);
    std::vector<double> D = square_matrix_product_col(invE, AE);
    return D;
}
std::vector<double> getDsquared(int order)
{
    std::vector<double> D = getD(order);
    std::vector<double> D2 = square_matrix_product_col(D, D);
    return D2;
}
std::vector<double> getDuffingR(int order)
{
    int n = 2 * order + 1;
    std::vector<double> R(n * n, 0.0);
    std::vector<double> u(n, 1.0); // input
    std::vector<double> E = getE(order);

    std::vector<double> u3(n, 0);
    for (size_t i = 0; i < u.size(); i++)
    {
        u3[i] = u[i] * u[i] * u[i];
    }
    return R;
}

std::vector<double> getIdentity(int order){
    if (order < 0) {
        std::vector<double> O;
        return O;
    }
    int n = 2*order+1;
    std::vector<double> Id(n*n,0.0);
    for (int row=0;row<n;row++) {
        Id[row*(n+1)] = 1.0;
    }
    return Id;
}