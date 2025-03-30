
#include <iostream>
#include <vector>
#include <utility>
#include <cmath>
#include <cassert>
#include "triuSol.hpp"
#include "getSquareUpperTri.hpp"
#include "util.hpp"  // for axpy

bool triSolTest1x1() {
// use QR to solve a linear system: test triangular linear system solver
constexpr double tol = 1.e-14;
double d = 1.1;
double ud = 0.0;
std::pair<int, int> sz1 = {1,1};
std::vector<std::vector<double>>  R1 = getSquareUpperTri(sz1, d, ud);
std::vector<std::vector<double>> invR1 =  getSquareUpperTriInv(sz1, d, ud);
std::vector<double>  c1 = diag (sz1, R1);
if ( std::abs( d - c1[0] )> tol) {
    std::cout << " diag test failed " << d << "   "  << c1[0] << "\n";
    return false;
}
double min_diag  = min_abs_diag( sz1, R1 );
if ( std::abs( min_diag * invR1[0][0]-1.0 )> tol ) {
    std::cout << " min diag test failed " << d << "   "  << min_diag << "\n";
    return false;
}

std::vector<double> sol1 = {1.0};
std::vector<double> rhs1 = apply_triu(sz1, R1, sol1);// d
std::vector<double> s1 = apply_inv_triu(sz1, R1, rhs1);// 1
std::vector<double> r1 = apply_triu(sz1, invR1, sol1);
return  (  std::abs( s1[0] - sol1[0]) <  tol  && 
    std::abs( r1[0]*d - 1.0 ) < tol );
}

bool triSolTestOrder(int n)
{
    // use QR to solve a linear system: test triangular linear system solver
    constexpr double tol = 1.e-14;
    double d = 1.1;
    double ud = 0.0;
    std::pair<int, int> sz = {n, n};
    std::vector<std::vector<double>> R = getSquareUpperTri(sz, d, ud);
    auto sol = std::vector<double>(n, 1.0);
    std::vector<double> rhs = apply_triu(sz, R, sol);
    std::vector<double> x = apply_inv_triu(sz, R, rhs);
    axpy(sol, -1.0, x);
    auto err = normInf(x);
    if (err < tol)
    {
        return true;
    }
    // rhs = R * sol
    //  rhs = R * x
    std::vector<std::vector<double>> invR = getSquareUpperTriInv(sz, d, ud);
    std::vector<double> s = apply_triu(sz, invR, rhs);
    axpy(sol, -1.0, s);
    auto u_err = normInf(s);
    if (u_err > tol)
    {
        std::cout << n << " apply_triu test failed " << u_err << "\n";
    }
    x = apply_inv_triu(sz, R, rhs);
    s = apply_triu(sz, invR, rhs);
    axpy(x, -1.0, s);
    double x_err = normInf(s);
    if (x_err > tol)
    {
        std::cout << n << " apply_triu test failed " << x_err << "\n";
    }
    else
    {
        std::cout << n << " triSol test failed " << err << "\n";
    }
    return false;
}

bool triSolTest() {
    if (!triSolTest1x1() ) {
        std::cout << " 1 by 1 tri sol test failed\n";
        return false;
    }
    for (int n = 2; n < 10; n++) {
        if ( !triSolTestOrder(n) ) {
            return false;
        }
    }
    return true;
}

