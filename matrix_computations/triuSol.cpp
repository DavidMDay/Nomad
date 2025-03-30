#include "triuSol.hpp"

#include <vector>
#include <utility>
#include <cmath>
#include <cassert>
// upper triangular square matrix
// min_abs_diag returns the minimum absolute value of the diagonal elements
// of an upper triangular matrix. If the matrix is empty, it returns 0.0.
// The function takes a pair of integers representing the size of the matrix
// and a 2D vector representing the matrix itself.

double min_abs_diag(const std::pair<int, int> sz,
                          const std::vector<std::vector<double>> &A){
auto d = diag(sz,A);
if (d.empty() ) {
       return 0.0;	
}
auto min = d[0]; 
for (const auto& x : d) { 
	 if ( std::abs(x) < std::abs( min ) ) {
		 min = x;
	 }
}
return min;
}

std::vector<double>
diag(const std::pair<int, int> sz,
                          const std::vector<std::vector<double>> &A) {
  auto s = static_cast<int>(std::min( sz.first, sz.second ));
  std::vector<double> d(s,0.0);
  for (int row = 0; row < s; row++) {
	   d[row] = A[row][row];
  }
  return d;
}

// n-by-n upper triangular   R    R is triangular and square
std::vector<double>
apply_triu(std::pair<int, int> sz, const std::vector<std::vector<double>> &R,
   const std::vector<double> &x) {
  auto n = sz.second;
  assert( sz.first == n );
  auto b = std::vector<double>(n, 0.0);
  for (int row = 0; row < n ; row++) {
    for (int col = row; col < n; col++) {
      b[row] += R[row][col] * x[col];
    }
  }
  return b;
}

// upper triangular square matrix 
std::vector<double>
apply_inv_triu(std::pair<int, int> sz, const std::vector<std::vector<double>> &R,
   const std::vector<double> &rhs) {
  auto n = sz.second;
  assert( sz.first == n );
  if (n==0) {
    auto sol = std::vector<double>(0, 0.0);
    return sol;
  }
  assert( R[0][0] != 0.0 );
  auto sol = std::vector<double>(n, 0.0);
  if (n==1) {
    sol[0] = rhs[0]/R[0][0];
    return sol;
  }
  for (int row = n-1; 0 <= row ; row--) {
    assert( row >= 0 );
    assert( R[row][row] != 0.0 );
    auto b = rhs[row];
    for (int col = row+1; col < n; col++) {
      b -= R[row][col] * sol[col];
    }
    sol[row] = b/R[row][row];
  }
  return sol;
}
