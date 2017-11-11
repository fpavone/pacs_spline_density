#ifndef HH_CREATE_MATRIX_HH
#define HH_CREATE_MATRIX_HH
#include <Eigen/Dense>
#include <vector>
#include "bspline.hpp"

// g = number of knots in the interval (u,v). g+2 knots includinf the end points. g >= 0
// k = degree of the splines. k > 0
// G=g+k+1 =  dimension of the vector space of polynomial splines of degree k defined
//          on a k on thfinite interval [u, v] with the sequence of knots  Deltaλ
// b=(b_−k,...,b_g) = vector of B-splinec oefficients of sk(x)
// Bk+1(x),i = −k,...,g = B-splines of degree k (basis)
// cp = vector of control points. Dimension n
template <int n, int G>
void fill_C(Eigen::Matrix<double,n,G> &C, std::vector<double> const &cp, int k, std::vector<double> const &knots) {
  double t = 0;
  Eigen::ArrayXd N = Eigen::ArrayXd::Constant(G,0.0);
  for (unsigned int i = 0; i <= n; i++) {
    t = cp[i];
    int fs = bspline::findspan(n,k,t,knots);
    bspline::basisfun (fs, t, n, knots, N);
    for (unsigned int j = 0; j < G; j++) {
      C(i,j) = N(j);
    }
  }
}
#endif
