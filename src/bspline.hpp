#ifndef HH_BSPLINE_HH
#define HH_BSPLINE_HH
#include <vector>
#include <Eigen/Dense>
namespace bspline{
  using vect = std::vector<double>;

  int
  findspan (int n, int p, double t, const vect &U);

  void
  basisfun (int i, double t, int p, const vect &U, Eigen::ArrayXd &N);
}

#endif
