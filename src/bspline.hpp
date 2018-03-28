#ifndef HH_BSPLINE_HH
#define HH_BSPLINE_HH
#include <vector>
#include <Eigen/Dense>
namespace bspline{
  using vect = std::vector<double>;

    /*
    p Spline degree
    t Parametric point
    U Knot sequence  int
    */
  unsigned int
  findspan (int p, double t, const vect &U);
    /*
    i (Input) Knot span (from findspan())
    t (Input) Parametric point
    p (Input) Spline degree
    U (Input) Knot sequence
    N (Output) Vector of the functions of the basis (p+1 dimensional)
    */
  void
  basisfun (unsigned int i, double t, int p, const vect &U, Eigen::ArrayXd &N);
}

#endif
