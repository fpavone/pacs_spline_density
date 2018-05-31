#ifndef HH_BSPLINE_HH
#define HH_BSPLINE_HH
#include <vector>
#include <Eigen/Dense>
namespace bspline{
  using vect = std::vector<double>;

  /*!
	@brief	Find the knot span of the parametric point u.

	@note	This is NOT	Algorithm A2.1 from 'The NURBS BOOK' pg68
			as that algorithm only works for nonperiodic
			knot vectors, nonetheless the results should
			be EXACTLY the same if U is nonperiodic

	@param n Number of control points - 1
	@param p Spline degree
	@param t Parametric point
	@param U Knot sequence
	@return	Knot span

	@todo	This implementation has linear, rather than log complexity
  */
  unsigned int
  findspan
  (int p, double t, const vect &U);

  /*!
  	@brief	Compute the functions of the basis

  	@note	Algorithm A2.2 from 'The NURBS BOOK' pg70.

  	@param i (Input) Knot span (from findspan())
  	@param t (Input) Parametric point
  	@param p (Input) Spline degree
  	@param U (Input) Knot sequence
  	@param N (Output) Vector of the functions of the basis (p+1 dimensional)
  */
  void
  basisfun
  (unsigned int i, double t, int p, const vect &U, Eigen::ArrayXd &N);
}

#endif
