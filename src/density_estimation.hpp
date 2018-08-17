#ifndef PROGETTO_DENSITY_ESTIMATION_HPP
#define PROGETTO_DENSITY_ESTIMATION_HPP

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <iterator>
#include <vector>
#include "bspline.hpp"
#include "sandia_rules.hpp"
#include "gauss_points_weights.hpp"

constexpr double tol = 1e-04;

// NOTE: in myDensity::solve() there is the following note:
//        NOTE: QR should do this automatically without giving any warning (CHECK)

/*!
@brief	Parameter class manager
@note These parameters are related to the problem itself, not to each single row.
*/
class parametersManager
{
protected:

  unsigned int k;  /*! Spline degree */
  unsigned int l;  /*! order of derivative in penalization term */
  double alpha;  /*! penalization parameter */

  unsigned int n;   /*! Number of control points */
  unsigned int g;   /*! Number of knots - 2 (number of internal ones) */
  unsigned int G;   /*! Number of knots including additional ones G = g+k+1 */

  std::vector<double> knots; /*! spline knots */
  double u, v;    /*! [u,v] support of splines */

  std::vector<double> xcp; /*! control points */

public:
  parametersManager
  (const unsigned int kk, const unsigned int ll, const double opt_param):
    k(kk), l(ll), alpha(opt_param) {};

  /*!
	@brief Read knots from C-array style and store them in the std::vector<double> member
  */
  void
  readKnots
  (const double * inputKnots, const unsigned int & size);

  /*!
	@brief Read control points from C-array style and store them in the std::vector<double> member
  @param cancel If cancel = j, it does not consider j-th column.
        Useful for cross-validation. Disabled by default.
  */
  void
  readXcp
  (const double * inputXcp, const unsigned int & size, const int & cancel = -1);

  void
  set_alpha
  (const double & opt_param)
  {
    alpha = opt_param;
  }

  unsigned int
  get_G() const { return G; };

  double
  get_u() const { return u; };

  double
  get_v() const { return v;};

  unsigned int
  get_k() const{ return k;};

  unsigned int
  get_n() const{ return n;};
};

/*!
@brief	Parameter class manager
@note These are the whole set of parameters
*/
class myDensity: public parametersManager
{
private:

  Eigen::MatrixXd C;   /*! Collocation matrix - dimensions nxG */
  Eigen::MatrixXd M;   /*! dimensions GxG */

  Eigen::SparseMatrix<double> S;
  Eigen::SparseMatrix<double> DK; /*! dimensions GxG */
  Eigen::MatrixXd P; /*! matrix of the problem to solve - dimensions GxG */

  Eigen::VectorXd p; /*! constant term vector of the problem to solve - dimensions G */
  Eigen::VectorXd c; /*! solution of the problem: c = P^(-)p - dimensions G */
  Eigen::VectorXd b; /*! B-spline coefficients - dimensions G */

  Eigen::VectorXd weights; /*! Weights associated to data, all set to value 1.
                              More in general they allow to take in account the various data accuracy */
  std::vector<double> lambda;  /*! extended vector of knots - with extra ones
                                 dimension: g + 2k + 2 = G + k + 1 */
  std::vector<double> lambda_der; /*! extended vector of knots of the l-th derivative */

  void
  fill_C
  (const std::vector<double>& cp);

  void
  fill_M
  ();

  void
  fill_DK
  ();

  /*!  Compute the S_l matrix for the penalization term */
  void
  fill_S
  ();

  void
  set_lambda
  (const std::vector<double>& knots);

  void
  set_lambda_der
  (const std::vector<double> & knots);

public:

  myDensity
  (const parametersManager & input): parametersManager(input)
  {};

  /*!
  	@brief Compute the matrixes of the problem
  */
  void
  set_matrix
  ();

  /*!
  	@brief Compute the matrix P of the problem
  */
  void
  set_system
  ();

  /*!
  	@brief Compute the value of the functional to be minimized.
  	@param ycp Value of data in the control points.
  */
  double
  eval_J
  (const std::vector<double>& ycp);

  /*!
  	@brief Compute the solution of the linear system.
  	@param bspline Vector of Bspline coefficients computed from the solution of the system.
  	@param ycp Value of data in the control points.
  */
  void
  solve
  (Eigen::Block<Eigen::Matrix<double, -1, -1>, 1, -1, false> bspline, const std::vector<double>& ycp);

  void
  print_all
  () const;

  void
  print_sol
  () const;

  std::vector<double>
  get_lambda() const { return lambda; };

};

#endif //PROGETTO_DENSITY_ESTIMATION_HPP
