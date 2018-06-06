#ifndef PROGETTO_DENSITY_ESTIMATION_HPP
#define PROGETTO_DENSITY_ESTIMATION_HPP

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <iterator>
// #include <string>
#include <vector>
#include "bspline.hpp"
#include "sandia_rules.hpp"
#include "gauss_points_weights.hpp"

constexpr double tol = 1e-04;

// NOTE: maybe useful a diagonal matrix W instead of weights.asDiagonal

// NOTE: in myDensity::solve() there is the following note:
//        NOTE: QR should do this automatically without giving any warning (CHECK)


class myParameters
{
protected:

  unsigned int k;  // Spline degree
  unsigned int l;  // order of derivative in penalization term
  double alpha;  // penalization parameter

  unsigned int n;   // Number of control points
  unsigned int g;   // Number knots - 2
  unsigned int G;   // Number of knots including additional ones G = g+k+1

  std::vector<double> knots; // spline knots
  double u, v;    // [u,v] support of spline

  std::vector<double> xcp;

public:
  myParameters
  (const unsigned int kk, const unsigned int ll, const double opt_param):
    k(kk), l(ll), alpha(opt_param) {};

  void
  readKnots
  (const double * inputKnots, const unsigned int & size); // read knots by copy

  void
  readXcp
  (const double * inputXcp, const unsigned int & size);

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
};


class myDensity: public myParameters
{
private:

  Eigen::MatrixXd C;   // Collocation matrix - nxG
  Eigen::MatrixXd M;   // GxG

  Eigen::SparseMatrix<double> S;
  Eigen::SparseMatrix<double> DK; // GxG
  Eigen::MatrixXd P; // matrix of the problem we have to solve - GxG

  Eigen::VectorXd p; // known vector of the problem we have to solve - G
  Eigen::VectorXd c; // solution of the problem: c = P^(-)p - G
  Eigen::VectorXd b; // B-spline coefficients - G

  Eigen::VectorXd weights;
  std::vector<double> lambda;  // extended vector of knots - with extra ones
                               // dimension: g + 2k + 2 = G + k + 1
  std::vector<double> lambda_der;

  void
  fill_C
  (const std::vector<double>& cp);

  void
  fill_M
  ();

  void
  fill_DK
  ();

  // Compute the S_l matrix for the penalization term
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
  (const myParameters & input): myParameters(input)
  {};

  void
  set_matrix
  ();

  void
  set_system
  ();

  double
  eval_J
  (const std::vector<double>& ycp);

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
