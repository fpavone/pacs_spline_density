#ifndef PROGETTO_DENSITY_ESTIMATION_HPP
#define PROGETTO_DENSITY_ESTIMATION_HPP

#include <Eigen/Eigenvalues>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>
#include "bspline.hpp"
#include "sandia_rules.hpp"
#include "gauss_points_weights.hpp"
#include "find_type.hpp"

constexpr double tol = 1e-04;

// NOTE: maybe useful a diagonal matrix W instead of weights.asDiagonal

// NOTE: in myDensity::solve() there is the following note:
//        NOTE: QR should do this automatically without giving any warning (CHECK)


class myParameters {
protected:

  unsigned int k;  // Spline degree
  unsigned int l;  // order of derivative in penalization term
  double alpha;  // penalization parameter

  unsigned int n;   // Number of control points
  unsigned int g;   // Number knots - 2
  unsigned int G;   // Number of knots including additional ones

  std::vector<double> knots; // spline knots
  double u, v;    // [u,v] support of spline

  std::vector<double> xcp;


  // NOTE: generalization for different data input
  unsigned int nclasses;
  double min, max;
  std::vector<double> intervals;

public:
  myParameters
  (const unsigned int kk, const unsigned int ll, const double opt_param):
    k(kk), l(ll), alpha(opt_param) {};

  void
  createKnots
  (const unsigned int & size, const double & uu, const double & vv);

  void
  getKnots
  (const double * inputKnots, const unsigned int & size); // read knots by copy

  void
  getXcp
  (const double * inputXcp, const unsigned int & size);

  unsigned int getG() const{return G;}; //inline ??

  double get_u() const { return u; };

  double get_v() const { return v;};

  unsigned int get_k() const{ return k;};
};


class myDensity: public myParameters {
// G = G + k + 1
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

    void fill_C
      (const std::vector<double>& cp);

    void fill_M
      (); // it uses lambda
      // NOTE: better to void fill_M() and use lambda member as in fill_DK?

    void fill_DK
      (); // it uses lambda

    // Compute the S_l matrix for the penalization term
    void fill_S
      (); // it uses lambda

    void set_lambda
      (const std::vector<double>& knots);

    void set_lambda_der
      (const std::vector<double> & knots);

public:

    myDensity
    (const myParameters & input): myParameters(input) {};

    void
    set_matrix
    ();

    void
    set_density
    (const std::vector<double>& ycp);

    void
    solve(Eigen::Block<Eigen::Matrix<double, -1, -1>, 1, -1, false> bspline);

    void
    print_all
    () const;

    void print_sol() const
    {
      std::cout << "\n Matrix P: " << '\n' << P << '\n';
      // std::cout << "SOLUTION c = P^(-)p:" << '\n' << c << '\n';
       Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(P);
      std::cout << "\n P positive definite? \n " << es.eigenvalues() << std::endl;
      std::cout << "\n Eigenvectors: \n " << es.eigenvectors() << std::endl;

      // double relative_error = (P*c - p).norm() / p.norm(); // norm() is L2 norm
      // std::cout << "The relative error is:\n" << relative_error << std::endl;
    //  std::cout << "B-SPLINE COEFFICIENTS b = DKc" << '\n' << b << '\n';
    //  std::cout << "PROVA P*c = p..?" << '\n' << Eigen::VectorXd(P*c) << '\n';
    };

    std::vector<double> get_lambda() const{ return lambda; };

    // void save_matrix() const
    // {
    //   P.makeCompressed();
    //   Eigen::saveMarket(P, "density.mtx");
    //   Eigen::saveMarketVector(p, "density_b.mtx");
    // };
};

#endif //PROGETTO_DENSITY_ESTIMATION_HPP
