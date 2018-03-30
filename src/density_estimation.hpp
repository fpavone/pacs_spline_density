#ifndef PROGETTO_DENSITY_ESTIMATION_HPP
#define PROGETTO_DENSITY_ESTIMATION_HPP

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include "bspline.hpp"
#include "sandia_rules.hpp"

// NOTE: maybe useful a diagonal matrix W instead of weights.asDiagonal

class Density {
// G = G + k + 1
private:
    unsigned int k;   // Spline degree
    unsigned int n;   // Number of control points
    unsigned int G;   // Number of knots including additional ones

    double alpha;  // penalization parameter
    double u, v;    // [u,v] support of spline
    double l;       // order of derivative in penalization term

    Eigen::MatrixXd C;   // Collocation matrix - nxG
    Eigen::MatrixXd M;   // GxG

    Eigen::SparseMatrix<double> S;
    Eigen::SparseMatrix<double> DK; // GxG
    Eigen::SparseMatrix<double> P; // matrix of the problem we have to solve - GxG

    Eigen::VectorXd p; // known vector of the problem we have to solve - G
    Eigen::VectorXd c; // solution of the problem: c = P^(-)p - G
    Eigen::VectorXd b; // B-spline coefficients - G

    Eigen::VectorXd weights;
    std::vector<double> lambda;  // extended vector of knots - with extra ones
                                 // dimension: g + 2k + 2 = G + k + 1

    void fill_C
      (const std::vector<double>& cp, const std::vector<double>& knots);

    void fill_M
      (const std::vector<double>& knots); // it uses lambda
      // NOTE: better to void fill_M() and use lambda member as in fill_DK?

    void fill_DK
      (); // it uses lambda

    // Compute the S_l matrix for the penalization term
    void fill_S
      (); // it uses lambda

    void set_lambda
      (const std::vector<double>& knots);

public:

    Density(const std::vector<double>& knots, const std::vector<double>& xcp,
      const std::vector<double>& ycp, double kk, double g, unsigned int ll, double opt_param = 1.0):
      k(kk), n(xcp.size()), G(g+k+1), u(knots[0]), v(*(knots.end()-1)), l(ll), alpha(opt_param)
    {
std::cout << "fill_C.." << '\n';
      // weights.assign(n,1.0);
      weights = Eigen::VectorXd::Constant(n,1.0);
      set_lambda(knots);
      fill_C(xcp, lambda);
std::cout << C << std::endl;
std::cout << "fill_M.." << '\n';
      fill_M(lambda);
std::cout << M << std::endl;
std::cout << "fill_DK.." << '\n';
      fill_DK();
std::cout << Eigen::MatrixXd(DK) << '\n';
std::cout << "fill_S.." << '\n';
      fill_S();
std::cout << Eigen::MatrixXd(S) << '\n';
      P = (1 / alpha * (DK).transpose() * M * (DK) + (C * DK).transpose() * weights.asDiagonal() * C * DK).sparseView();
      Eigen::VectorXd newycp(ycp.size());
      for (unsigned int i = 0; i < ycp.size() ; ++i) {
          newycp[i] = ycp[i];
      }
      p = DK.transpose()* C.transpose() * weights.asDiagonal() * newycp;
std::cout << "Constructor done - p:" << '\n' << p << '\n';
    }

    // Density(const std::vector<double>& knots, const std::vector<double>& xcp, const std::vector<double>& ycp,
    //   double kk, double g, double opt_param):
    //   alpha(opt_param)
    // {
    //   Density(knots, xcp, ycp, kk,g);
    // }
    //
    // Density(const std::vector<double>& knots, const std::vector<double>& xcp, const std::vector<double>& ycp,
    //    double kk, double g, unsigned int ll)
    // {
    //   // assert(ll<G);
    //   l = ll;
    //   Density(knots, xcp, ycp, kk, g);
    // }
    //
    // Density(const std::vector<double>& knots, const std::vector<double>& xcp, const std::vector<double>& ycp,
    //    double kk, double g, double opt_param, unsigned int ll):
    //   alpha(opt_param)
    // {
    //   // assert(ll<G);
    //   l = ll;
    //   Density(knots, xcp, ycp, kk, g);
    // }

    void print_all() const
    {
        std::cout << "MATRIX C:" << '\n' << C << std::endl;
        std::cout << lambda[0] << std::endl;
        std::cout << "MATRIX M:" << '\n' << M << std::endl;
        std::cout << "MATRIX DK:" << '\n' << Eigen::MatrixXd(DK) << std::endl;
        std::cout << "MATRIX W:" << '\n' << Eigen::MatrixXd(weights.asDiagonal()) << std::endl;
    }

    void solve()
    {   /* NAIVE SOLVER */
      Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
      // Compute the ordering permutation vector from the structural pattern of P
      solver.analyzePattern(P);
      // Compute the numerical factorization
      solver.factorize(P);
      //Use the factors to solve the linear system
      c = solver.solve(p);
      b = DK*c;
    };

    void print_sol() const
    {
      std::cout << "SOLUTION c = P^(-)p:" << '\n' << c << '\n';
      std::cout << "B-SPLINE COEFFICIENTS b = DKc" << '\n' << b << '\n';
    };
};

#endif //PROGETTO_DENSITY_ESTIMATION_HPP
