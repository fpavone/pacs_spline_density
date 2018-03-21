#ifndef PROGETTO_DENSITY_ESTIMATION_HPP
#define PROGETTO_DENSITY_ESTIMATION_HPP

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include "bspline.hpp"
#include "sandia_rules.hpp"

class Density {

private:
    unsigned int k; //Spline degree
    unsigned int n; //Number of control points
    unsigned int G;

    double alpha = 1.0;
    double u, v;  // [u,v] support of spline
    double l;     // order of derivative in penalization term

    Eigen::MatrixXd C;   // nxG
    Eigen::MatrixXd M;    // GxG
    Eigen::MatrixXd W;    // nxn
    Eigen::SparseMatrix<double> S;

    Eigen::SparseMatrix<double> DK; // GxG
    Eigen::SparseMatrix<double> P; // matrix of the problem we have to solve GxG
    Eigen::VectorXd p; // known vector of the problem we have to solve G
    Eigen::VectorXd c; // c = P^(-)p G
    Eigen::VectorXd b; // B-spline coefficients G

    std::vector<double> weights;
    std::vector<double> lambda;  // extended vector of knots - with extra ones
                                // dimension: g + 2k + 2 = G + k + 1 ??

    void fill_C
      (const std::vector<double>& cp, const std::vector<double>& knots);

    void fill_M
      (const std::vector<double>& knots);

    void fill_DK
      (const std::vector<double>& knots);

    void fill_W
      (const std::vector<double> &weights);

    // Compute the S_l matrix for the penalization term
    void fill_S
      (const std::vector<double> & knots);

    void set_lambda
      (const std::vector<double> & knots);

public:

    Density(const std::vector<double>& knots, const std::vector<double>& cp, double kk, double g):
      k(kk), n(cp.size()), G(g+k+1), u(knots[0]), v(*knots.end())
    {
std::cout << "fill_C:" << '\n';
      weights.assign(n,1.0);
      set_lambda(knots);
      fill_C(cp, knots);
std::cout << "fill_M:" << '\n';
      fill_M(knots);
std::cout << "fill_DK:" << '\n';
      fill_DK(knots);
std::cout << "fill_W:" << '\n';
      fill_W(weights);
      P = (1 / alpha * (DK).transpose() * M * (DK) + (C * DK).transpose() * W * C * DK).sparseView();
      Eigen::VectorXd newcp(cp.size());
      for (int i = 0; i < cp.size() ; ++i) {
          newcp[i] = cp[i];
      }
      p = DK.transpose()* C.transpose() * W * newcp;
std::cout << "Constructor done:" << '\n' << p << '\n';
    }


    Density(const std::vector<double>& knots, const std::vector<double>& cp, double kk, double g, double opt_param):
      alpha(opt_param)
    {
      Density(knots, cp, kk,g);
    }

    Density(const std::vector<double>& knots, const std::vector<double>& cp, double kk, double g, unsigned int ll)
    {
      assert(ll<G);
      l = ll;
      Density(knots, cp, kk, g);
    }

    Density(const std::vector<double>& knots, const std::vector<double>& cp, double kk, double g, double opt_param, unsigned int ll):
      alpha(opt_param)
    {
      assert(ll<G);
      l = ll;
      Density(knots, cp, kk, g);
    }

    void print_all() const
    {
        std::cout << "MATRIX C:" << '\n' << C << '\n';
        std::cout << "MATRIX M:" << '\n' << M << '\n';
        std::cout << "MATRIX DK:" << '\n' << DK << '\n';
        std::cout << "MATRIX W:" << '\n' << W << '\n';
    }

    void solve()
    {   /* NAIVE SOLVER */
      Eigen::SparseLU<Eigen::SparseMatrix<double>>   solver;
      // Compute the ordering permutation vector from the structural pattern of A
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
      std::cout << "B-SPLINE COEFFICIENTS b = CKc" << '\n' << b << '\n';
    };
};

#endif //PROGETTO_DENSITY_ESTIMATION_HPP
