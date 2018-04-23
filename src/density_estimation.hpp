#ifndef PROGETTO_DENSITY_ESTIMATION_HPP
#define PROGETTO_DENSITY_ESTIMATION_HPP

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>
#include "bspline.hpp"
#include "sandia_rules.hpp"
#include "gauss_points_weights.hpp"
//#include <unsupported/Eigen/SparseExtra> // to save the matrix

// NOTE: maybe useful a diagonal matrix W instead of weights.asDiagonal
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
  myParameters(const unsigned int kk, const unsigned int ll, const double opt_param):
    k(kk), l(ll), alpha(opt_param) {};

  void createKnots()
  {
    std::cout << "ATTENTO NODI NON CREATI: DA FINIRE.." << std::endl;
    g = knots.size()-2;
    G = g+k+1;
    u = knots.front();
    v = knots.back();
  }

  void readKnots(const std::string & fileK)
  {
    std::string line;

    // Read knots
    std::ifstream file_stream(fileK, std::ios_base::in);
    if (!file_stream)
    {
      std::cout << "Could not open file " << fileK << std::endl;
    }

    getline(file_stream, line, '\n');
    std::stringstream line_stream(line);
    std::istream_iterator<double> start(line_stream), end;
    std::vector<double> numbers(start, end);
    for (const auto x:numbers)
      knots.push_back(x);

    g = knots.size()-2;
    G = g+k+1;
    u = knots.front();
    v = knots.back();
  }

  void readXcp(const std::string & fileD)
  {
    std::string line;

    std::ifstream file_stream(fileD, std::ios_base::in);
    if (!file_stream) {
      std::cout << "Could not open file " << fileD << std::endl;
    }

    // Read xcp
    getline(file_stream, line, '\n');
    std::stringstream line_stream(line);
    std::istream_iterator<double> start(line_stream), end;
    std::vector<double> numbers(start, end);
    for (const auto x:numbers)
      xcp.push_back(x);

    n = xcp.size();
  }
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

    myDensity(const myParameters & input): myParameters(input) {};

    void set_density(const std::vector<double>& ycp)
    {
//std::cout << "fill_C.." << '\n';
      // weights.assign(n,1.0);
      weights = Eigen::VectorXd::Constant(n,1.0);
      set_lambda(knots);
      set_lambda_der(knots);
      fill_C(xcp);
//std::cout << C << std::endl;
//std::cout << "fill_M.." << '\n';
      fill_M();
//std::cout << M << std::endl;
//std::cout << "fill_DK.." << '\n';
      fill_DK();
//std::cout << Eigen::MatrixXd(DK) << '\n';
//std::cout << "fill_S.." << '\n';
      fill_S();
      Eigen::IOFormat HeavyFmt(Eigen::FullPrecision, 0, ", ", ";\n", "[", "]", "[", "]");
        std::cout << P.format(HeavyFmt) << '\n';
//std::cout << Eigen::MatrixXd(S) << '\n';
//std::cout << "P: " << '\n';
      // sarebbe meglio fare P.noalias()= .. per evitare copie inutili
      P.noalias() = 1.0 / alpha * (DK).transpose() * S.transpose() * M * S * (DK) +
            (C * DK).transpose() * weights.asDiagonal() * C * DK;
      Eigen::VectorXd newycp(ycp.size());
      for (unsigned int i = 0, nn = ycp.size(); i < nn ; ++i) {
          newycp[i] = ycp[i];
      }
      p.noalias() = DK.transpose()* C.transpose() * weights.asDiagonal() * newycp;
    }

    void print_all() const
    {
        std::cout << "MATRIX C:" << '\n' << C << std::endl;
        std::cout << "MATRIX M:" << '\n' << M << std::endl;
        std::cout << "MATRIX DK:" << '\n' << Eigen::MatrixXd(DK) << std::endl;
        std::cout << "MATRIX W:" << '\n' << Eigen::MatrixXd(weights.asDiagonal()) << std::endl;
    }

    Eigen::VectorXd
    solve()
    {   /* NAIVE SOLVER */
      //Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
      // Compute the ordering permutation vector from the structural pattern of P
      //solver.analyzePattern(P);
      // Compute the numerical factorization
      //solver.factorize(P);
      //Use the factors to solve the linear system
      //c = solver.solve(p);
      c = P.fullPivHouseholderQr().solve(p);
      b = DK*c;
      return b;
    };

    void print_sol() const
    {
      std::cout << "Matrix P: " << '\n' << P << '\n';
      std::cout << "SOLUTION c = P^(-)p:" << '\n' << c << '\n';
    //  std::cout << "B-SPLINE COEFFICIENTS b = DKc" << '\n' << b << '\n';
    //  std::cout << "PROVA P*c = p..?" << '\n' << Eigen::VectorXd(P*c) << '\n';
    };

    // void save_matrix() const
    // {
    //   P.makeCompressed();
    //   Eigen::saveMarket(P, "density.mtx");
    //   Eigen::saveMarketVector(p, "density_b.mtx");
    // };
};

#endif //PROGETTO_DENSITY_ESTIMATION_HPP
