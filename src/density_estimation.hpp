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
//#include <unsupported/Eigen/SparseExtra> // to save the matrix
constexpr double tol = 1e-04;

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

  void createKnots(const unsigned int & size, const double & uu, const double & vv)
  {
    /* Creating equispaced knots between uu and vv */
    u = uu;
    v = vv;
    double step = (v - u)/(size-1);

    for(unsigned int i = 0; i < size; i++)
      knots.push_back(u + i*step);

    g = knots.size()-2;
    G = g+k+1;
  }

  void getKnots(const double * inputKnots, const unsigned int & size) // read knots by copy
  {
    for(std::size_t i=0; i < size; i++)
      knots.push_back(inputKnots[i]);

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

  void getXcp(const double * inputXcp, const unsigned int & size){
    for(std::size_t i=0; i < size; i++)
      xcp.push_back(inputXcp[i]);

    n = xcp.size();
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

  unsigned int getG() {return G;} //inline ??
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

    void set_matrix()
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
//std::cout << P.format(HeavyFmt) << '\n';
//std::cout << Eigen::MatrixXd(S) << '\n';
//std::cout << "P: " << '\n';
      // sarebbe meglio fare P.noalias()= .. per evitare copie inutili
      P.noalias() = 1.0 / alpha * (DK).transpose() * S.transpose() * M * S * (DK) +
            (C * DK).transpose() * weights.asDiagonal() * C * DK;
    }

    void set_density(const std::vector<double>& ycp)
    {
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

    void
    solve(auto bspline)
    {

      Eigen::FullPivHouseholderQR<Eigen::MatrixXd> solverQR;
      solverQR.compute(P);

      unsigned int dimKer = solverQR.dimensionOfKernel();
      double relative_error = 0.0;

      // // if(solverLDLT.info() == Eigen::Success){
      // //   c = solverLDLT.solve(p);
      // //   std::cout << "Solver LDLT.." << '\n' << c << '\n';
      // //   relative_error = (P*c - p).norm() / p.norm(); // norm() is L2 norm
      // //   std::cout << "\nThe relative error is:\n" << relative_error << std::endl;
      // // }
      // // else std::cout << "\nNot positive/negative semidefinite.. \n" << std::endl;
      std::cout << "\nSolver FullPivHouseholderQR.." << '\n' << c << '\n';

      switch(dimKer){
        case 0:
        {
          c = solverQR.solve(p);
          break;
        }
        case 1:
        /*
        we have to find a vector of the kernel and find minimal norm solution,
        in order to do this we exploit property of matrix Q in QR decomposition: last n-r column
        of Q are basis for Ker(P) where r=rank(P)
        */
        {
          c = solverQR.solve(p);
          Eigen::VectorXd kernel = solverQR.matrixQ().col(G-1);
          double scale = c.dot(kernel)/kernel.dot(kernel);
          c = c - scale*kernel;
          break;
        }
        case 2:
        {
          c = solverQR.solve(p);
          Eigen::VectorXd ker1= solverQR.matrixQ().col(G-1);
          Eigen::VectorXd ker2 = solverQR.matrixQ().col(G-2);
          double scale1 = 0.5*c.dot(ker1 - ker2*(ker1.dot(ker2)/ker2.dot(ker2)))/(ker1.dot(ker1)
                              - ker1.dot(ker2)*ker1.dot(ker2)/ker2.dot(ker2));
          double scale2 = -ker1.dot(ker2)/ker2.dot(ker2)*scale1 + 0.5*c.dot(ker2)/ker2.dot(ker2);
          c = c - scale1*ker1 - scale2*ker2;
          break;
        }
        default:   // case >1
        {
          std::cerr << "\n WARNING: kernel dimension of the problem exceeds 2," << '\n';
          std::cerr << " using Andrey Tychonoff regularization.." << std::endl;
          double tychlambda2 = 0.01;
          solverQR.compute(P.transpose()*P + tychlambda2*Eigen::MatrixXd::Identity(G,G));
          c = solverQR.solve(P.transpose()*p);
          break;
        }
      }

      relative_error = (P*c - p).norm() / p.norm(); // norm() is L2 norm
      std::cout << "\nThe relative error is:\n" << relative_error << std::endl;

      if(relative_error > tol)
      /*
        Least-square solution: we look for a solution in the col space projecting the b (in Ax=b)
        NOTE: QR should do this automatically without giving any warning (CHECK)
      */
        std::cerr << "\n WARNING: found least-square solution..." << std::endl;

      std::cout << "\n SOLUTION:\n" << c << '\n';

      // std::cout << '\n' << "MATRIX Q: \n" << '\n';
      // std::cout << solverQR.matrixQ() << '\n';
      // std::cout << '\n' << "KERNEL DIM of P: " << solverQR.dimensionOfKernel() << '\n';
      bspline.noalias() = DK*c;
      return;
    };

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

    // void save_matrix() const
    // {
    //   P.makeCompressed();
    //   Eigen::saveMarket(P, "density.mtx");
    //   Eigen::saveMarketVector(p, "density_b.mtx");
    // };
};

#endif //PROGETTO_DENSITY_ESTIMATION_HPP
