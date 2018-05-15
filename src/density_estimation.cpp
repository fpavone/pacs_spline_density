#include <Eigen/Eigenvalues>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <iterator>
#include <vector>
#include "bspline.hpp"
#include "sandia_rules.hpp"
#include "gauss_points_weights.hpp"
#include "density_estimation.hpp"


/************* myParameters class ***************/

void
myParameters::createKnots
(const unsigned int & size, const double & uu, const double & vv)
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

void
myParameters::getKnots
(const double * inputKnots, const unsigned int & size) // read knots by copy
{
  for(std::size_t i=0; i < size; i++)
    knots.push_back(inputKnots[i]);

  g = knots.size()-2;
  G = g+k+1;
  u = knots.front();
  v = knots.back();


  std::cout << "knots:\n " << "\n";

  for(const auto &x:knots)
    std::cout << x << "\n";

}

void
myParameters::getXcp
(const double * inputXcp, const unsigned int & size)
{
  for(std::size_t i=0; i < size; i++)
    xcp.push_back(inputXcp[i]);

  n = xcp.size();

  std::cout << "xcp:\n " << "\n";

  for(const auto &x:xcp)
    std::cout << x << "\n";
}



/************* myDensity class ***************/

void
myDensity::fill_C
(const std::vector<double>& cp)
{
    C.resize(n,G);
    Eigen::ArrayXd N;
    for (unsigned int i = 0; i < n; i++) {
      N = Eigen::ArrayXd::Constant(G, 0.0);
      int fs = bspline::findspan(k, cp[i], lambda);
      bspline::basisfun(fs, cp[i], k, lambda, N);
      C.row(i) = N;
    }
}

void
myDensity::fill_M
()
{
    M.resize(G-l,G-l);
    M.setZero();
    Eigen::ArrayXd N = Eigen::ArrayXd::Constant(G-l, 0.0);
    std::vector<double> x(n);
    std::vector<double> w(n);
    integral::grule(x, w);
    for (unsigned int ll = 0; ll < n; ++ll) {
        x[ll] = (v - u) / 2 * x[ll] + (v + u) / 2;
        w[ll] = (v - u) / 2 * w[ll];
    }

    int fs;
    for (unsigned int i = 0; i < n; ++i) {
        N.setZero();
        fs = bspline::findspan(k-l, x[i], lambda_der);
        bspline::basisfun(fs, x[i], k-l, lambda_der, N);
        for (unsigned int j = 0; j < G-l; ++j) {
            for (unsigned int y = 0; y < G-l; ++y) {

                M(j, y) += w[i] * N(j) * N(y);
            }
        }
    }
    //I need Gauss point weights
}



void
myDensity::fill_DK
()
{
  DK.resize(G,G);
//  DK.reserve(Eigen::VectorXi::Constant(G,2));
  DK.insert(0, 0) = (double)(k+1)/(lambda[k+1] - lambda[0]); // DK.insert(0, 0) = (k+1)/(lambda[k+2] - knots[0]);
  DK.insert(0, G-1) = -(double)(k+1)/(lambda[k+1] - lambda[0]); //DK.insert(0, G-1) = -(double)(k+1)/(lambda[G + k] - lambda[G-1]);
  for (std::size_t i = 1; i < G; i++) {
      DK.insert(i,i-1) = -(double)(k+1)/(lambda[k+1+i] - lambda[i]);  // DK.insert(i,i-1) = -1/(lambda[k+2+i] - knots[i]);
      DK.insert(i,i) = (double)(k+1)/(lambda[k+1+i] - lambda[i]); // DK.insert(i,i) = 1/(lambda[k+2+i] - knots[i]);
  }
  DK.makeCompressed();
}


// Compute the S_l matrix for the penalization term
void
myDensity::fill_S
()
{
  for (std::size_t j = l; j >= 1; j--)
  {
    Eigen::SparseMatrix<double> DL(G - j, G + 1 - j);
    /*
    There is a delay between the indexing of lambda in the reference paper
    [J. Machalova et al] and the indexing in the code. Precisely:
    index_code = index_ref + k
    This is why the following loop start from j instead of j - k.
    */
    for (std::size_t i = j; i <= (G - 1) ; i++)
    {
      DL.insert(i-j,i-j) = -(double)(k + 1 - j)/(lambda[i+k+1-j] - lambda[i]);
      DL.insert(i-j,i-j+1) = (double)(k + 1 - j)/(lambda[i+k+1-j] - lambda[i]);
    }
    if( j == l )
    {
      S.resize(G - l, G + 1 - l);
      S = DL;
    }
    else
    {
      S = S*DL;
    }
  }
}


void
myDensity::set_lambda
(const std::vector<double> & knots)
{
  lambda.assign(k, knots[0]);
  lambda.insert(lambda.begin() + k, knots.begin(), knots.end());
  lambda.insert(lambda.end(), k ,knots.back());
}


void
myDensity::set_lambda_der
        (const std::vector<double> & knots)
{
    lambda_der.assign(k-l, knots[0]);
    lambda_der.insert(lambda_der.begin() + k-l, knots.begin(), knots.end());
    lambda_der.insert(lambda_der.end(), k-l ,knots.back());
}

void
myDensity::set_matrix
()
{
  weights = Eigen::VectorXd::Constant(n,1.0);
  set_lambda(knots);
  set_lambda_der(knots);
  fill_C(xcp);
  fill_M();
  fill_DK();
  fill_S();
  Eigen::IOFormat HeavyFmt(Eigen::FullPrecision, 0, ", ", ";\n", "[", "]", "[", "]");
  P.noalias() = 1.0 / alpha * (DK).transpose() * S.transpose() * M * S * (DK) +
        (C * DK).transpose() * weights.asDiagonal() * C * DK;
}

void
myDensity::set_density
(const std::vector<double>& ycp)
{
  Eigen::VectorXd newycp(ycp.size());
  for (unsigned int i = 0, nn = ycp.size(); i < nn ; ++i) {
      newycp[i] = ycp[i];
  }
  p.noalias() = DK.transpose()* C.transpose() * weights.asDiagonal() * newycp;
}

void
myDensity::solve
(Eigen::Block<Eigen::Matrix<double, -1, -1>, 1, -1, false> bspline)
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

void
myDensity::print_all
() const
{
    std::cout << "MATRIX C:" << '\n' << C << std::endl;
    std::cout << "MATRIX M:" << '\n' << M << std::endl;
    std::cout << "MATRIX DK:" << '\n' << Eigen::MatrixXd(DK) << std::endl;
    std::cout << "MATRIX W:" << '\n' << Eigen::MatrixXd(weights.asDiagonal()) << std::endl;
}
