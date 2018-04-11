#include "density_estimation.hpp"



void
Density::fill_C
(const std::vector<double>& cp)
{
    C.resize(n,G);
    Eigen::ArrayXd N;
    for (unsigned int i = 0; i < n; i++) {
      N = Eigen::ArrayXd::Constant(G, 0.0);
      unsigned int fs = bspline::findspan(k, cp[i], lambda);
      bspline::basisfun(fs, cp[i], k, lambda, N);
      C.row(i) = N;
    }
}

void
Density::fill_M
()
{
    M.resize(G-l,G-l);
    M.setZero();
    Eigen::ArrayXd N = Eigen::ArrayXd::Constant(G-l, 0.0);
    double x[n];
    double w[n];
    webbur::legendre_compute(n, x, w);
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
Density::fill_DK
()
{
  DK.resize(G,G);
//  DK.reserve(Eigen::VectorXi::Constant(G,2));
  DK.insert(0, 0) = (double)(k+1)/(lambda[k+1] - lambda[0]); // DK.insert(0, 0) = (k+1)/(lambda[k+2] - knots[0]);
  DK.insert(0, G-1) = -(double)(k+1)/(lambda[G + k] - lambda[G-1]); // DK.insert(0, G-1) = -(k+1)/(knots[k+1] - knots[0]);
  for (std::size_t i = 1; i < G; i++) {
      DK.insert(i,i-1) = -(double)(k+1)/(lambda[k+1+i] - lambda[i]);  // DK.insert(i,i-1) = -1/(lambda[k+2+i] - knots[i]);
      DK.insert(i,i) = (double)(k+1)/(lambda[k+1+i] - lambda[i]); // DK.insert(i,i) = 1/(lambda[k+2+i] - knots[i]);
  }
  DK.makeCompressed();
}


// Compute the S_l matrix for the penalization term
void
Density::fill_S
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
Density::set_lambda
(const std::vector<double> & knots)
{
  lambda.assign(k, knots[0]);
  lambda.insert(lambda.begin() + k, knots.begin(), knots.end());
  lambda.insert(lambda.end(), k ,knots.back());
}
void
Density::set_lambda_der
        (const std::vector<double> & knots)
{
    lambda_der.assign(k-l, knots[0]);
    lambda_der.insert(lambda_der.begin() + k-l, knots.begin(), knots.end());
    lambda_der.insert(lambda_der.end(), k-l ,knots.back());
}
