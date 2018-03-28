#include "density_estimation.hpp"



void
Density::fill_C
(const std::vector<double>& cp, const std::vector<double>& knots)
{
    C.resize(n,G);
    double t = 0.0;
    Eigen::ArrayXd N = Eigen::ArrayXd::Constant(G, 0.0);
    for (unsigned int i = 0; i < n; i++) {
      t = cp[i];
      int fs = bspline::findspan(k, t, knots);
      bspline::basisfun(fs, t, k, knots, N);
      C.row(i) = N;
    }
}

void
Density::fill_M
(const std::vector<double>& knots)
{
    M.resize(G,G);

    M.setZero();
    Eigen::ArrayXd N = Eigen::ArrayXd::Constant(G, 0.0);
    double x[n];
    double w[n];
    webbur::legendre_compute(n, x, w);
    for (unsigned int l = 0; l < n; ++l) {
        x[l] = (v - u) / 2 * x[l] + (v + u) / 2;
        w[l] = (v - u) / 2 * w[l];
    }
    int fs;
    double t;
    for (unsigned int i = 0; i < n; ++i) {
        N.setZero();
        t = x[i];
        fs = bspline::findspan(n, k, t, knots);
        bspline::basisfun(fs, t, n, knots, N);
        for (unsigned int j = 0; j < G; ++j) {
            for (int y = 0; y < G; ++y) {
                M(j, y) += w[i] * N(j) * N(y);
            }
        }
    }
    //I need Gauss point weights
}



void
Density::fill_DK
(const std::vector<double>& knots)
{
  DK.resize(G,G);
  DK.reserve(Eigen::VectorXi::Constant(G,2));
  DK.insert(0, 0) = (k+1)/(lambda[k+2] - knots[0]);
  DK.insert(0, G-1) = -(k+1)/(knots[k+2] - knots[0]);
  for (std::size_t i = 1; i < G; i++) {
      DK.insert(i,i-1) = -1/(lambda[k+2+i] - knots[i]);
      DK.insert(i,i) = 1/(lambda[k+2+i] - knots[i]);
  }
  DK.makeCompressed();
}


// Compute the S_l matrix for the penalization term
// void
// Density::fill_S
// (const std::vector<double> & knots)
// {
//   int l=2;
//   for (size_t j = l; it >= 1; j--){
//     Eigen::SparseMatrix<double> DL(G - it, G + 1 - it);
//     /*
//     There is a delay between the indexing of lambda in the reference paper
//     [J. Machalova et al] and the indexing in the code. Precisely:
//     index_code = index_ref + k
//     This is why the following loop start from j instead of j - k.
//     */
//     for (size_t i = j; i < G - 1 ; i++) {
//       DL.insert(i,i) = -(k + 1 - j)/(lambda[i+k+1-j] - lambda[i]);
//       DL.insert(i,i+1) = (k + 1 - j)/(lambda[i+k+1-j] - lambda[i]);
//     }
// std::cout << "DL j = " << it << '\n' << Eigen::MatrixXd(S) << std::endl;
//     if( it == l ) S = DL;
//     else{
//       S = S*DL;
//       S.resize(G - l, G + 1 - it);
//     }
//   }
// std::cout << "Printing S:" << '\n' << Eigen::MatrixXd(S) << std::endl;
// }



void
Density::set_lambda
(const std::vector<double> & knots)
{
  lambda.assign(k, knots[0]);
  lambda.insert(lambda.begin() + k, knots.begin(), knots.end());
  lambda.insert(lambda.end(), k ,knots.back());
}
