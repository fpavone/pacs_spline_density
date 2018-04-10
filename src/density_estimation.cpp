#include "density_estimation.hpp"
#include "classData.hpp"


void
Density::fill_C
(const std::vector<double>& cp)
{
    C.resize(obj.n,obj.G);
    Eigen::ArrayXd N;
    for (unsigned int i = 0; i < obj.n; i++) {
      N = Eigen::ArrayXd::Constant(obj.G, 0.0);
      int fs = bspline::findspan(obj.k, cp[i], lambda);
      bspline::basisfun(fs, cp[i], obj.k, lambda, N);
      C.row(i) = N;
    }
}

void
Density::fill_M
()
{

    M.resize(obj.G,obj.G);
    M.setZero();
    Eigen::ArrayXd N = Eigen::ArrayXd::Constant(G, 0.0);
    double x[obj.n];
    double w[obj.n];
    webbur::legendre_compute(obj.n, x, w);
    for (unsigned int l = 0; l < obj.n; ++l) {
        x[l] = (obj.v - obj.u) / 2 * x[l] + (obj.v + obj.u) / 2;
        w[l] = (obj.v - obj.u) / 2 * w[l];
    }

    int fs;
    for (unsigned int i = 0; i < obj.n; ++i) {
        N.setZero();
        fs = bspline::findspan(obj.k, x[i], lambda);
        bspline::basisfun(fs, x[i], obj.k, lambda, N);
        for (unsigned int j = 0; j < obj.G; ++j) {
            for (int y = 0; y < obj.G; ++y) {
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
  DK.resize(obj.G,obj.G);
//  DK.reserve(Eigen::VectorXi::Constant(G,2));
  DK.insert(0, 0) = (double)(obj.k+1)/(lambda[obj.k+1] - lambda[0]); // DK.insert(0, 0) = (k+1)/(lambda[k+2] - knots[0]);
  DK.insert(0, G-1) = -(double)(obj.k+1)/(lambda[obj.G + obj.k] - lambda[obj.G-1]); // DK.insert(0, G-1) = -(k+1)/(knots[k+1] - knots[0]);
  for (std::size_t i = 1; i < obj.G; i++) {
      DK.insert(i,i-1) = -(double)(obj.k+1)/(lambda[obj.k+1+i] - lambda[i]);  // DK.insert(i,i-1) = -1/(lambda[k+2+i] - knots[i]);
      DK.insert(i,i) = (double)(obj.k+1)/(lambda[obj.k+1+i] - lambda[i]); // DK.insert(i,i) = 1/(lambda[k+2+i] - knots[i]);
  }
  DK.makeCompressed();
}


// Compute the S_l matrix for the penalization term
void
Density::fill_S
()
{
  for (std::size_t j = obj.l; j >= 1; j--)
  {
    Eigen::SparseMatrix<double> DL(obj.G - j, obj.G + 1 - j);
    /*
    There is a delay between the indexing of lambda in the reference paper
    [J. Machalova et al] and the indexing in the code. Precisely:
    index_code = index_ref + k
    This is why the following loop start from j instead of j - k.
    */
    for (std::size_t i = j; i <= (obj.G - 1) ; i++)
    {
      DL.insert(i-j,i-j) = -(double)(obj.k + 1 - j)/(lambda[i+obj.k+1-j] - lambda[i]);
      DL.insert(i-j,i-j+1) = (double)(obj.k + 1 - j)/(lambda[i+obj.k+1-j] - lambda[i]);
    }
    if( j == obj.l )
    {
      S.resize(obj.G - obj.l, obj.G + 1 - obj.l);
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
  lambda.assign(obj.k, knots[0]);
  lambda.insert(lambda.begin() + obj.k, knots.begin(), knots.end());
  lambda.insert(lambda.end(), obj.k ,knots.back());
}
