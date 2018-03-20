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

    Eigen::Matrix<double, n, G> C;
    Eigen::Matrix<double, G, G> M;
    Eigen::Matrix<double, n, n> W;
    Eigen::SparseMatrix<double> S;

    Eigen::SparseMatrix<double> DK (G,G);
    Eigen::SparseMatrix<double> P (G,G); // matrix of the problem we have to solve GxG
    Eigen::VectorXd p; // known vector of the problem we have to solve G
    Eigen::VectorXd c; // c = P^(-)p G
    Eigen::VectorXd b; // B-spline coefficients G

    std::vector<double> weights;
    std::vector<double> lambda;  // extended vector of knots - with extra ones
                                // dimension: g + 2k + 2 = G + k + 1 ??

    void fill_C(const std::vector<double>& cp, const std::vector<double>& knots)
    {
        double t = 0.0;
        Eigen::ArrayXd N = Eigen::ArrayXd::Constant(G, 0.0);
        for (unsigned int i = 0; i < n; i++) {
            t = cp[i];
            int fs = bspline::findspan(n, k, t, knots);
            bspline::basisfun(fs, t, n, knots, N);
            for (unsigned int j = 0; j < G; j++) {
                C(i, j) = N(j);
            }
        }
    }

    void fill_M(const std::vector<double>& knots)
    {
        M.setZero();
        Eigen::ArrayXd N = Eigen::ArrayXd::Constant(G, 0.0);
        double x[n];
        double w[n];
        webbur::legendre_compute(n, x, w);
        for (int l = 0; l < n; ++l) {
            x[l] = (v - u) / 2 * x[l] + (v + u) / 2;
            w[l] = (v - u) / 2 * w[l];
        }
        int fs;
        double t;
        for (int i = 0; i < n; ++i) {
            N.setZero();
            t = x[i];
            fs = bspline::findspan(n, k, t, knots);
            bspline::basisfun(fs, t, n, knots, N);
            for (int j = 0; j < G; ++j) {
                for (int y = 0; y < G; ++y) {
                    M(j, y) += w[i] * N(j) * N(y);
                }
            }
        }
        //I need Gauss point weights
    }

    void fill_DK(const std::vector<double>& knots)
    {
        K.insert(0, 0) = (k+1)/(lambda[k+2] - knots[0]);
        K.insert(0, G-1) = -(k+1)/(knots[k+2] - knots[0]);
        for (std::size_t i = 1; i < G; i++) {
            K.insert(i-1,i) = -1/(lambda[k+2+i] - knots[i]);
            K.insert(i,i) = 1/(lambda[k+2+i] - knots[i]);
        }
    }

    void fill_W(const std::vector<double> &weights)
    {
        W.setZero();
        for (size_t i = 0; i < n; i++) {
            W(i, i) = weights[i];
        }
    }

    void fill_S(const std::vector<double> & knots)
    // Compute the S_l matrix for the penalization term
    {
      int l=2;
      for (size_t j = l; it >= 1; j--){
        Eigen::SparseMatrix<double> DL(G - it, G + 1 - it);
        /*
        There is a delay between the indexing of lambda in the reference paper
        [J. Machalova et al] and the indexing in the code. Precisely:
        index_code = index_ref + k
        This is why the following loop start from j instead of j - k.
        */
        for (size_t i = j; i < G - 1 ; i++) {
          DL.insert(i,i) = -(k + 1 - j)/(lambda[i+k+1-j] - lambda[i]);
          DL.insert(i,i+1) = (k + 1 - j)/(lambda[i+k+1-j] - lambda[i]);
        }
std::cout << "DL j = " << it << '\n' << Eigen::MatrixXd(S) << std::endl;
        if( it == l ) S = DL;
        else{
          S = S*DL;
          S.resize(G - l, G + 1 - it);
        }
      }
std::cout << "Printing S:" << '\n' << Eigen::MatrixXd(S) << std::endl;
    }

    void set_lambda(const std::vector<double> & knots)
    {
      lambda = std::vector(G + k + 1, knots[0]);
      lambda.insert(lambda.begin() + k + 1, knots.begin(), knots.end());
      lambda.insert(lambda.end() - k - 1 , knots.back());
std::cout << "lambda: " << '\n' << lambda << std::endl;
    }

public:

    Density(const std::vector<double>& knots, const std::vector<double>& cp, double kk):
      k(kk), n(cp.size()), G(g+k+1), u(knots[0]), v(*knots.end())
    {
      weigths.assign(n,1.0);
      fill_C(cp, knots);
      fill_M(knots);
      fill_DK(knots);
      fill_W(weights);
      P = (1 / alpha * (D * K).transpose() * M * (D * K) + (C * D * K).transpose() * W * C * D * K).sparseView();
      Eigen::VectorXd newcp(cp.size());
      for (int i = 0; i < cp.size() ; ++i) {
          newcp[i] = cp[i];
      }
      // prima era c =, penso che però sia p (Federico):
      p = K.transpose() * D.transpose() * C.transpose() * W * newcp;
      // as a consequence, i print p and not c (Gianluca)
      std::cout << "p:" << '\n' << p << '\n';
    }


    Density(const std::vector<double>& knots, const std::vector<double>& cp, double opt_param):
      alpha(opt_param)
    {
      Density(knots, cp);
    }

    Density(const std::vector<double>& knots, const std::vector<double>& cp, unsigned int ll):
    {
      assert(ll<G);
      l = ll;
      Density(knots, cp);
    }

    Density(const std::vector<double>& knots, const std::vector<double>& cp, double opt_param, unsigned int ll):
    {
      assert(ll<G);
      l = ll;
      Density(knots, cp, opt_param);
    }

    void print_all() const
    {
        std::cout << "MATRIX C:" << '\n' << C << '\n';
        std::cout << "MATRIX M:" << '\n' << C << '\n';
        std::cout << "MATRIX D:" << '\n' << D << '\n';
        std::cout << "MATRIX K:" << '\n' << K << '\n';
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
      b = D*K*c;
    };

    void print_sol() const
    {
      std::cout << "SOLUTION c = P^(-)p:" << '\n' << c << '\n';
      std::cout << "B-SPLINE COEFFICIENTS b = CKc" << '\n' << b << '\n';
    };
};

#endif //PROGETTO_DENSITY_ESTIMATION_HPP
