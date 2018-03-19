#ifndef PROGETTO_DENSITY_ESTIMATION_HPP
#define PROGETTO_DENSITY_ESTIMATION_HPP

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include "bspline.hpp"
#include "sandia_rules.hpp"

template<int k, int n, int G>
class Density {

private:
    double alpha = 0.0;
    double u;
    double v;

    Eigen::Matrix<double, n, G> C;
    Eigen::Matrix<double, G, G> M;
    Eigen::Matrix<double, G, G> D;
    Eigen::Matrix<double, G, G> K;
    Eigen::Matrix<double, n, n> W;
    Eigen::SparseMatrix<double> S;

    Eigen::SparseMatrix<double> P; // matrix of the problem we have to solve GxG
    Eigen::VectorXd p; // known vector of the problem we have to solve G
    Eigen::VectorXd c; // c = P^(-)p G
    Eigen::VectorXd b; // B-spline coefficients G


    void fill_C(std::vector<double> const &cp, std::vector<double> const &knots) {
// std::cout<<"Fill C: start.. "<<std::endl;
// std::cout<<"vecctor cp:"<<std::endl;
// for(auto x:cp) std::cout<<x<<'\t'; std::cout<<std::endl;
// std::cout<<"vecctor knots:"<<std::endl;
        double t = 0;
        Eigen::ArrayXd N = Eigen::ArrayXd::Constant(G, 0.0);
// for(auto x:knots) std::cout<<x<<'\t'; std::cout<<std::endl;
// std::cout<<"G: "<<G <<" N"<<'\n'<<N<<std::endl;
        for (unsigned int i = 0; i < n; i++) {
// std::cout<<"i.. "<<i<<std::endl;
            t = cp[i];
// std::cout<<"fs.. "<<std::endl;
            int fs = bspline::findspan(n, k, t, knots);
// std::cout<<"basisfun.. "<<std::endl;
            bspline::basisfun(fs, t, n, knots, N);
            for (unsigned int j = 0; j < G; j++) {
                C(i, j) = N(j);
            }
        }
// std::cout<<"Fill C: ..end "<<std::endl;
    }

    void fill_M(std::vector<double> const &knots) {
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

    void fill_D(std::vector<double> const &knots) {
        D.setZero();
        for (size_t i = 0; i < G; i++) {
            D(i, i) = (k + 1) / (knots[i] - knots[i - k + 1]);
        }
    }

    void fill_K() {
        K.setZero();
        K(0, G - 1) = -1;
        K(G - 1, G - 1) = 1;
        for (size_t i = 0; i < G - 1; i++) {
            K(i, i) = 1;
            K(i + 1, i) = -1;
        }
    }

    void fill_W(std::vector<double> const &weights) {
        W.setZero();
        for (size_t i = 0; i < n; i++) {
            W(i, i) = weights[i];
        }
    }

    void fill_S(std::vector<double> const & knots) {
      int l=2;
      for (size_t it = l; it >= 1; it--){
        Eigen::SparseMatrix<double> DL(G - it, G + 1 - it);
        for (size_t i = 0; i < G - it ; i++) {
          DL.insert(i,i) = -(k + 1 - it)/(knots[i+k+1-it] - knots[i]);
          DL.insert(i,i+1) = (k + 1 - it)/(knots[i+k+1-it] - knots[i]);
        }
        std::cout << "DL it = " << it << '\n' << Eigen::MatrixXd(S) << std::endl;
        if( it == l ) S = DL;
        else{
          S = S*DL;
          S.resize(G - l, G + 1 - it);
        }
      }
      std::cout << "Printing S:" << '\n' << Eigen::MatrixXd(S) << std::endl;
    }


public:
    Density(double opt_param, std::vector<double> const knots, std::vector<double> const cp,
            std::vector<double> const weights) {
//        n = cp.size() - 1;
//        g = knots.size() - 2;
//        G = g + k + 1;
        alpha = opt_param;
        u = knots[0];
        v = *knots.end();
        fill_C(cp, knots);
        fill_M(knots);
        fill_D(knots);
        fill_K();
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
        fill_S(knots);
    }

    void print_all() const {
        std::cout << "MATRIX C:" << '\n' << C << '\n';
        std::cout << "MATRIX M:" << '\n' << C << '\n';
        std::cout << "MATRIX D:" << '\n' << D << '\n';
        std::cout << "MATRIX K:" << '\n' << K << '\n';
        std::cout << "MATRIX W:" << '\n' << W << '\n';
    }

    void solve() {   /* NAIVE SOLVER */
      Eigen::SparseLU<Eigen::SparseMatrix<double>>   solver;
      // Compute the ordering permutation vector from the structural pattern of A
      solver.analyzePattern(P);
      // Compute the numerical factorization
      solver.factorize(P);
      //Use the factors to solve the linear system
      c = solver.solve(p);
      b = D*K*c;
    };

    void print_sol() const{
      std::cout << "SOLUTION c = P^(-)p:" << '\n' << c << '\n';
      std::cout << "B-SPLINE COEFFICIENTS b = CKc" << '\n' << b << '\n';
    };
};

#endif //PROGETTO_DENSITY_ESTIMATION_HPP
