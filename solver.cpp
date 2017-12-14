#include "solver.hpp"
#include <Eigen/Core>
#include<Eigen/SparseLU>
#include<Eigen/SparseQR>
#include <Eigen/IterativeLinearSolvers>
using namespace Eigen;

// NOTE QUANDO USARE QR O LU????

void Solver::setup(){
  if(Inv == IsInvertible::T && Method == MethodType::DIRECT) Option = SolverOption::LU;
  if(Inv == IsInvertible::F && Method == MethodType::DIRECT) Option = SolverOption::QR;
  else Option = SolverOption::LSCG;
}

//gestire errori
// NOTE È UNA BUONA IDEA CREARE I SOLVER DENTRO UNA FUNZIONE? ALLA FINE VENGONO DISTRUTTI

VectorXd Solver::solve(const VectorXd & b){
switch(Option){
  case SolverOption::LU: {
    SparseLU<SparseMatrix<scalar, ColMajor>, COLAMDOrdering<Index> >   solver;
    // Compute the ordering permutation vector from the structural pattern of A
    solver.analyzePattern(A);
    // Compute the numerical factorization
    solver.factorize(A);
    //Use the factors to solve the linear system
    Sol = solver.solve(b);
  }
  case SolverOption::QR: {
    SparseLU<SparseMatrix<scalar, ColMajor>, COLAMDOrdering<Index> >   solver;
    // Compute the ordering permutation vector from the structural pattern of A
    solver.analyzePattern(A);
    // Compute the numerical factorization
    solver.factorize(A);
    //Use the factors to solve the linear system
    Sol = solver.solve(b);
  }
  case SolverOption::LSCG: {
    LeastSquaresConjugateGradient<SparseMatrix<double> > lscg;
    lscg.compute(A);
    Sol = lscg.solve(b);
  }
}
};
