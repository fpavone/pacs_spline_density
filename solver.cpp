#include "solver.hpp"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/QR>
// #include <Eigen/IterativeLinearSolvers>
using namespace Eigen;


void Solver::setup(){
  if(Method == MethodType::DEFAULT){
    if(A.rows() > NTreshold)
      Method = MethodType::ITERATIVE;
    else
      Method = MethodType::DIRECT;
  }
  if(Inv == IsInvertible::F && Method == MethodType::DIRECT) Option = SolverOption::QR;
  else Option = SolverOption::LSCG;
}

//gestire errori
VectorXd Solver::solve(const VectorXd & b){
switch(Option){
  case SolverOption::QR: {
    ColPivHouseholderQR<MatrixXd> dec(A);   // È buona cosa creare il  solver dentro una funzione?
                                            // terminata la funzione viene distrutto!
    Sol = dec.solve(b);
  }
  case SolverOption::LSCG: {
    // LeastSquaresConjugateGradient<MatrixXd> lscg;
    // lscg.compute(A);
    // Sol = lscg.solve(b);
  }
}
};
