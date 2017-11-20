#include "solver.hpp"
#include <Eigen/Dense>
#include <Eigen/LU>

void Solver::setup(){
  if(Method == MethodType::DEFAULT){
    if(A.rows() > NTreshold)
      Method = MethodType::ITERATIVE;
    else
      Method = MethodType::DIRECT;
  }
  if(Inv == F && MethodType::DIRECT) Solver = SolverOption::QR;
  else Solver = SolverOption::LSCG;
}

//gestire errori
VectorXd Solver::solve(const VectorXd & b){
switch(Solver){
  case QR: {
    ColPivHouseholderQR<MatrixXd> dec(A);
    Sol = dec.solve(b);
  }
  case LSCG: {
    LeastSquaresConjugateGradient<MatrixXd> lscg;
    lscg.compute(A);
    Sol = lscg.solve(b);
  }
}
};
