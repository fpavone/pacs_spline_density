#ifndef __MY__SOLVER_
#define __MY__SOLVER_

#include <Eigen/Core>
#include <Eigen/SparseLU>
#include <Eigen/SparseQR>
#include <Eigen/IterativeLinearSolvers>

using namespace Eigen;


enum class MethodType {DIRECT, ITERATIVE};
enum class IsInvertible {T,F};
enum class SolverOption{QR,LU,LSCG};



template <typename DerivedMat>
class Solver
{

private:
  const DerivedMat & A;
  VectorXd Sol;

  MethodType Method;
  IsInvertible Inv;
  SolverOption Option;

public:
  Solver(const DerivedMat & AA, MethodType MethodInput = MethodType::DIRECT,
    IsInvertible InvInput = IsInvertible::F):
    A(AA), Method(MethodInput), Inv(InvInput) {

      if(Inv == IsInvertible::T && Method == MethodType::DIRECT)
        Option = SolverOption::LU;
      if(Inv == IsInvertible::F && Method == MethodType::DIRECT)
        Option = SolverOption::QR;
      else Option = SolverOption::LSCG;

    }

  void solve(const VectorXd & b);

  inline VectorXd getsolution(){
    return Sol;
  }
};


//gestire errori
// NOTE È UNA BUONA IDEA CREARE I SOLVER DENTRO UNA FUNZIONE? ALLA FINE VENGONO DISTRUTTI

template <typename DerivedMat>
void Solver<DerivedMat>::solve(const DerivedVect & b)
{
  switch(Option)
  {
    case SolverOption::LU: {
      SparseLU< DerivedMat>   solver;
      // Compute the ordering permutation vector from the structural pattern of A
      solver.analyzePattern(A);
      // Compute the numerical factorization
      solver.factorize(A);
      //Use the factors to solve the linear system
      Sol = solver.solve(b);
    }
    case SolverOption::QR: {
      SparseLU<DerivedMat>   solver;
      // Compute the ordering permutation vector from the structural pattern of A
      solver.analyzePattern(A);
      // Compute the numerical factorization
      solver.factorize(A);
      //Use the factors to solve the linear system
      Sol = solver.solve(b);
    }
    case SolverOption::LSCG: {
      LeastSquaresConjugateGradient<DerivedMat> lscg;
      lscg.compute(A);
      Sol = lscg.solve(b);
    }
  }
};

#endif
