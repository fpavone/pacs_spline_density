#ifndef __MY__SOLVER_
#define __MY__SOLVER_
#include <Eigen/Core>
// #include<Eigen/SparseLU>
// #include<Eigen/SparseQR>
// #include<Eigen/IterativeLinearSolvers>

using namespace Eigen;
// NOTE CONTROLLARE TUTTE LE & e I CONST!!!!

// Solving linear system of the type Ax = b
enum class MethodType {DIRECT, ITERATIVE};
enum class IsInvertible {T,F};
enum class SolverOption{QR,LU,LSCG};
// Se il metodo è automatico, guardare dimensione n e decidere se
// farlo iterativo o diretto

// NOTE AGGIUSTARE: COME PASSARE UN OGGETTO EIGEN PER REFERENCE
template <typename DerivedMat,typename DerivedVect>
class Solver
{
private:
  const EigenBase<DerivedMat>& A;
  VectorXd Sol;
  MethodType Method;
  IsInvertible Inv;
  static const int NTreshold = 1000; // treshold tra metodo diretto e iterativo
  SolverOption Option;
public:
  Solver(const EigenBase<DerivedMat>& AA, MethodType MethodInput = MethodType::DIRECT,
    IsInvertible InvInput = IsInvertible::F):
    A(AA), Method(MethodInput), Inv(InvInput) {}; // CONTROLLARE #RIGHE = #COLONNE?
  // altri costruttori
  void setup(); // Per impostare il solver da utilizzare guardando il Method
  // e Inv (ed eventualmente il condizionamento)
  // Si potrebbe anche aggiungere caso se è simmetrica o definita/semidefinita
  // setup() dovrebbe funzionare come una factory, per adesso faccio classe enum dei solver
  void solve(const EigenBase<DerivedVect>& b);
  inline VectorXd getsolution(){return Sol;};
};


template <typename DerivedMat,typename DerivedVect>
void Solver<DerivedMat,DerivedVect>::setup()
{
  if(Inv == IsInvertible::T && Method == MethodType::DIRECT) Option = SolverOption::LU;
  if(Inv == IsInvertible::F && Method == MethodType::DIRECT) Option = SolverOption::QR;
  else Option = SolverOption::LSCG;
};

template <typename DerivedMat,typename DerivedVect>
void Solver<DerivedMat,DerivedVect>::solve(const EigenBase<DerivedVect>& b)
{
  switch(Option)
  {
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

#endif
