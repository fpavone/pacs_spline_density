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
class Solver{
private:
  const MatrixXd & A;
  VectorXd Sol;
  MethodType Method;
  IsInvertible Inv;
  static const int NTreshold = 1000; // treshold tra metodo diretto e iterativo
  SolverOption Option;
public:
  Solver(const MatrixXd & AA, MethodType MethodInput = MethodType::DIRECT,
    IsInvertible InvInput = IsInvertible::F):
    A(AA), Sol(AA.rows()), Method(MethodInput), Inv(InvInput) {}; // CONTROLLARE #RIGHE = #COLONNE?
  // altri costruttori
  void setup(); // Per impostare il solver da utilizzare guardando il Method
  // e Inv (ed eventualmente il condizionamento)
  // Si potrebbe anche aggiungere caso se è simmetrica o definita/semidefinita
  // setup() dovrebbe funzionare come una factory, per adesso faccio classe enum dei solver
  void solve(const VectorXd & b);
  inline VectorXd getsolution(){return Sol;};
};

#endif
