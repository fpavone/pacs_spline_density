#include "density_estimation.hpp"
#include "classData.hpp"
#include "zeros.hpp"
#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <sstream>
#include <iterator>
#include <RcppEigen.h>
#include <Rcpp.h>
#include <ctime>
#ifdef _OPENMP
  #include <omp.h>
#else
  #define omp_get_thread_num() 0
#endif
#include "timer.h"
#include <chrono>



using namespace Rcpp;


extern "C"{
SEXP smoothingSplines_(SEXP k_, SEXP l_, SEXP alpha_, SEXP data_, SEXP Xcp_, SEXP knots_, SEXP numPoints_, SEXP prior_, SEXP nCPU_, SEXP fast_)
{
  cns::timer<> t1,t4,t5;
  t1.start();

  #ifdef _OPENMP
    omp_set_dynamic(0);         // Explicitly disable dynamic teams
    omp_set_num_threads(INTEGER(nCPU_)[0]);
  #endif

  bool furious = INTEGER(fast_)[0];

  // Read parameters
  unsigned int k = INTEGER(k_)[0];     // Spline degree
  unsigned int l = INTEGER(l_)[0];
  double alpha = REAL(alpha_)[0];  // penalization parameter
  unsigned int numPoints = INTEGER(numPoints_)[0]; // number of points of the grid for the density
  unsigned int prior_num = INTEGER(prior_)[0];

  PRIOR prior = PRIOR::DEFAULT;

  switch(prior_num)
  {
    case 1: // "perks"
      prior = PRIOR::PERKS;
      break;
    case 2: // "jeffreys"
      prior = PRIOR::JEFFREYS;
      break;
    case 3: // "bayes_laplace"
      prior = PRIOR::BAYES_LAPLACE;
      break;
    case 4:  // "sq"
      prior = PRIOR::SQ;
      break;
    default: {};
  }

  dataManager obj;
  myDensity dens(parametersManager(k,l,alpha));

  // Read xcp
  double *Xcp = REAL(Xcp_);
  unsigned int Xcpsize = LENGTH(Xcp_);
  dens.readXcp(Xcp,Xcpsize);

  // Read knots
  double *knots = REAL(knots_);
  unsigned int knotsSize = LENGTH(knots_);
  dens.readKnots(knots,knotsSize);

  // Read data
  Eigen::Map<Eigen::MatrixXd> data(as<Eigen::Map<Eigen::MatrixXd>> (data_));
  unsigned int nrow = data.rows();
  furious = furious || (nrow < 100); // if not useful, progress bar will not be shown
  t4.start();
  dens.set_matrix();
  dens.set_system();

  Eigen::MatrixXd bsplineMat(nrow,dens.get_G());
  Eigen::MatrixXd yvalueMat(nrow,numPoints);
  Eigen::MatrixXd yvalueMatClr(nrow,numPoints);

  int barWidth = 70; // for progress bar plot
  int count = 0;
  int pos = 0;

  #pragma omp parallel private(obj) firstprivate(dens)// default(shared)
  {
  #pragma omp for
    for(std::size_t i = 0; i < nrow; i++)
    {
      obj.readData(data.row(i), prior);
      obj.transfData();
      obj.pacs(dens, bsplineMat.row(i));
      obj.plotData(dens, numPoints, bsplineMat.row(i), yvalueMat.row(i));
      obj.plotData_Clr(dens, numPoints, bsplineMat.row(i), yvalueMatClr.row(i));

      // Progress bar
      if(!furious)
      {
        #pragma omp critical
        {
          std::cout << "[";
          pos = barWidth * (double)count/nrow;
          for (int j = 0; j < barWidth; ++j)
          {
            if (j < pos) std::cout << "=";
            else if (j == pos) std::cout << ">";
            else std::cout << " ";
          }
          std::cout << "] " << int((double)count/(nrow-1) * 100.0) << "%\r";
          std::cout.flush();
          count++;
        }
      }
    }
  }
  t4.stop();
  t5.start();
  List result = List::create(Named("bspline") = bsplineMat,
                             Named("Y") = yvalueMat,
                             Named("Y_clr") = yvalueMatClr,
                             Named("Xcp") = Xcp_,
                             Named("NumPoints") = numPoints_);
  t5.stop();
  t1.stop();
  std::cout << "\nIt took "<< t1.elapsed<std::chrono::milliseconds>() <<" milliseconds. " << std::endl;
  std::cout << "\ncomputations "<< t4.elapsed<std::chrono::milliseconds>() <<" milliseconds. " << std::endl;
  std::cout << "\nwriting "<< t5.elapsed<std::chrono::milliseconds>() <<" milliseconds. "  << '\n\n'<< std::endl;

  return wrap(result);
};


/**************************************************************/
/**************************************************************/
/******************* VALIDATION FUNCTION **********************/
/**************************************************************/
/**************************************************************/



SEXP smoothingSplinesValidation_(SEXP k_, SEXP l_, SEXP alpha_, SEXP data_, SEXP Xcp_, SEXP knots_, SEXP prior_, SEXP nCPU_)
{
  cns::timer<> t,t2;
  t.start();

  #ifdef _OPENMP
    omp_set_dynamic(0);         // Explicitly disable dynamic teams
    omp_set_num_threads(INTEGER(nCPU_)[0]);
  #endif

  // Read parameters
  unsigned int k = INTEGER(k_)[0];     // Spline degree
  unsigned int l = INTEGER(l_)[0];
  double * alpha = REAL(alpha_);  // penalization parameters vector
  unsigned int alpha_size= LENGTH(alpha_);
  unsigned int prior_num = INTEGER(prior_)[0];

  PRIOR prior = PRIOR::DEFAULT;

  switch(prior_num)
  {
    case 1: // "perks"
      prior = PRIOR::PERKS;
      break;
    case 2: // "jeffreys"
      prior = PRIOR::JEFFREYS;
      break;
    case 3: // "bayes_laplace"
      prior = PRIOR::BAYES_LAPLACE;
      break;
    case 4:  // "sq"
      prior = PRIOR::SQ;
      break;
    default: {};
  }

  dataManager obj;
  myDensity dens(parametersManager(k,l,alpha[0]));

  // Read xcp
  double *Xcp = REAL(Xcp_);
  unsigned int Xcpsize = LENGTH(Xcp_);

  // Read knots
  double *knots = REAL(knots_);
  unsigned int knotsSize = LENGTH(knots_);
  dens.readKnots(knots,knotsSize);

  // Read data
  Eigen::Map<Eigen::MatrixXd> data(as<Eigen::Map<Eigen::MatrixXd>> (data_));

  unsigned int nrow = data.rows();
  unsigned int ncol = data.cols();

  int barWidth = 70; // for bar progress plot
  int count = 0;
  int pos = 0;
  unsigned int min_index = 0;
  Eigen::VectorXd Jvalues = Eigen::VectorXd::Zero(alpha_size);
  Eigen::VectorXd CVerror = Eigen::VectorXd::Zero(alpha_size,0.0);
  double CVopt = -1.0;
  t2.start();
  #pragma omp parallel private(obj) firstprivate(dens)// default(shared)
  {
    Eigen::MatrixXd threadBsplineMat(nrow,dens.get_G());
    Eigen::ArrayXd N;
    int span;
    long double fvalue;

    #pragma omp for
    for(std::size_t z = 0; z < alpha_size; z++)
    {
      dens.set_alpha(alpha[z]);
      for(std::size_t i = 0; i < nrow; i++)
      {
        for(std::size_t j = 1; j < ncol-1; j++)    // we do not leave out extremals, otherwise we have to change the knots
        {
          // LOO-CV, leave out column j, fit spline and compute error
          dens.readXcp(Xcp,Xcpsize,j);
          dens.set_matrix();
          dens.set_system();
          obj.readData(data.row(i),prior,j);
          obj.transfData();
          obj.pacs(dens, threadBsplineMat.row(i));
          Jvalues(z)+=dens.eval_J(obj.getNumbers())/nrow;

          // computing error using spline coefficients, fvalue is the predicted value
          span = bspline::findspan(k,Xcp[j],dens.get_lambda());
          N = Eigen::ArrayXd::Constant(dens.get_G(), 0.0);
          bspline::basisfun(span,Xcp[j],k,dens.get_lambda(), N);
          fvalue = obj.compute_fvalue(threadBsplineMat.row(i), N);
          CVerror(z) += fabs(fvalue - data(i,j))/nrow;
        }
      }
    }
  }
  for(std::size_t z = 0; z < nrow; z++)
  {
    if(CVerror(z)<CVopt || CVopt<0)
    {
      CVopt = CVerror(z);
      min_index = z;
    }
  }

  t2.stop();

  // Printing out J[alpha]
  for(std::size_t k = 0; k < alpha_size; k++)
  {
    Rcout << std::string((k==min_index)?" ***":" ") << "\talpha = "<< alpha[k] << "\tJ = " << Jvalues(k)
          << "\t CV-error: " << CVerror(k) << '\n';
  }



  List result = List::create(Named("alpha") = alpha_,
                             Named("Jvalues") = Jvalues,
                             Named("CVerror") = CVerror);


  t.stop();
  std::cout << "\nIt took "<< t.elapsed<std::chrono::milliseconds>() <<" milliseconds. " << std::endl;
  std::cout << "\ncomputations "<< t2.elapsed<std::chrono::milliseconds>() <<" milliseconds. " << std::endl;

  return wrap(result);
};
}
