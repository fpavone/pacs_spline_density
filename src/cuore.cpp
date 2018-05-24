#include "density_estimation.hpp"
#include "classData.hpp"
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <iterator>
#include <RcppEigen.h>
#include <Rcpp.h>
#include <thread>
#include <chrono>


using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]

extern "C"{
SEXP mymain(SEXP k_, SEXP l_, SEXP alpha_, SEXP data_, SEXP Xcp_, SEXP knots_, SEXP numPoints_)
{

  // Read parameters
  unsigned int k = INTEGER(k_)[0];     // Spline degree
  unsigned int l = INTEGER(l_)[0];
  double alpha = REAL(alpha_)[0];  // penalization parameter
  unsigned int numPoints = INTEGER(numPoints_)[0]; // number of points of the grid for the density

  myData obj;
  myDensity dens(myParameters(k,l,alpha));

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

  dens.set_matrix();
  // dens.print_all();

  Eigen::MatrixXd bsplineMat(nrow,dens.get_G());
  Eigen::MatrixXd yvalueMat(nrow,numPoints);
  Eigen::MatrixXd yvalueMatClr(nrow,numPoints);

  int barWidth = 70; // for bar progress plot

  for(std::size_t i = 0; i < nrow; i++)
  {
    obj.readData(data.row(i));
    obj.transfData();
    obj.pacs(dens, bsplineMat.row(i));
    obj.plotData_parallel(dens, numPoints, bsplineMat.row(i), yvalueMat.row(i));
    obj.plotData_parallel_Clr(dens, numPoints, bsplineMat.row(i), yvalueMatClr.row(i));

    std::cout << "[";
    int pos = barWidth * (double)i/nrow;
    for (int j = 0; j < barWidth; ++j) {
      if (j < pos) std::cout << "=";
      else if (j == pos) std::cout << ">";
      else std::cout << " ";
    }
    std::cout << "] " << int((double)i/(nrow-1) * 100.0) << "%\r";
    std::cout.flush();
    //std::this_thread::sleep_for(std::chrono::seconds(1));
  }
  std::cout << std::endl;

  List result = List::create(Named("bspline") = bsplineMat,
                             Named("Y") = yvalueMat,
                             Named("Y_clr") = yvalueMatClr,
                             Named("Numbers") = obj.getNumbers(),
                             Named("Xcp") = Xcp_,
                             Named("NumPoints") = numPoints_);

  return wrap(result);
};
}
