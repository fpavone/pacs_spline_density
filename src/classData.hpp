#ifndef STORE_DATA_1909_HPP
#define STORE_DATA_1909_HPP

#include <math.h>
#include <Eigen/Dense>
#include <iostream>
#include <iterator>
#include <vector>
#include <algorithm>
#include <functional>
#include "density_estimation.hpp"
#include "zeros.hpp"
//#include "find_type.hpp"

// NOTE: suppose data are stored as vector of vectors named data
// NOTE: data are supposed not to be ordered (otherwise counting occurrences is easier)
// NOTE: data are assumed to be already without extreme observations

// NOTE: how to divide data in subinterval-counts is a good challenge
// NOTE: if intervals are equispaced, it's easier;
// NOTE: geometric mean is computed in a straigth way: values are not so big,
// there's no reason to have overflow (i hope no underflow too, not so many classes)

//NOTE: rewrite coda input parameters
//NOTE: coda::BM modify by reference numbers
//NOTE: better to keep in memory nclasses instead of computing numbers.size in for loop

//NOTE: in what was "AUTO TYPE" with & or not??
//      for getData & is ok
//      for pacs & gives compile-time problems..


class myData
{
private:

  std::vector<double> numbers; //where data are stored (one row at a time)
  std::vector<double> grid;
  // unsigned int nclasses;

public:

  void
  readData
  (const Eigen::Block<Eigen::Map<Eigen::Matrix<double, -1, -1>, 0, Eigen::Stride<0, 0> >, 1, -1, false> & row,
         PRIOR prior, const int & cancel = -1);

  void
  transfData
  ();

  std::vector<double>
  getNumbers
  ();

  void
  pacs
  (myDensity & dens, Eigen::Block<Eigen::Matrix<double, -1, -1>, 1, -1, false> bspline);
  // bspline is the row of output matrix

  void
  antitData
  (Eigen::Block<Eigen::Matrix<double, -1, -1>, 1, -1, false> x);

  void
  fillGrid
  (double start, double end, unsigned int numPoints);

  void
  plotData_parallel
  (const myDensity & dens, unsigned long int numPoints,
    Eigen::Block<Eigen::Matrix<double, -1, -1>, 1, -1, false> bspline,
    Eigen::Block<Eigen::Matrix<double, -1, -1>, 1, -1, false> yplot);

  void
  plotData_parallel_Clr
  (const myDensity & dens, unsigned long int numPoints,
    Eigen::Block<Eigen::Matrix<double, -1, -1>, 1, -1, false> bspline,
    Eigen::Block<Eigen::Matrix<double, -1, -1>, 1, -1, false> yplot);

  long double
  compute_fvalue
  (Eigen::Block<Eigen::Matrix<double, -1, -1>, 1, -1, false> vec1, Eigen::ArrayXd vec2);

};

#endif //STORE_DATA_1909_HPP
