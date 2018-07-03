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

/*!
@brief	Data manager

@note	Suppose data are stored as vector of vectors named data and data are not ordered
@note Data are assumed to be already without extreme observations
*/

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

  std::vector<double> numbers; /** vector where data are stored (one row at a time)*/
  std::vector<double> grid;    /** vector of input knots*/
  // unsigned int nclasses;

public:
  /*!
	@brief	Read one row of data
  @details Read one row and apply Bayesian-multiplicative treatment of count zeros if necessary.
  @see BM()
	@param row Input row of data
	@param prior Prior for BM treatment
  @see PRIOR
	@param cancel
  */
  void
  readData
  (const Eigen::Block<Eigen::Map<Eigen::Matrix<double, -1, -1>, 0, Eigen::Stride<0, 0> >, 1, -1, false> & row,
         PRIOR prior, const int & cancel = -1);

  /*!
	@brief	Transform data using centered log-ratio (clr) function
  @details Apply clr transformation to data:
            \f$ y_i = \frac{z_i}{g(z_1,...,z_n)} \f$
            where g() is the geometric mean, \f$z_i\f$ are the elements in the row of data and \f$y_i\f$ are the transformed elements.
  */
  void
  transfData
  ();

  /*!
	@brief	Get processed data row
  */
  std::vector<double>
  getNumbers
  ();

  /*!
  @brief
  @details
  @param dens (Input) myDensity object
  @see myDensity
  @param bspline (Output) Row of output matrix with coefficients of the bspline
  */
  void
  pacs
  (myDensity & dens, Eigen::Block<Eigen::Matrix<double, -1, -1>, 1, -1, false> bspline);

  /*!
	@brief	Anti-transform data using the inverse of centered log-ratio (clr) function
  @details Apply \f$ clr^{-1} \f$ transformation:
            \f$ z_i = \frac{exp(y_i)}{n \sum\limits_{k=1}^{n}exp(y_k)} \f$
  @param x Vector to anti-transform
  */
  void
  antitData
  (Eigen::Block<Eigen::Matrix<double, -1, -1>, 1, -1, false> x);

  /*!
  @brief	Generate equispaced grid in [start,end]
  @param start Leftend point of the interval
  @param end Rightend point of the interval
  @param numPoints Number of points to generate in the interval
  */
  void
  fillGrid
  (double start, double end, unsigned int numPoints);

  /*!
  @brief	Generate points to plot
  @param numPoints (input) Number of points to plot
  @param bspline (input) Coefficients of bspline
  @param yplot (output) Anti-transformed values of the bspline evaluated in the the points generated for plot
  */
  void
  plotData_parallel
  (const myDensity & dens, unsigned long int numPoints,
    Eigen::Block<Eigen::Matrix<double, -1, -1>, 1, -1, false> bspline,
    Eigen::Block<Eigen::Matrix<double, -1, -1>, 1, -1, false> yplot);

  /*!
  @brief	Generate points to plot
  @param numPoints (input) Number of points to plot
  @param bspline (input) Coefficients of bspline
  @param yplot (output) clr transformed values of the bspline evaluated in the the points generated for plot
  */
  void
  plotData_parallel_Clr
  (const myDensity & dens, unsigned long int numPoints,
    Eigen::Block<Eigen::Matrix<double, -1, -1>, 1, -1, false> bspline,
    Eigen::Block<Eigen::Matrix<double, -1, -1>, 1, -1, false> yplot);

  /*!
  @brief Compute the value of the bspline in a point
  @param vec1 Coefficients of the bspline
  @param vec2 basis splines evaluated in the point
  @return Value of the bspline evaluated in the point point
  */
  long double
  compute_fvalue
  (Eigen::Block<Eigen::Matrix<double, -1, -1>, 1, -1, false> vec1, Eigen::ArrayXd vec2);

};

#endif //STORE_DATA_1909_HPP
