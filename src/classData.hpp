#ifndef STORE_DATA_1909_HPP
#define STORE_DATA_1909_HPP

#include <math.h>
#include <Eigen/Dense>
#include <iostream>
#include <iterator>
#include <vector>
#include <algorithm>
#include <functional>
//#include "zeros.hpp"
#include "density_estimation.hpp"
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


class myData {
private:
  std::vector<double> numbers; //where data are stored (one row at a time)
//  unsigned int nclasses;

  // Eigen::VectorXd bspline;

public:

  void
  getData
  (const Eigen::Block<Eigen::Map<Eigen::Matrix<double, -1, -1>,
                                  0, Eigen::Stride<0, 0> >, 1, -1, false> & row);

  void
  transfData
  ();

  std::vector<double>
  getNumbers   //temporanea 
  ()
  {
    return numbers;
  }

  void
  pacs
  (myDensity & dens, Eigen::Block<Eigen::Matrix<double, -1, -1>, 1, -1, false> bspline);
  // bspline is the row of output matrix

  void
  antitData (Eigen::Block<Eigen::Matrix<double, -1, -1>, 1, -1, false> x)
  {
    // anti clr transformation
    double a = exp(x.array()).sum();
    for(unsigned int i = 0; i < x.size(); i++)
     x[i] = exp(x[i])/a;
  }

  std::vector<double> grid;

  void fillGrid(double start, double end, unsigned int numPoints){
    double step = (end - start)/numPoints;
    grid.resize(numPoints);
    grid[0] = start;
    for (unsigned int i = 1; i < numPoints-1 ; ++i) {
      grid[i] = grid[i-1] + step;
    }
    grid[numPoints-1] = end;
  };

  void
  plotData_parallel
  (const myDensity & dens, unsigned long int numPoints,
    Eigen::Block<Eigen::Matrix<double, -1, -1>, 1, -1, false> bspline,
    Eigen::Block<Eigen::Matrix<double, -1, -1>, 1, -1, false> yplot)
  {

    double start = dens.get_u();
    double end = dens.get_v();
    unsigned int degree = dens.get_k();
    unsigned int G = dens.getG();
    // std::vector<double> knots;
    // set_lambda(degree, knots, pars.get_knots());
    const std::vector<double> knots = dens.get_lambda(); // NOTE: use a reference &?

    fillGrid(start, end, numPoints);

    Eigen::ArrayXd N;
    // yplot.resize(bspline.size());
  //  for (int row = 0; row < bspline.size(); ++row) {
      for (int i = 0; i < grid.size(); ++i) {
        int j = bspline::findspan(degree,grid[i],knots);
        N = Eigen::ArrayXd::Constant(G, 0.0);
        bspline::basisfun(j, grid[i], degree, knots, N);
        long double fvalue = compute_fvalue(bspline, N);
        yplot(i)=fvalue;
      }

      antitData(yplot);

      // for (int l = 0; l < yplot.size() ; ++l) {
      //   antitData(yplot[h]);
      // }
  //  }
  }

  void
  plotData_parallel_Clr
  (const myDensity & dens, unsigned long int numPoints,
    Eigen::Block<Eigen::Matrix<double, -1, -1>, 1, -1, false> bspline,
    Eigen::Block<Eigen::Matrix<double, -1, -1>, 1, -1, false> yplot)
  {

    double start = dens.get_u();
    double end = dens.get_v();
    unsigned int degree = dens.get_k();
    unsigned int G = dens.getG();
    // std::vector<double> knots;
    // set_lambda(degree, knots, pars.get_knots());
    const std::vector<double> knots = dens.get_lambda(); // NOTE: use a reference &?

    fillGrid(start, end, numPoints);

    Eigen::ArrayXd N;
    // yplot.resize(bspline.size());
  //  for (int row = 0; row < bspline.size(); ++row) {
      for (int i = 0; i < grid.size(); ++i) {
        int j = bspline::findspan(degree,grid[i],knots);
        N = Eigen::ArrayXd::Constant(G, 0.0);
        bspline::basisfun(j, grid[i], degree, knots, N);
        long double fvalue = compute_fvalue(bspline, N);
        yplot(i)=fvalue;
      }

    //  antitData(yplot);

      // for (int l = 0; l < yplot.size() ; ++l) {
      //   antitData(yplot[h]);
      // }
  //  }
  }

  long double
  compute_fvalue
  (Eigen::Block<Eigen::Matrix<double, -1, -1>, 1, -1, false> vec1, Eigen::ArrayXd vec2 )
  {
    long double res = 0.0;
    unsigned int n = vec1.size();
    if(vec2.size() != n){
      //gestisci errore
      exit(1);
    }
    for (int i = 0; i < n; ++i) {
      res += vec1[i]*vec2[i];
    }
    return res;
  }


};

#endif //STORE_DATA_1909_HPP
