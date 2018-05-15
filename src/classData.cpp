#include <math.h>
#include <Eigen/Dense>
#include <iostream>
#include <iterator>
#include <vector>
#include <algorithm>
#include <functional>
#include "zeros.hpp"
#include "density_estimation.hpp"
#include "classData.hpp"


void
myData::transfData
()
{
  // clr transformation of prop_data in transf_data
  double a = 1.0;

  // computing geometric mean
  for (const auto& y:numbers)
    a *= y;
  a = pow(a, 1.0/numbers.size());   // nclasses = x.size()

  // clr transformation
  for (auto& y:numbers)
    y = log(y/a);
};

void
myData::getData
(const Eigen::Block<Eigen::Map<Eigen::Matrix<double, -1, -1>,
                                0, Eigen::Stride<0, 0> >, 1, -1, false> & row)
{
  numbers.clear();
  BM(numbers,row);

  std::cout << "numbers:\n " << "\n";

  for(const auto &x:numbers)
    std::cout << x << "\n";
};

void
myData::pacs
(myDensity & dens, Eigen::Block<Eigen::Matrix<double, -1, -1>, 1, -1, false> bspline)  // bspline is the row of output matrix
{
  dens.set_density(numbers);
  dens.solve(bspline);
  dens.print_sol();
};
