#ifndef ZEROS_JACK_FRUSCIANTE_1903
#define ZEROS_JACK_FRUSCIANTE_1903

#include <math.h>
#include <Eigen/Dense>
#include <iostream>
#include <iterator>
#include <vector>
#include <algorithm>
#include <functional>
#include "find_type.hpp"

namespace help{
std::vector<double>
divide(const std::vector<double> & vect, const double & D);

double
sum(const std::vector<double> & vect);

std::vector<double>
uniform(const unsigned int & n);

double
geom_mean(const std::vector<double> & vect);
}

void
BM(std::vector<double> & numbers,
  const Eigen::Block<Eigen::Map<Eigen::Matrix<double, -1, -1>, 0, Eigen::Stride<0, 0> >, 1, -1, false> & data);

#endif //ZEROS_JACK_FRUSCIANTE_1903
