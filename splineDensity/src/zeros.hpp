#ifndef ZEROS_JACK_FRUSCIANTE_1903
#define ZEROS_JACK_FRUSCIANTE_1903

#include <math.h>
#include <Eigen/Dense>
#include <iostream>
#include <iterator>
#include <vector>
#include <algorithm>
#include <functional>

namespace help
{
  std::vector<double>
  divide
  (const std::vector<double> & vect, const double & D);

  double
  sum
  (const std::vector<double> & vect);

  std::vector<double>
  uniform
  (const unsigned int & n);

  double
  geom_mean
  (const std::vector<double> & vect);
}



/*
  BM function:
    Input vector must be not scaled!
    - D: dimension of the vector
    - n: total mass of the vector
    - s: strength of the prior information
    Additionally to input vector, you may want to pass:
    1) strength of the prior
    2) strength + prior estimates
    3) choose a method from PRIOR class (default is PRIOR::DEFAULT)
    example: BM(input, coda::PRIOR::BAYES::LAPLACE;
    Reference: "Bayesian-multiplicative treatment of count zeros in compositional data sets"
    Authors: Josep-Antoni Martín-Fernández, Karel Hron, Matthias Templ,
              Peter Filzmoser and Javier Palarea-Albaladejo
    Periodical: Statistical Modelling 2015; 15(2): 134–158
  */

enum class PRIOR{PERKS,          // s = 1, t = 1/D
                 JEFFREYS,       // s = D/2, t = 1/D
                 BAYES_LAPLACE,  // s = D, t = 1/D
                 SQ,             // s = sqrt(n), t = 1/D
                 DEFAULT};       // if sqrt(n)> D --> SQ, otherwise BAYES_LAPLACE

void
BM
(std::vector<double> & numbers,
 const Eigen::Block<Eigen::Map<Eigen::Matrix<double, -1, -1>, 0, Eigen::Stride<0, 0> >, 1, -1, false> & data,
 PRIOR p = PRIOR::DEFAULT);

void
BM
(std::vector<double> & numbers,
  const Eigen::Block<Eigen::Map<Eigen::Matrix<double, -1, -1>, 0, Eigen::Stride<0, 0> >, 1, -1, false> & data,
  const double & s, const bool is_strength_inverse = false);

#endif //ZEROS_JACK_FRUSCIANTE_1903
