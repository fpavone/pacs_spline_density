#ifndef __COUNT_ZEROS_TREATMENT__
#define __COUNT_ZEROS_TREATMENT__

#include <vector>

using dataframe = std::vector<std::vector<double>>;

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



namespace coda{
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


  std::vector<double>
  BM(const std::vector<double> & in, const std::vector<double> & t,
    const double & s = 1.0, const bool & is_strength_inverse = false);

  std::vector<double>
  BM(const std::vector<double> & in, const double & s);

  std::vector<double>
  BM(const std::vector<double> & in, coda::PRIOR p = coda::PRIOR::DEFAULT);

  dataframe
  BM(const dataframe & in, coda::PRIOR p = coda::PRIOR::DEFAULT);

  dataframe
  BM(const dataframe & in, const dataframe & t, const std::vector<double> & s, const bool & is_strength_inverse = false);

  /*
  GBM function:
    Prior estimate are computed through leave-one-out cross-validation.
    You need to pass a vector of vector of samples.
    Every sample is supposed to have the same dimension.

    Strength is computed as 1/g, where g is the geometric mean of the LOO estimates


    Reference: "Bayesian-multiplicative treatment of count zeros in compositional data sets"
    Authors: Josep-Antoni Martín-Fernández, Karel Hron, Matthias Templ,
              Peter Filzmoser and Javier Palarea-Albaladejo
    Periodical: Statistical Modelling 2015; 15(2): 134–158
  */

  dataframe
  GBM(const dataframe & in);
}



#endif
