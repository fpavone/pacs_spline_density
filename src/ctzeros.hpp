#ifndef __COUNT_ZEROS_TREATMENT__
#define __COUNT_ZEROS_TREATMENT__

#include <vector>

namespace help{
  std::vector<double>
  divide(const std::vector<double> & vect, const double & D);

  double
  sum(const std::vector<double> & vect);

  std::vector<double>
  uniform(const unsigned int & n);
}


namespace coda{
  /* Input vector must be not scaled!
    - D: dimension of the vector
    - n: total mass of the vector
    - s: strength of the prior information
  */
  enum class PRIOR{PERKS,          // s = 1, t = 1/D
                   JEFFREYS,       // s = D/2, t = 1/D
                   BAYES_LAPLACE,  // s = D, t = 1/D
                   SQ,             // s = sqrt(n), t = 1/D
                   DEFAULT};       // if sqrt(n)> D --> SQ, otherwise BAYES_LAPLACE


  std::vector<double>
  BM(const std::vector<double> & in, const std::vector<double> & t, const double & s = 1.0);

  std::vector<double>
  BM(const std::vector<double> & in, const double & s = 1.0);

  std::vector<double>
  BM(const std::vector<double> & in, coda::PRIOR p = coda::PRIOR::DEFAULT);

  // GBM: prior information computed by cross-validation
  // void
  // GBM(const std::vector<double> & in, std::vector<double> & out, const double & s = 1.0)
  // {
  //
  // };

}



#endif
