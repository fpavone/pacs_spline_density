#ifndef __COUNT_ZEROS_TREATMENT__
#define  __COUNT_ZEROS_TREATMENT__

#include <vector>

namespace coda{

  std::vector<double> uniform(unsigned int n){
    std::vector<double> tmp;
    tmp.insert(tmp.begin(),n,1.0/(double)(n))
    return tmp;
  };

  class BM{

  private:
    const double & s; // Strength of prior information
    std::vector<double> & prior; // Dirichlet prior estimate

    const std::vector<double> & x; // Normalized input vector
    std::vector<double> & r; // Output vector

  public:
    // constructor: - riceve prior
    //              - usa una delle prior di default
    //              - fa scegliere al programma la prior --> SQ o Laplace
    BayesMultiplicative(const std::vector<double> & in, std::vector<double> & out)
       const double & strength = 1.0, const std::vector<double> & user_prior = uniform(in.size())):
      x(in), s(strength), r(out), prior(user_prior) {};

    void treat();
  };

  class GBM{

  private:
    const double & s; // Strength of prior information
    std::vector<double> & prior; // Dirichlet prior estimate

    const std::vector<double> & x; // Normalized input vector
    std::vector<double> & r; // Output vector

  public:
    // constructor: unica possibilità, la prior viene calcolata con crossvalidation
    BayesMultiplicative(const std::vector<double> & in, std::vector<double> & out)
       const double & strength = 1.0):
      x(in), s(strength), r(out) {};

    void set_prior();
    void treat();
  };
}



#endif __COUNT_ZEROS_TREATMENT__
