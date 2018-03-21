#ifndef __COUNT_ZEROS_TREATMENT__
#define  __COUNT_ZEROS_TREATMENT__

#include <stdio>
#include <vector>
#include <numeric>
#include <Eigen>

namespace help{
  double total_sum(const std::vector<double> & vect)
  {
    return
    std::accumulate<std::vector<double>::iterator, double>(vect.begin(),vect.end(),0.0);
  };

  std::vector<double> divide(const std::vector<double> &vect, const double & k)
  {
    std::vector<double> result(vect.size());
    for(auto it = vect.begin(); it != vect.end(); it++)
    {
      result.push_back(*it/k);
    }
    return result;
  };
}


namespace coda{


  template < class T >// T is supposed to be either std::vector<double> or Eigen::Vector
  class BayesMultiplicative{

  private:
    double s; // Strength of prior information
    T t; // Dirichlet prior hyperparameters

    T x; // Normalized input vector
    T r; // Output vector

  public:
    BayesMultiplicative(const T & c, const bool & normalized):
      x( normalized ? c : help::divide(&c,help::total_sum(&c)) ){}
  };
}



#endif __COUNT_ZEROS_TREATMENT__
