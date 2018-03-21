#ifndef __COUNT_ZEROS_TREATMENT__
#define  __COUNT_ZEROS_TREATMENT__

#include <vector>
//#include <numeric>

// namespace help{
//   double sum(const std::vector<double> & vect)
//   {
//     return
//     std::accumulate<std::vector<double>::iterator, double>(vect.begin(),vect.end(),0.0);
//   };
//
//   std::vector<double> divide(const std::vector<double> &vect, const double & k)
//   {
//     std::vector<double> result(vect.size());
//     for(auto it = vect.begin(); it != vect.end(); it++)
//     {
//       result.push_back(*it/k);
//     }
//     return result;
//   };
// }


namespace coda{

  std::vector<double> uniform(unsigned int n){
    std::vector<double> tmp;
    tmp.insert(tmp.begin(),n,1.0/(double)(n))
    return tmp;
  };

  class BayesMultiplicative{

  private:
    const double & s; // Strength of prior information
    std::vector<double> & prior; // Dirichlet prior estimate

    const std::vector<double> & x; // Normalized input vector
    std::vector<double> & r; // Output vector

  public:
    BayesMultiplicative(const std::vector<double> & in, std::vector<double> & out)
       const double & strength = 1.0, const std::vector<double> & user_prior = uniform(in.size())):
      x(in), s(strength), r(out), prior(user_prior) {};

    void treat();
  };
}



#endif __COUNT_ZEROS_TREATMENT__
