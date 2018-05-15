#ifndef ZEROS_JACK_FRUSCIANTE_1903
#define ZEROS_JACK_FRUSCIANTE_1903

#include <math.h>
//#include <Eigen/Dense>
#include <iostream>
#include <iterator>
#include <vector>
#include <algorithm>
#include <functional>

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
BM(std::vector<double> & numbers, const auto data)
{
  double s = data.size();
  const bool  is_strength_inverse = false;
  std::vector<double> t = help::uniform(data.size());

  // assert(s>=0 && " Error (BM): strength must be >=0..");
  // assert(in.size()==t.size() && " Error(BM): different sizes of input..");

  double n = data.sum();
  double tmp = 0.0;

  double t_tot = 0.0;

  for(std::size_t i = 0; i < data.size(); i++)
  {   //computing term of the summation
    // assert(in[i]>=0 && " Error (BM): input must be >=0..");
    // assert(t[i]>=0 && " Error (BM): input must be >=0..");

    if(data(i)==0) tmp += t[i];

    t_tot += t[i];
  }

  // assert()  NOTE: need to check that t_tot is "equal" to 1
  for(std::size_t i = 0; i < data.size(); i++)
  {  //applying BM method
    if(data[i] == 0)
    {
      if(is_strength_inverse == false) numbers.push_back(t[i]*s/(n + s));
      else numbers.push_back(t[i]/(s*n + 1));
    }

    else
    {
      if(is_strength_inverse == false) numbers.push_back(data[i]*(1 - s*tmp/(n+s))/n);
      else numbers.push_back(data[i]*(1 - tmp/(s*n+1))/n);
    }

  }
  return;
}

#endif //ZEROS_JACK_FRUSCIANTE_1903
