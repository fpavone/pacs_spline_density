#include <math.h>
//#include <Eigen/Dense>
#include <iostream>
#include <iterator>
#include <vector>
#include <algorithm>
#include <functional>
#include "zeros.hpp"

std::vector<double>
help::divide(const std::vector<double> & vect, const double & D)
{
  //assert(D>0 && " Error: dividing by zero..");
  std::vector<double> out;
  for(auto it:vect) out.push_back(it/D);
  return out;
};

double
help::sum(const std::vector<double> & vect)
{
  return std::accumulate(vect.begin(),vect.end(),0.0);
};

std::vector<double>
help::uniform(const unsigned int & n){
  std::vector<double> tmp;
  tmp.insert(tmp.begin(),n,1.0/(double)(n));
  return tmp;
};

double
help::geom_mean(const std::vector<double> & vect)
{
  double out = 1.0;
  for(const auto it:vect){
    out=out*it;
  }
  return pow(out, 1.0/vect.size());
};
