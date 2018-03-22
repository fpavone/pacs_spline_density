#include <vector>
#include <cassert>
#include <cmath>
#include <numeric>
#include "ctzeros.hpp"

std::vector<double>
help::divide(const std::vector<double> & vect, const double & D)
{
  assert(D>0 && " Error: dividing by zero..");
  std::vector<double> out;
  for(auto it:vect) out.push_back(it/D);
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



std::vector<double>
coda::BM(const std::vector<double> & in, const std::vector<double> & t, const double & s)
{
  assert(s>0 && " Error (BM): strength must be >0..");
  assert(in.size()==t.size() && " Error(BM): different sizes of input..");

  std::vector<double> out;
  double n = help::sum(in);
  double tmp = 0;

  double t_tot = 0;

  for(std::size_t i = 0; i < in.size(); i++){   //computing term of the summation
    assert(in[i]>=0 && " Error (BM): input must be >=0..");
    assert(t[i]>=0 && " Error (BM): input must be >=0..");

    if(in[i]==0) tmp += t[i];

    t_tot += t[i];
  }

  // assert()  controllare che t_tot sia "uguale" a 1

  for(std::size_t i = 0; i < in.size(); i++){  //applying BM method
    if(in[i] == 0)
        out.push_back(t[i]*s/(n + s));
    else
        out.push_back(in[i]*(1 - s*tmp/(n+s))/n);
  }
  return out;
};


std::vector<double>
coda::BM(const std::vector<double> & in, const double & s)
{
  return coda::BM(in, help::uniform(in.size()), s);
};


std::vector<double>
coda::BM(const std::vector<double> & in, coda::PRIOR p)
{
  if(p == coda::PRIOR::DEFAULT){
    if( help::sum(in) > (double) (in.size()*in.size()) )  // n > D^2
        p = coda::PRIOR::SQ;
    else
        p = coda::PRIOR::BAYES_LAPLACE;
  }

  double s;

  switch(p){
    case coda::PRIOR::PERKS :
    {
      s = 1;
      return coda::BM(in, s);
    }
    case coda::PRIOR::JEFFREYS :
    {
      s = (double)(in.size())/2;
      return coda::BM(in, s);
    }
    case coda::PRIOR::BAYES_LAPLACE :
    {
      s = (double)(in.size());
      return coda::BM(in, s);
    }
    case coda::PRIOR::SQ :
    {
      s = sqrt(help::sum(in));
      return coda::BM(in, s);
    }
  }
};
