#include <vector>
#include <iostream>
#include <cassert>
#include <cmath>
#include <numeric>
#include "ctzeros.hpp"

using dataframe = std::vector<std::vector<double>>;

/****** help namespace definitions *******/

std::vector<double>
help::divide(const std::vector<double> & vect, const double & D)
{
  assert(D>0 && " Error: dividing by zero..");
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
  for(auto it:vect){
    out=out*it;
  }
  return pow(out, 1.0/vect.size());;
}


/****** coda namespace definitions *******/

std::vector<double>
coda::BM(const std::vector<double> & in, const std::vector<double> & t,
  const double & s, const bool & is_strength_inverse)
{
  assert(s>=0 && " Error (BM): strength must be >=0..");
  assert(in.size()==t.size() && " Error(BM): different sizes of input..");

  std::vector<double> out;
  double n = help::sum(in);
  double tmp = 0;

  double t_tot = 0;

  for(std::size_t i = 0; i < in.size(); i++)
  {   //computing term of the summation
    assert(in[i]>=0 && " Error (BM): input must be >=0..");
    assert(t[i]>=0 && " Error (BM): input must be >=0..");

    if(in[i]==0) tmp += t[i];

    t_tot += t[i];
  }

  // assert()  NOTE: need to check that t_tot is "equal" to 1
  for(std::size_t i = 0; i < in.size(); i++)
  {  //applying BM method
    if(in[i] == 0)
    {
      if(is_strength_inverse == false) out.push_back(t[i]*s/(n + s));
      else out.push_back(t[i]/(s*n + 1));
    }

    else
    {
      if(is_strength_inverse == false) out.push_back(in[i]*(1 - s*tmp/(n+s))/n);
      else out.push_back(in[i]*(1 - tmp/(s*n+1))/n);
    }

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
  if(p == coda::PRIOR::DEFAULT)
  {
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

dataframe
coda::BM(const dataframe & in, coda::PRIOR p){
  const std::size_t N = in.size();
  const std::size_t D = in[0].size();

  dataframe out(N, std::vector<double>(D,0.0));

  for(std::size_t i=0; i < in.size(); i++)
  {
    out[i] = BM(in[i], p);
  }

  return out;
};

dataframe
coda::BM(const dataframe & in, const dataframe & t, const std::vector<double> & s, const bool & is_strength_inverse)
{
  const std::size_t N = in.size();
  const std::size_t D = in[0].size();

  dataframe out(N, std::vector<double>(D,0.0));

  for(std::size_t i = 0; i < N; i++)
  {
    out[i] = coda::BM(in[i],t[i],s[i],is_strength_inverse);
  }
  return out;
};

dataframe
coda::GBM(const dataframe & in)  // "in" is NxM
{
  const std::size_t N = in.size();
  const std::size_t D = in[0].size();

  dataframe alpha(N, std::vector<double>(D,0.0)); // 0 initialization
  dataframe t(N, std::vector<double>(D,0.0));

  std::vector<double> strength(N,0.0);

  for(std::size_t i = 0; i < N; i++)
  {
    for(std::size_t j = 0; j < D; j++)
    {
      for(std::size_t k = 0; k < N; k++){
        if(k != i) alpha[i][j] += in[k][j];
      }
    }
  }

  for(std::size_t i = 0; i < N; i++)
  {
    t[i] = help::divide(alpha[i],help::sum(alpha[i]));
    strength[i] = help::geom_mean(t[i]);
  }

  return coda::BM(in,t,strength, true);
};
