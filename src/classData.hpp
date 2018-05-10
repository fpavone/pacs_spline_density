#ifndef STORE_DATA_1909_HPP
#define STORE_DATA_1909_HPP

#include <math.h>
#include <Eigen/Dense>
// #include <Eigen/Sparse>
// #include <Eigen/Core>
#include <iostream>
#include <iterator>
#include <vector>
#include <algorithm>
#include <functional>
#include <string>
// #include "ctzeros.hpp"
#include "density_estimation.hpp"

// NOTE: suppose data are stored as vector of vectors named data
// NOTE: data are supposed not to be ordered (otherwise counting occurrences is easier)
// NOTE: data are assumed to be already without extreme observations

// NOTE: how to divide data in subinterval-counts is a good challenge
// NOTE: if intervals are equispaced, it's easier;
// NOTE: geometric mean is computed in a straigth way: values are not so big,
// there's no reason to have overflow (i hope no underflow too, not so many classes)

//NOTE: rewrite coda input parameters
//NOTE: coda::BM modify by reference numbers
//NOTE: better to keep in memory nclasses instead of computing numbers.size in for loop

std::vector<double>
divide(const std::vector<double> & vect, const double & D)
{
  assert(D>0 && " Error: dividing by zero..");
  std::vector<double> out;
  for(auto it:vect) out.push_back(it/D);
  return out;
};

double
sum(const std::vector<double> & vect)
{
  return std::accumulate(vect.begin(),vect.end(),0.0);
};

std::vector<double>
uniform(const unsigned int & n){
  std::vector<double> tmp;
  tmp.insert(tmp.begin(),n,1.0/(double)(n));
  return tmp;
};

double
geom_mean(const std::vector<double> & vect)
{
  double out = 1.0;
  for(const auto it:vect){
    out=out*it;
  }
  return pow(out, 1.0/vect.size());;
};

void
BM(std::vector<double> & numbers, const auto data)
{
  double s = 1.0;
  const bool  is_strength_inverse = false;
  std::vector<double> t = uniform(data.size());

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



class myData {
private:
  // std::vector<std::vector<double>> data;
  // std::vector<std::vector<double>> prop_data; // count data with zero replacement
  // std::vector<std::vector<double>> transf_data; // clr-transformed data

  std::vector<double> numbers; //where data are stored (one row at a time)
//  unsigned int nclasses;

  Eigen::VectorXd bspline;

public:

  void getData(const auto row)
  {
    numbers.clear();
    BM(numbers,row);
  }

  // void readData(const std::string & fileD)
  // {
  //   std::string line;
  //
  //   std::ifstream file_stream(fileD, std::ios_base::in);
  //   if (!file_stream) {
  //     std::cout << "Could not open file " << fileD << std::endl;
  //   }
  //
  //   // Read xcp and do nothing
  //   getline(file_stream, line, '\n');
  //
  //   // Read prop_data
  //   while (getline(file_stream, line, '\n')) // Split up each line of the file.
  //   {
  //     std::stringstream line_stream(line);
  //     // Get a vector of numbers directly from a stream iterator.
  //     std::istream_iterator<double> start(line_stream), end;
  //     std::vector<double> numbers(start, end);
  //     // Taking account zero values replacement with BM
  //     prop_data.push_back(coda::BM(numbers));
  //   }
  // }

  void
  transfData ()
  {
    // clr transformation of prop_data in transf_data
    double a = 1.0;

    // computing geometric mean
    for (const auto& y:numbers)
      a *= y;
    a = pow(a, 1.0/numbers.size());   // nclasses = x.size()

    // clr transformation
    for (auto& y:numbers)
      y = log(y/a);
  }

  void
  pacs(myDensity & dens, auto bspline)  // bspline is the row of output matrix
  {
    dens.set_density(numbers);
    dens.solve(bspline);
    // dens.print_sol();
  }

  // NOTE: generalization for different data input

  // unsigned int
  // set_nclasses () const
  // {
  //   // Adopting Sturge Formula
  //   double k=0.0;
  //   for(const auto& x:data)
  //     k += 1.0 + ceil(log2(x.size()));
  //   return (int) k/data.size();
  // }
  //
  // void
  // set_intervals()
  // {
  //   double l = (max-min)/nclasses;
  //   unsigned int i = 0;
  //   while (i < nclasses+1)
  //     intervals.push_back(min + (i++)*l);
  // }

  // bool
  // loadData (const std::string& filepath)
  // {
  //   std::ifstream file_stream(filepath, std::ios_base::in);
  //   if (!file_stream) {
  //     std::cout << "Could not open file " << filepath << std::endl;
  //     return false;
  //   }
  //
  //   std::string line;
  //   // Split up each line of the file.
  //   while (getline(file_stream, line, '\n')) {
  //     std::stringstream line_stream(line);
  //
  //     // Get a vector of numbers directly from a stream iterator.
  //     std::istream_iterator<double> start(line_stream), end;
  //     std::vector<double> numbers(start, end);
  //     data.push_back(number);
  //   }
  //   return true;
  // }

  // void
  // propData ()
  // {
  //   // divide data in intervals, and count occurrences, and divide to
  //   // obtain proportions
  //   std::vector<double> v(nclasses,0.0);
  //   unsigned int it = 0;
  //   for(const auto& x:data)
  //   {
  //     for (const auto& y:x)
  //     {
  //       it = 0;
  //       // while (y > intervals[++it]) {};           v[--it] ++;
  //       while (y > intervals[it+1]) { it++; }
  //       v[it] ++;
  //     }
  //     std::transform(v.begin(), v.end(), v.begin(),
  //              std::bind1st(std::divides<double>(),x.size()));
  //     prop_data.push_back(v);
  //   }
  //
  //   //        how to deal with zero counts
  //   // Because of count character of values in prop_data, also some zero values
  //   // occurred that would make further processing by clr transformation impossible.
  //   // For this reason, their imputation using a model-based procedure was performed
  //   prop_data = coda::BM(prop_data);
  // }

  // void
  // antitData (std::vector<double>& v)
  // {
  //   // anti clr transformation
  //   double a = 0.0;
  //   for (const auto& x:v)
  //     a += exp(x);
  //   for (auto& x:v)
  //     x = exp(x)/a;
  //
  // }

};

#endif //STORE_DATA_1909_HPP
