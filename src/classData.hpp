#ifndef STORE_DATA_1909_HPP
#define STORE_DATA_1909_HPP

#include <math.h>
// #include <Eigen/Dense>
// #include <Eigen/Sparse>
#include <Eigen/Core>
#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>
#include <string>
#include "ctzeros.hpp"
#include "density_estimation.hpp"

// NOTE: suppose data are stored as vector of vectors named data
// NOTE: data are supposed not to be ordered (otherwise counting occurrences is easier)
// NOTE: data are assumed to be already without extreme observations

// NOTE: how to divide data in subinterval-counts is a good challenge
// NOTE: if intervals are equispaced, it's easier;
// NOTE: geometric mean is computed in a straigth way: values are not so big,
// there's no reason to have overflow (i hope no underflow too, not so many classes)



class Mother {

  friend class Density;

private:

  unsigned int k;  // Spline degree
  unsigned int l;  // order of derivative in penalization term
  double alpha;  // penalization parameter

  unsigned int n;   // Number of control points
  unsigned int g;   // Number knots - 2
  unsigned int G;   // Number of knots including additional ones

  std::vector<double> knots; // spline knots
  double u, v;    // [u,v] support of spline

  std::vector<double> xcp;

  // std::vector<std::vector<double>> data;
  std::vector<std::vector<double>> prop_data; // count data with zero replacement
  std::vector<std::vector<double>> transf_data; // clr-transformed data

  std::vector<Eigen::VectorXd> bspline;

  // NOTE: generalization for different data input
  unsigned int nclasses;
  double min, max;
  std::vector<double> intervals;

public:
  Mother(unsigned int kk, unsigned int ll, double opt_param):
    k(kk), l(ll), alpha(opt_param) {};

  void createKnots()
  {
    std::cout << "ATTENTO NODI NON CREATI: DA FINIRE.." << std::endl;
    g = knots.size()-2;
    G = g+k+1;
    u = knots.front();
    v = knots.back();
  }

  void readData(const std::string & fileD)
  {
    std::string line;

    std::ifstream file_stream(fileD, std::ios_base::in);
    if (!file_stream) {
      std::cout << "Could not open file " << fileD << std::endl;
    }

    // Read xcp
    getline(file_stream, line, '\n');
    std::stringstream line_stream(line);
    std::istream_iterator<double> start(line_stream), end;
    std::vector<double> numbers(start, end);
    for (const auto x:numbers)
      xcp.push_back(x);

    // Read prop_data
    while (getline(file_stream, line, '\n')) // Split up each line of the file.
    {
      std::stringstream line_stream(line);
      // Get a vector of numbers directly from a stream iterator.
      std::istream_iterator<double> start(line_stream), end;
      std::vector<double> numbers(start, end);
      // Taking account zero values replacement with BM
      prop_data.push_back(coda::BM(numbers));
    }

    n = xcp.size();
  }

  void readData(const std::string & fileD, const std::string & fileK)
  {
    readData(fileD);

    // Read knots
    std::ifstream file_stream(fileK, std::ios_base::in);
    if (!file_stream)
    {
      std::cout << "Could not open file " << fileK << std::endl;
    }

    getline(file_stream, line, '\n');
    std::stringstream line_stream(line);
    std::istream_iterator<double> start(line_stream), end;
    std::vector<double> numbers(start, end);
    for (const auto x:numbers)
      knots.push_back(x);

    g = knots.size()-2;
    G = g+k+1;
    u = knots.front();
    v = knots.back();
  }

  void
  transfData ()
  {
    // clr transformation of prop_data in transf_data
    double a = 1.0;
    std::vector<double> temp;

    for (auto& x:prop_data)
    {
      // computing geometric mean
      for (const auto& y:x)
        a *= y;
      a = pow(a, 1.0/nclasses);   // nclasses = x.size()

      // clr transformation
      for (const auto& y:x)
        temp.push_back(log(y/a));

      transf_data.push_back(temp);
      temp.clear();
    }
  }

  void
  pacs()
  {
    Density dens;
    for(const auto & it:transf_data)
    {
      dens.set_density(it);
      bspline.push_back(dens.solve());
    }
  }

  // NOTE: generalization for different data input

  unsigned int
  set_nclasses () const
  {
    // Adopting Sturge Formula
    double k=0.0;
    for(const auto& x:data)
      k += 1.0 + ceil(log2(x.size()));
    return (int) k/data.size();
  }

  void
  set_intervals()
  {
    double l = (max-min)/nclasses;
    unsigned int i = 0;
    while (i < nclasses+1)
      intervals.push_back(min + (i++)*l);
  }

  bool
  loadData (const std::string& filepath)
  {
    std::ifstream file_stream(filepath, std::ios_base::in);
    if (!file_stream) {
      std::cout << "Could not open file " << filepath << std::endl;
      return false;
    }

    std::string line;
    // Split up each line of the file.
    while (getline(file_stream, line, '\n')) {
      std::stringstream line_stream(line);

      // Get a vector of numbers directly from a stream iterator.
      std::istream_iterator<double> start(line_stream), end;
      std::vector<double> numbers(start, end);
      data.push_back(number);
    }
    return true;
  }

  void
  propData ()
  {
    // divide data in intervals, and count occurrences, and divide to
    // obtain proportions
    std::vector<double> v(nclasses,0.0);
    unsigned int it = 0;
    for(const auto& x:data)
    {
      for (const auto& y:x)
      {
        it = 0;
        // while (y > intervals[++it]) {};           v[--it] ++;
        while (y > intervals[it+1]) { it++; }
        v[it] ++;
      }
      std::transform(v.begin(), v.end(), v.begin(),
               std::bind1st(std::divides<double>(),x.size()));
      prop_data.push_back(v);
    }

    //        how to deal with zero counts
    // Because of count character of values in prop_data, also some zero values
    // occurred that would make further processing by clr transformation impossible.
    // For this reason, their imputation using a model-based procedure was performed
    prop_data = coda::BM(prop_data);
  }

  void
  antitData (std::vector<double>& v)
  {
    // anti clr transformation
    double a = 0.0;
    for (const auto& x:v)
      a += exp(x);
    for (auto& x:v)
      x = exp(x)/a;

  }

};

#endif //STORE_DATA_1909_HPP
