#ifndef STORE_DATA_1909_HPP
#define STORE_DATA_1909_HPP

#include <math.h>
// #include <Eigen/Dense>
// #include <Eigen/Sparse>
#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>
#include "ctzeros.hpp"

// NOTE: suppose data are stored as vector of vectors named data
// NOTE: data are supposed not to be ordered (otherwise counting occurrences is easier)
// NOTE: data are assumed to be already without extreme observations

// NOTE: how to divide data in subinterval-counts is a good challenge
// NOTE: if intervals are equispaced, it's easier;
// NOTE: geometric mean is computed in a straigth way: values are not so big,
// there's no reason to have overflow (i hope no underflow too, not so many classes)



class Mother {

    // std::vector<std::vector<double>> data;
    // std::vector<std::vector<double>> prop_data;
    std::vector<std::vector<double>> transf_data;

    unsigned int k;  // Spline degree
    unsigned int l;
    double alpha;  // penalization parameter
    bool knots_given;

    unsigned int nclasses;
    double min, max;
    std::vector<double> intervals;

    Mother(unsigned int kk, unsigned int ll, double opt_param, bool given):
      k(kk), l(ll), alpha(opt_param), knots_given(given);

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
    transfData ()
    {
      // clr transformation of prop_data in transf_data
      double a = 1.0;
      for (auto& x:prop_data)
      {
        // computing geometric mean
        for (const auto& y:x)
          a *= y;
        a = pow(a, 1.0/nclasses);   // nclasses = x.size()

        // clr transformation
        for (const auto& y:x)
          transf_data.push_back(log(y/a));
      }
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
