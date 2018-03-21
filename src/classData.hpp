#ifndef STORE_DATA_1909_HPP
#define STORE_DATA_1909_HPP

#include <math.h>
// #include <Eigen/Dense>
// #include <Eigen/Sparse>
#include <iostream>

// NOTE: suppose data are stored as vector of vectors named data

class Mother {

  std::vector<std::vector<double>> data;
  std::vector<std::vector<double>> transf_data;

  unsigned int
  Sturge () const
  {
    double k=0.0;
    for(const auto& x:data)
      k += 1.0 + ceil(log2(x.size()));
    return (int) k/data.size();
  }

  void
  transform_data ()
  {
    for(const auto& x:data)
      

  }


};

#endif //STORE_DATA_1909_HPP
