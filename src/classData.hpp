#ifndef STORE_DATA_1909_HPP
#define STORE_DATA_1909_HPP

#include <math.h>
// #include <Eigen/Dense>
// #include <Eigen/Sparse>
#include <iostream>

// NOTE: suppose data are stored as vector of vectors named data
// NOTE: data are assumed to be already without extreme observations


class Mother {

  std::vector<std::vector<double>> data;
  std::vector<std::vector<double>> prop_data;
  std::vector<std::vector<double>> transf_data;

  unsigned int nclasses;

  unsigned int
  Sturge () const
  {
    double k=0.0;
    for(const auto& x:data)
      k += 1.0 + ceil(log2(x.size()));
    return (int) k/data.size();
  }

  bool loadData(const std::string& filepath)
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
    std::vector<double> v(nclasses,0.0);
    for(const auto& x:data){
      for (const auto& v:x)
        v[] ++;
      std::transform(v.begin(), v.end(), v.begin(),
               std::bind1st(std::divides<double>(),x.size()));
      prop_data.push_back(v);
    }



  }


};

#endif //STORE_DATA_1909_HPP
