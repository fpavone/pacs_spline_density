// #include "create_matrix.hpp"
#include "density_estimation.hpp"
#include "ctzeros.hpp"
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <iterator>
#include "GetPot"

// hope
// hope bis

int main(int argc, char* argv[]) {

    // Read parameters
    GetPot commandline(argc, argv);
    const std::string fileName = commandline.follow("input/parameters", 2, "-p", "--pars");
    GetPot filePara ( fileName.c_str() );

    unsigned int k = filePara("k", 3);     // Spline degree
    unsigned int l = filePara("l", 2);
    double alpha = filePara("alpha", 1.0);  // penalization parameter
    bool knots_given = filePara("knots_given", 0);
    std::string line;

    std::vector<double> knots;
    std::vector<double> xcp;
    std::vector<std::vector<double>> ycp;

    // Read data
    const std::string fileD = commandline.follow("input/data", 1, "-d");
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

    // Read ycp
    while (getline(file_stream, line, '\n')) // Split up each line of the file.
    {
      std::stringstream line_stream(line);
      // Get a vector of numbers directly from a stream iterator.
      std::istream_iterator<double> start(line_stream), end;
      std::vector<double> numbers(start, end);
      ycp.push_back(numbers);
    }

    // Read knots if knots_given
    if(knots_given)
    {
      const std::string fileK = commandline.follow("input/knots", 2, "-k", "--knots");
      std::ifstream file_stream(fileK, std::ios_base::in);
      if (!file_stream) {
        std::cout << "Could not open file " << fileK << std::endl;
      }

      getline(file_stream, line, '\n');
      std::stringstream line_stream(line);
      std::istream_iterator<double> start(line_stream), end;
      std::vector<double> numbers(start, end);
      for (const auto x:numbers)
        knots.push_back(x);
    }
    else
    {

    }

    // Checking data are read correctly
    // for (const auto x:xcp)
    //   std::cout<<x<<"  ";
    // std::cout<<std::endl<<"ycp: "<<std::endl;
    // for (const auto x:ycp){
    //   for (const auto y:x)
    //     std::cout<<y<<"  ";
    //   std::cout<<"\n";
    // }

    // Testing first row
    std::vector<std::vector<double>> prop_data = coda::BM(ycp);
    std::vector<std::vector<double>> transf_data;
    std::vector<double> temp;
    double a = 1.0;
    for (auto& x:prop_data)
    {
      // computing geometric mean
      for (const auto& y:x)
        a *= y;
      a = pow(a, 1.0/x.size());   // nclasses = x.size()

      // clr transformation
      for (const auto& y:x)
        temp.push_back(log(y/a));

      transf_data.push_back(temp);
      temp.clear();
    }

    // std::cout<<"prop_data: "<<std::endl;
    // for (const auto x:prop_data){
    //   for (const auto y:x)
    //     std::cout<<y<<"  ";
    //   std::cout<<"\n";
    // }
    // std::cout <<"\n\n\n\n\n" << std::endl;
    //
    // std::cout<<"transf: "<<std::endl;
    // for (const auto x:transf_data){
    //   for (const auto y:x)
    //     std::cout<<y<<"  ";
    //   std::cout<<"\n";
    // }

    Density MyDensity(knots, xcp, transf_data[0], k, l, alpha);
    // MyDensity.print_all();
    // MyDensity.solve();
    // MyDensity.print_sol();
    return 0;
}
