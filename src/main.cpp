// #include "create_matrix.hpp"
// #include "density_estimation.hpp"
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <iterator>
#include "GetPot"

int main(int argc, char* argv[]) {

    // Read parameters
    GetPot commandline(argc, argv);
    const std::string fileName = commandline.follow("../input/parameters", 2, "-p", "--pars");
    GetPot filePara ( fileName.c_str() );

    unsigned int k = filePara("k", 3);     // Spline degree
    unsigned int l = filePara("l", 2);
    double alpha = filePara("alpha", 1.0);  // penalization parameter
    bool knots_given = filePara("knots_given", false);

    // Read data
    const std::string fileD = commandline.follow("../input/data", 1, "-d");
    std::ifstream file_stream(fileD, std::ios_base::in);
    if (!file_stream) {
      std::cout << "Could not open file " << fileD << std::endl;
    }

    std::string line;
    std::vector<double> xcp;
    std::vector<std::vector<double>> ycp;

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

    // Checking data are read correctly
    for (const auto x:xcp)
      std::cout<<x<<"  ";
    std::cout<<std::endl<<"ycp: "<<std::endl;
    for (const auto x:ycp){
      for (const auto y:x)
        std::cout<<y<<"  ";
      std::cout<<"\n";
    }



    // //1st ROW EXAMPLE
    // std::vector<double> const knots{0.0, 30000, 70000, 110709}; //Knot sequence
    //
    // unsigned int g = knots.size()-2;  //number of knots excluding the end points
    // // const int n = xcp.size(); //corretto: cp.size()     //Number of control points
    // // const int G = g + k + 1;
    //
    // Density MyDensity(knots, xcp, ycp, k, g, l);
    // MyDensity.print_all();
    // MyDensity.solve();
    // MyDensity.print_sol();
    return 0;
}
