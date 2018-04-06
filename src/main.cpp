// #include "create_matrix.hpp"
#include "density_estimation.hpp"
#include "ctzeros.hpp"
#include "classData.hpp"
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <iterator>
#include "GetPot"

int main(int argc, char* argv[]) {

    // Read parameters
    GetPot commandline(argc, argv);
    const std::string fileName = commandline.follow("input/parameters", 2, "-p", "--pars");
    GetPot filePara ( fileName.c_str() );

    unsigned int k = filePara("k", 3);     // Spline degree
    unsigned int l = filePara("l", 2);
    double alpha = filePara("alpha", 1.0);  // penalization parameter
    bool knots_given = filePara("knots_given", 0);

    Mother obj(k,l,alpha);

    // Read data
    const std::string fileD = commandline.follow("input/data", 1, "-d");

    // Read knots if knots_given
    if(knots_given)
    {
      const std::string fileK = commandline.follow("input/knots", 2, "-k", "--knots");
      obj.readData(fileD,fileK);
    }
    else
    {
      obj.createKnots();
      obj.readData(fileD);
    }

    obj.transfData();
    obj.pacs();

    // Checking data are read correctly
    // for (const auto x:xcp)
    //   std::cout<<x<<"  ";
    // std::cout<<std::endl<<"ycp: "<<std::endl;
    // for (const auto x:ycp){
    //   for (const auto y:x)
    //     std::cout<<y<<"  ";
    //   std::cout<<"\n";
    // }

    return 0;
}
