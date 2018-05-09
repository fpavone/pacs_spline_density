// #include "create_matrix.hpp"
#include "density_estimation.hpp"
#include "ctzeros.hpp"
#include "classData.hpp"
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <iterator>
#include <Eigen/Dense>
// #include "GetPot"
// #include <R.h>
// #include <Rinternals.h>
// #include <Rdefines.h>
#include <Rcpp.h>
//#include <RcppEigen.h>

using namespace Rcpp;

using Rmat as NumericMatrix;
using Rvec as NumericVector;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]

//NOTE: maybe it is better to deal with knots in R and pass only knots_
//NOTE: decide where to check if data_ is really a matrix
//NOTE: return a list

SEXP mymain(SEXP k_, SEXP l_, SEXP alpha_, SEXP knots_given_, Rmat data_, Rvec Xcp_, Rvec knots_)
{

    // Read parameters

    unsigned int k = INTEGER(k_)[0];     // Spline degree
    unsigned int l = INTEGER(l_)[0];
    double alpha = REAL(alpha_)[0];  // penalization parameter
    bool knots_given = INTEGER(knots_given_)[0];

    myData obj;
    myParameters pars(k,l,alpha);

    // Read data
    double *Xcp = REAL(Xcp_);
    unsigned int Xcpsize = LENGTH(Xcp_);
    pars.getXcp(Xcp,Xcpsize);

    // Read knots if knots_given
    if(knots_given)
    {
      double *knots = REAL(knots_);
      unsigned int knotsSize = LENGTH(knots_);
      pars.getKnots(knots,knotsSize);
      // const std::string fileK = commandline.follow("input/knots", 2, "-k", "--knots");
      //pars.readKnots(fileK);
    }
    else
    {
    //  pars.createKnots();
    }

    Eigen::Map<Eigen::MatrixXd> data(as<Eigen::Map<Eigen::MatrixXd>> (data_));

    unsigned int nrow = data.rows();

    myDensity dens(pars);
    dens.set_matrix();

    Eigen::MatrixXd bsplineMat(nrow,pars.getG());

    for(std::size_t i = 0; i < nrow; i++)
    {
      obj.getData(data.row(i));
      obj.transfData();
      obj.pacs(dens, bsplineMat.row(i));   //NOTE: check if eigen has problem with references
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
    List result = List::create(bsplineMat);

    return wrap(result);
};


// int main(int argc, char* argv[]) {
//
//     // Read parameters
//     GetPot commandline(argc, argv);
//     const std::string fileName = commandline.follow("input/parameters", 2, "-p", "--pars");
//     GetPot filePara ( fileName.c_str() );
//
//     unsigned int k = filePara("k", 3);     // Spline degree
//     unsigned int l = filePara("l", 2);
//     double alpha = filePara("alpha", 1.0);  // penalization parameter
//     bool knots_given = filePara("knots_given", 0);
//
//     myData obj;
//     myParameters pars(k,l,alpha);
//
//     // Read data
//     const std::string fileD = commandline.follow("input/data", 1, "-d");
//     obj.readData(fileD);
//     pars.readXcp(fileD);
//
//     // Read knots if knots_given
//     if(knots_given)
//     {
//       const std::string fileK = commandline.follow("input/knots", 2, "-k", "--knots");
//       pars.readKnots(fileK);
//     }
//     else
//     {
//     //  pars.createKnots();
//     }
//
//     obj.transfData();
//
//     obj.pacs(pars);
//
//     // Checking data are read correctly
//     // for (const auto x:xcp)
//     //   std::cout<<x<<"  ";
//     // std::cout<<std::endl<<"ycp: "<<std::endl;
//     // for (const auto x:ycp){
//     //   for (const auto y:x)
//     //     std::cout<<y<<"  ";
//     //   std::cout<<"\n";
//     // }
//
//     return 0;
// }
