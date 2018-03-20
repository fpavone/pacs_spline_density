// #include "create_matrix.hpp"
#include "density_estimation.hpp"
#include <vector>
#include <iostream>

int main() {

    //TOY EXAMPLE
    std::vector<double> const knots{-1.0, -0.5, 0.5, 1.0}; //Knot sequence
    std::vector<double> const cp{-0.5, -0.2, 0.2, 0.5};

    const int n = 4;//corretto: cp.size()     OLD: cp.size()-1; //Number of control points - 1
    const int k = 2; //Spline degree
    const int g = 2;//knots.size()-2;  //number of knots excluding the end points
    const int G = g + k + 1;

    Density MyDensity(knots, cp, k, g);
    MyDensity.print_all();
    MyDensity.solve();
    MyDensity.print_sol();
    return 0;
}
