// #include "create_matrix.hpp"
#include "density_estimation.hpp"
#include <vector>
#include <iostream>

int main() {

    //TOY EXAMPLE
    std::vector<double> const knots{-1.0, -0.5, 0.5, 1.0}; //Knot sequence
    std::vector<double> const cp{.0, 0.8844, 0.8976, 0.0};
    std::vector<double> const weights{1.0, 1.0, 1.0, 1.0, 1.0};

    const int n = cp.size()-1; //Number of control points - 1
    const int k = 2; //Spline degree
    const int g = knots.size()-2;  //number of knots excluding the end points
    const int G = g + k + 1;

    Density<k, n, G> MyDensity(1.0, knots, cp, weights);
    MyDensity.print_all();
    return 0;
}
