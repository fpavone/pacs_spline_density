// #include "create_matrix.hpp"
#include "density_estimation.hpp"
#include <vector>
#include <iostream>

int main() {

    //1st ROW EXAMPLE
    std::vector<double> const knots{0.0, 30000, 70000, 110709}; //Knot sequence
    std::vector<double> const ycp{0.587, 2.331, 2.154, 1.271, 0.331, −0.550, −1.437, −1.997, −2.690};
    std::vector<double> const xcp{6574, 19591, 32608, 45625, 58641, 71658, 84675, 97692, 110709};
    const int n = xcp.size(); //corretto: cp.size()     OLD: cp.size()-1; //Number of control points - 1
    const int k = 3; //Spline degree
    unsigned int l = 2;
    const int g = knots.size()-2;  //number of knots excluding the end points
    const int G = g + k + 1;

    Density MyDensity(knots, xcp, ycp, k, g, l);
    MyDensity.print_all();
    MyDensity.solve();
    MyDensity.print_sol();
    return 0;
}
