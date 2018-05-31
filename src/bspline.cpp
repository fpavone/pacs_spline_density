#include "bspline.hpp"
#include <iostream>
#include <cassert>

using vect = std::vector<double>;


// Find the knot span of the parametric point u.
//
// INPUT:
//
//   n - number of control points - 1 //NOTE!!!!!!
//   p - spline degree
//   t - parametric point
//   U - knot sequence
//
// RETURN:
//
//   ret - knot span
//
// Note: This is NOT
// Algorithm A2.1 from 'The NURBS BOOK' pg68
// as that algorithm only works for nonperiodic
// knot vectors, nonetheless the results should
// be EXACTLY the same if U is nonperiodic


unsigned int
bspline::findspan
(int p, double t, const vect& U)
// ret is the index of the last knot at the left of the point t
{
	unsigned int n = U.size();
	unsigned int ret = 0;
	if (t > U[U.size () - 1] || t < U[0])
	{
		std::cerr << "Value " << t
	            << " of t is outside the knot span by "
	            << U[U.size () - 1] - t << "\n";
	    exit(EXIT_FAILURE);
	}
	else
	{
		while ((ret < n) && (U[ret] <= t))
		 	ret++;
	}
	if(ret>n-p-2) return n-p-2;
	return ret-1;
};	//findspan



void
bspline::basisfun
(unsigned int i, double t, int p, const vect& U, Eigen::ArrayXd& N)
/*
Evaluates basis splines which include x=t in the span and save values in N
*/
{

	double saved, temp;

	// work space
	double left[p+1];
	double right[p+1];
	if(i == p && t == U[i]){
		N[0] = 1.0;
	}
	else if(i == U.size() ){
		N[U.size()-p-2]= 1.0;
	}
	else{
		vect P(p + 1,1.0);
		for (unsigned int j = 1; j <= p; j++) {
			left[j] = t - U[i + 1 - j];
			right[j] = U[i + j] - t;
			saved = 0.0;
			for (unsigned int r = 0; r < j; r++) {
				temp = P[r] / (right[r + 1] + left[j - r]);
				P[r] = saved + right[r + 1] * temp;
				saved = left[j - r] * temp;
			}
			P[j] = saved;
		}
		for (unsigned int k = 0; k <= p; k++) {
			N[i - p + k] = P[k];
		}
	}

};	//basisfun
