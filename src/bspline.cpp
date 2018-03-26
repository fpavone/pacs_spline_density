#include "bspline.hpp"
#include <iostream>
#include <cassert>

using vect = std::vector<double>;

int
bspline::findspan (int p, double t, const vect& U)
{
	n = U.size();
	int ret = 0;
	if (t > U[U.size () - 1] || t < U[0])
	{
		std::cerr << "Value " << t
	            << " of t is outside the knot span by "
	            << U[U.size () - 1] - t << "\n";
	    exit(EXIT_FAILURE);
	}
	else
	{
		while ((ret++ < n) && (U[ret] <= t)) { };
	}
	return (ret-1);
};	//findspan

void
bspline::basisfun (int i, double t, int p, const vect& U, Eigen::ArrayXd& N)
{
	int j,r;
	double saved, temp;

	// work space
	double left[p+1];
	double right[p+1];

	vect  P(p+1);
	P[0] = 1.0;
	for (j = 1; j <= p; j++) {
		left[j]  = t - U[i+1-j];
		right[j] = U[i+j] - t;
		saved = 0.0;
		for (r = 0; r < j; r++) {
			temp = P[r] / (right[r+1] + left[j-r]);
			P[r] = saved + right[r+1] * temp;
			saved = left[j-r] * temp;
		}
		P[j] = saved;
	}
	for (int k = 0; k <= p ; ++k) {
		N[i-p+k]=P[k];
	}
};	//basisfun
