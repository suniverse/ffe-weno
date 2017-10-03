#include <cmath>
#include "Interpolation.hpp"

int sgn(double x)
{
	return (x > 0) - (x < 0);
}

double minmod(double a, double b, double c)
{
	double s = fmin(std::fabs(a),std::fabs(b));
	s = fmin(s,std::fabs(c));
	s = s * std::fabs(sgn(a)+sgn(b))/4.*(sgn(a)+sgn(c));
	return s;
}



double PLM(double a, double b, double c, int i)
{
	// i<0 for left state, i>0 for right state

	const double theta = 1.5;
	double u1 = theta * (b-a);
	double u2 = (c-a)/2.;
	double u3 = theta * (c-b);
	double deltau = 0.5*minmod(u1,u2,u3);
	double u;
	if (i<0) u = b+deltau;
	if (i>=0) u = b-deltau;

	return u;
}
