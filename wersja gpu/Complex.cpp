#include "Complex.h"

clDoubleComplex make_clDoubleComplex(double r, double i)
{
	clDoubleComplex res;
	//res.real = r;
	//res.imag = i;
	res.x = r;
	res.y = i;
	return res;
}

double clCreal(clDoubleComplex x)
{
	//return x.real;
	return x.x;
}

double clCimag(clDoubleComplex x)
{
	//return x.imag;
	return x.y;
}

clDoubleComplex clCadd(clDoubleComplex x, clDoubleComplex y)
{
	//return make_clDoubleComplex(x.real + y.real, x.imag + y.imag);
	return make_clDoubleComplex(x.x + y.x, x.y + y.y);
}

clDoubleComplex clCsub(clDoubleComplex x, clDoubleComplex y)
{
	//return make_clDoubleComplex(x.real - y.real, x.imag - y.imag);
	return make_clDoubleComplex(x.x - y.x, x.y - y.y);
}

clDoubleComplex clCmul(clDoubleComplex x, clDoubleComplex y)
{
	return make_clDoubleComplex((clCreal(x) * clCreal(y)) -(clCimag(x) * clCimag(y)),(clCreal(x) * clCimag(y)) +(clCimag(x) * clCreal(y)));
}

clDoubleComplex clCmulD(clDoubleComplex x, double y)
{
	return make_clDoubleComplex(clCreal(x)*y, clCimag(x)*y);
}

clDoubleComplex clCdiv(clDoubleComplex x, clDoubleComplex y)
{
	clDoubleComplex quot;
	double s = (abs(clCreal(y))) + (abs(clCimag(y)));
	double oos = 1.0 / s;
	double ars = clCreal(x) * oos;
	double ais = clCimag(x) * oos;
	double brs = clCreal(y) * oos;
	double bis = clCimag(y) * oos;
	s = (brs * brs) + (bis * bis);
	oos = 1.0 / s;
	quot = make_clDoubleComplex(((ars * brs) + (ais * bis)) * oos,
		((ais * brs) - (ars * bis)) * oos);
	return quot;
}

clDoubleComplex clConj(clDoubleComplex x)
{
	return make_clDoubleComplex(clCreal(x), -clCimag(x));
}

double clCabs(clDoubleComplex x)
{
	double a = clCreal(x);
	double b = clCimag(x);
	double v, w, t;
	a = fabs(a);
	b = fabs(b);
	if (a > b) {
		v = a;
		w = b;
	}
	else {
		v = b;
		w = a;
	}
	t = w / v;
	t = 1.0 + t * t;
	t = v * sqrt(t);
	if ((v == 0.0) ||
		(v > 1.79769313486231570e+308) || (w > 1.79769313486231570e+308)) {
		t = v + w;
	}
	return t;
}