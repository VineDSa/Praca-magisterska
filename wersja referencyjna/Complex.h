#pragma once

#include <iostream>

//#include "cuda_runtime.h"
//#include "device_launch_parameters.h"
#include <cuComplex.h>

typedef cuDoubleComplex clDoubleComplex;

#define Z_MAKE(r, i)        make_cuDoubleComplex((r), (i))
#define Z_REAL(a)           cuCreal(a)
#define Z_IMAG(a)           cuCimag(a)
#define Z_ADD(a, b)         cuCadd((a), (b))
#define Z_SUB(a, b)         cuCsub((a), (b))
#define Z_MUL(a, b)         cuCmul((a), (b))
#define Z_DIV(a, b)         cuCdiv((a), (b))
#define Z_ABS(a)            cuCabs((a))
#define Z_ZERO              make_cuDoubleComplex(0.0, 0.0)
#define Z_ONE               make_cuDoubleComplex(1.0, 0.0)
#define Z_ONE_IMG           make_cuDoubleComplex(0.0, 1.0)
#define Z_HALF              make_cuDoubleComplex(0.5, 0.0)
#define Z_NEG_ONE           make_cuDoubleComplex(-1.0, 0.0)
#define Z_NEG_HALF          make_cuDoubleComplex(-0.5, 0.0)
/*#define Z_MAKE(r, i)        make_clDoubleComplex((r), (i))
#define Z_REAL(a)           clCreal(a)
#define Z_IMAG(a)           clCimag(a)
#define Z_ADD(a, b)         clCadd((a), (b))
#define Z_SUB(a, b)         clCsub((a), (b))
#define Z_MUL(a, b)         clCmul((a), (b))
#define Z_MULD(a, b)		clCmulD((a), (b))
#define Z_DIV(a, b)         clCdiv((a), (b))
#define Z_ABS(a)            clCabs((a))
#define Z_CONJ(a)			clConj(a)
#define Z_ZERO              make_clDoubleComplex(0.0, 0.0)
#define Z_ONE               make_clDoubleComplex(1.0, 0.0)
#define Z_HALF              make_clDoubleComplex(0.5, 0.0)
#define Z_NEG_ONE           make_clDoubleComplex(-1.0, 0.0)
#define Z_NEG_HALF          make_clDoubleComplex(-0.5, 0.0) */

using namespace std;

/*struct clDoubleComplex
{
	double real;
	double imag;
}; */

clDoubleComplex make_clDoubleComplex(double r, double i);
double clCreal(clDoubleComplex x);
double clCimag(clDoubleComplex x);
clDoubleComplex clCadd(clDoubleComplex x, clDoubleComplex y);
clDoubleComplex clCsub(clDoubleComplex x, clDoubleComplex y);
clDoubleComplex clCmul(clDoubleComplex x, clDoubleComplex y);
clDoubleComplex clCmulD(clDoubleComplex x, double y);
clDoubleComplex clCdiv(clDoubleComplex x, clDoubleComplex y);
clDoubleComplex clConj(clDoubleComplex x);
double clCabs(clDoubleComplex x);