#pragma once
#include <iostream>
#include "Complex.h"

using namespace std;

void zmatrix(clDoubleComplex *cz, float *x, float *y, float *z,
	 int *ipc, int *ips, int *ise, float *rad, int npls, double ka, int nsmax, int ienv,
	float sgnx, float sgny, float sgnz, float sgnenv);

void wim(clDoubleComplex *cz, float *x, float *y, float *z, 
	int *ipc, int *ips, int *ise, float *rad, int npls, double ka, double eta0, 
	int nsmax, int ienv, float sgnx, float sgny, float sgnz, float sgnenv);

void cwire1_re(float *x, float *y, float *z, int *ipc, int *ips, int *ise, float *rad,
	int ms, int n, int nsmax, int ifun, double rmcx, double rmcy, double rmcz, double ka,
	double *dcx, double *dcy, double *dcz, double *dl1, clDoubleComplex *cval);

void cwire2_re(float *x, float *y, float *z, int *ipc, int *ips, int *ise, float *rad,
	int ms, int n, int nsmax, int ifun, double rmcx, double rmcy, double rmcz, double ka,
	double *dcx1, double *dcy1, double *dcz1, double *dl1, clDoubleComplex *cval);

void cwire1_im(float *x, float *y, float *z, int *ipc, int *ips, int *ise, float *rad,
	int ms, int n, int nsmax, int ifun, double rmcx, double rmcy, double rmcz, double ka,
	double *dcx1, double *dcy1, double *dcz1, double *dl1, float sgnx, float sgny, float sgnz, clDoubleComplex *cval);

void cwire2_im(float *x, float *y, float *z, int *ipc, int *ips, int *ise, float *rad,
	int ms, int n, int nsmax, int ifun, double rmcx, double rmcy, double rmcz, double ka,
	double *dcx1, double *dcy1, double *dcz1, double *dl1, float sgnx, float sgny, float sgnz, clDoubleComplex *cval);

clDoubleComplex ckmn(double s, int ifun, double ka, double rmcx, double rmcy, double rmcz,
	double rnx, double rny, double rnz, double dcx, double dcy, double dcz, double an, double dl);

clDoubleComplex cknn(double s, int ifun, double ka, double dl, double an);

void delick(double bet, double *res);