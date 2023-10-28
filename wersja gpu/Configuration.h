#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include<stdio.h>
#include<math.h>

using namespace std;

void geoinp(int nnmax, int nwmax, int ntmax, int &nnod, float *xnod, float *ynod, float *znod, \
	int &nw, int *iwn, float *wrad, int *nsw, int &nt, int *iwt);

void config(int &npnt, int nnod, float *xnod, float *ynod, float *znod, float *x, float *y, float *z, int nw, \
	int *iwn, float *wrad, int *nsw, int *ise, int *isp, int *ips, int *ipc, int *iwe, float *rad, \
	int &npls, int &nplsw, int &nseg, int nsmax, int nwmax, int ienv);

void configwire(int &npnt, int nnod, float *xnod, float *ynod, float *znod, float *x, float *y, float *z, int nw, \
	int *iwn, float *wrad, int *nsw, int *ise, int *isp, int *ips, int *ipc, int *iwe, float *rad, \
	int &nplsw, int &nseg, int nsmax, int nwmax, int ienv);

void tablice(int npnt, int nseg, int nsmax, float *x, float *y, float *z, int nplsw, int *ise, int *isp, int *ips, int *ipc);

/*
void config(int nnod, float *xnod, float *ynod, float *znod, float *x, float *y, float *z, \
	int nt, int *iwt, int *ipt, int *itp, int *iti, int &npls, int &nplsb, int npbmax);

//void configbody(int npntw, int nnod, float *xnod, float *ynod, float *znod, float *x, float *y, float *z, \
	int nt, int *iwt, int *ipt, int *itp, int *iti, int &nplsb, int npbmax);
//void freetriangle(float *x, float *y, float *z, float *xnod, float *ynod, float *znod, \
	int *iwt, int *iti, int nt, int m, int &ii);
//void tcon(int k11, int k12, int k21, int k22, int l, int m, int n, int &ntv, int &ntp, int np, int nt, \
	int *iti, int *itp);
//void newbodypulse(int npbmax, int &np, int nt, int ntv, int ntp, int m, int *ipt);
//void tablice(int nnod, int nw, int nt, int npls, int nplsb, int npbmax, float *x, float *y, float *z, \
	int *ipt, int *itp, int *iti);
	*/
