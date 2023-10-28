#pragma once
#include "Zmatrix.h"
#include <chrono>
#include <iomanip>

//----------------------------------------------------------------------
struct Timer
{
	chrono::time_point<chrono::high_resolution_clock> start, end;
	chrono::duration<float> duration;
	char* info = "None";
	int process = 9999;
	Timer()
	{
		start = chrono::high_resolution_clock::now();
	}

	~Timer()
	{
		end = chrono::high_resolution_clock::now();
		duration = end - start;
		float s = duration.count();
		std::cout << std::endl << " Process - " << process << " " << info << " took " << std::setprecision(4) << s << "s" << std::endl;
	}
};

void zmatrix(clDoubleComplex *cz, float *x, float *y, float *z,
	int *ipc, int *ips, int *ise, float *rad, int npls, double ka, int nsmax, int ienv,
	float sgnx, float sgny, float sgnz, float sgnenv)
{
	Timer timer;
	timer.process = 0;
	timer.info = "zmatrix on CPU";

	double pi = acos(-1.0);
	double eta0 = (double)120.*pi;
	clDoubleComplex cjq = Z_MAKE(0.0, eta0 / (4.0*pi*ka));

	for (int i = 0; i < npls; i++)
		for (int j = 0; j < npls; j++)    
		{
			cz[i + j * npls] = Z_ZERO;
		}

	wim(cz, x, y, z, ipc, ips, ise, rad, npls, ka, eta0, nsmax, ienv, sgnx, sgny, sgnz, sgnenv);

	for (int i = 0; i < npls; i++)
		for (int j = 0; j < npls; j++)
		{
			cz[i + j * npls] = Z_MUL(cjq, cz[i + j * npls]);
		}
}

//----------------------------------------------------------------------

void wim(clDoubleComplex *cz, float *x, float *y, float *z,
	int *ipc, int *ips, int *ise, float *rad, int npls, double ka, double eta0, 
	int nsmax, int ienv, float sgnx, float sgny, float sgnz, float sgnenv)
{
	int mp, ms, ml, mr, ifun;
	double rmpx, rmpy, rmpz, rmlx, rmly, rmlz, rmrx, rmry, rmrz;
	double am, rmcx, rmcy, rmcz, lmcx, lmcy, lmcz;
	double tvecl, dcx, dcy, dcz, dl;
	clDoubleComplex cvval1, cvval2, csval1, csval2, cvval, csval, cpom;
	
	for (int m = 0; m < npls; m++)
	{
		//punkt obserwacji na pierwszym segmencie
		mp = abs(ipc[m]);
		rmpx = (double)x[mp-1];
		rmpy = (double)y[mp-1];
		rmpz = (double)z[mp-1];
		ms = ips[m];
		ml = ise[ms-1];
		rmlx = (double)x[ml-1];
		rmly = (double)y[ml-1];
		rmlz = (double)z[ml-1];
		am = rad[ms-1];						  //promien Sm+   
		rmcx = 0.5*(rmpx + rmlx);             //wspolrzedne punktu obserwacji na Sm +
		rmcy = 0.5*(rmpy + rmly);
		rmcz = 0.5*(rmpz + rmlz);
		lmcx = rmpx - rmcx;                   //wektor lm +
		lmcy = rmpy - rmcy;
		lmcz = rmpz - rmcz;

		for (int n = 0; n < npls; n++)
		{
			ifun = 1;
			cwire1_re(x, y, z, ipc, ips, ise, rad, ms, n, nsmax, ifun, rmcx, rmcy, rmcz, ka, &dcx, &dcy, &dcz, &dl, &cvval1);  // obliczanie calki na Sn+  pot.wek.			
			tvecl = dcx * lmcx + dcy*lmcy + dcz*lmcz;      //obliczane  !lm + *1n +
			cvval = Z_MUL(cvval1,Z_MAKE(tvecl,(double)0.e0));  //dcmplx(tvecl, 0.d0)*cvval1

			//!sk³adowa potencja³u A od obrazu dla cwire 1       PEC
			if (ienv == 2)
			{
				//call cwire1_im(ms, n, cvval1)
				cwire1_im(x, y, z, ipc, ips, ise, rad, ms, n, nsmax, ifun, rmcx, rmcy, rmcz, ka, &dcx, &dcy, &dcz, &dl, sgnx, sgny, sgnz, &cvval1);
				tvecl = dcx * lmcx + dcy * lmcy + dcz * lmcz;
				//cvval = cvval + sgnenv * dcmplx(tvecl, 0.d0)*cvval1
				tvecl = tvecl * (double)sgnenv;
				cvval = Z_ADD(cvval, Z_MUL(cvval1, Z_MAKE(tvecl, (double)0.e0)));				
			}		    

			ifun = 2;
			cwire2_re(x, y, z, ipc, ips, ise, rad, ms, n, nsmax, ifun, rmcx, rmcy, rmcz, ka, &dcx, &dcy, &dcz, &dl, &cvval2);  // obliczanie calki na Sn-  pot.wek.
			tvecl = dcx * lmcx + dcy * lmcy + dcz * lmcz;       //obliczane  !lm + *1n -
			//	cvval = cvval + dcmplx(tvecl, 0.d0)*cvval2   //A(rm + )
			cvval = Z_ADD(cvval, Z_MUL(cvval2, Z_MAKE(tvecl, (double)0.e0)));

			//sk³adowa potencja³u A od obrazu dla cwire 2       PEC
			if (ienv == 2)
			{
				//call cwire2_im(ms, n, cvval2)
				cwire2_im(x, y, z, ipc, ips, ise, rad, ms, n, nsmax, ifun, rmcx, rmcy, rmcz, ka, &dcx, &dcy, &dcz, &dl, sgnx, sgny, sgnz, &cvval2);
				tvecl = dcx * lmcx + dcy * lmcy + dcz * lmcz;
				//cvval = cvval + sgnenv * dcmplx(tvecl, 0.d0)*cvval2
				tvecl = tvecl * (double)sgnenv;
				cvval = Z_ADD(cvval, Z_MUL(cvval2, Z_MAKE(tvecl, (double)0.e0)));
			}

			ifun = 3;
			cwire1_re(x, y, z, ipc, ips, ise, rad, ms, n, nsmax, ifun, rmcx, rmcy, rmcz, ka, &dcx, &dcy, &dcz, &dl, &csval1);
			//csval = csval1 / dl;
			csval = Z_DIV(csval1, Z_MAKE(dl, (double)0.e0));
			//sk³adowa potencja³u skalarnego od cwire 1     PEC
			if (ienv==2)
			{
				cwire1_im(x, y, z, ipc, ips, ise, rad, ms, n, nsmax, ifun, rmcx, rmcy, rmcz, ka, &dcx, &dcy, &dcz, &dl, sgnx, sgny, sgnz, &csval1);
				//csval = csval + sgnenv * csval1 / dl
				dl = (double)sgnenv / dl;
				csval = Z_ADD(csval, Z_MUL(csval1, Z_MAKE(dl, (double)0.e0)));
			}
				
			cwire2_re(x, y, z, ipc, ips, ise, rad, ms, n, nsmax, ifun, rmcx, rmcy, rmcz, ka, &dcx, &dcy, &dcz, &dl, &csval2);
			//csval = csval - csval2 / dl       !O(rm + )
			csval = Z_SUB(csval, Z_DIV(csval2, Z_MAKE(dl, (double)0.e0)));

			//sk³adowa potencja³u skalarnego od cwire 2     PEC
			if (ienv == 2)
			{
				cwire2_im(x, y, z, ipc, ips, ise, rad, ms, n, nsmax, ifun, rmcx, rmcy, rmcz, ka, &dcx, &dcy, &dcz, &dl, sgnx, sgny, sgnz, &csval2);
				//csval = csval - sgnenv * csval2 / dl
				dl = (double)sgnenv / dl;
				csval = Z_SUB(csval, Z_MUL(csval2, Z_MAKE(dl, (double)0.e0)));
			}
						
			//cz(mm, nn) = cz(mm, nn) + ka * ka*cvval - csval			
			cpom = Z_MUL(cvval, Z_MAKE(ka*ka, (double)0.e0));
			cpom = Z_SUB(cpom, csval);
			cz[m + n * npls] = Z_ADD(cz[m + n * npls], cpom);
		}

		// punkt obserwacji na drugim segmencie
		ms = ips[m + nsmax];
		mr = ise[ms + nsmax - 1];
		rmrx = (double)x[mr - 1];
		rmry = (double)y[mr - 1];
		rmrz = (double)z[mr - 1];
		am = (double)rad[ms - 1];                //promien Sm -
		rmcx = (double)0.5e0*(rmrx + rmpx);             //wspolrzedne punktu obserwacji na Sm -
		rmcy = (double)0.5e0*(rmry + rmpy);
		rmcz = (double)0.5e0*(rmrz + rmpz);
		lmcx = rmcx - rmpx;                     //wektor lm -
		lmcy = rmcy - rmpy;
		lmcz = rmcz - rmpz;

		for (int n = 0; n < npls; n++)
		{
			ifun = 1;
			cwire1_re(x, y, z, ipc, ips, ise, rad, ms, n, nsmax, ifun, rmcx, rmcy, rmcz, ka, &dcx, &dcy, &dcz, &dl, &cvval1);  // obliczanie calki na Sn+  pot.wek.
			tvecl = dcx * lmcx + dcy * lmcy + dcz * lmcz;      //obliczane  !lm - *1n +
			cvval = Z_MUL(cvval1, Z_MAKE(tvecl, (double)0.e0));  //dcmplx(tvecl, 0.d0)*cvval1
								 
			if (ienv == 2)				//!sk³adowa potencja³u A od obrazu dla cwire 1       PEC
			{
				//call cwire1_im(ms, n, cvval1)
				cwire1_im(x, y, z, ipc, ips, ise, rad, ms, n, nsmax, ifun, rmcx, rmcy, rmcz, ka, &dcx, &dcy, &dcz, &dl, sgnx, sgny, sgnz, &cvval1);
				tvecl = dcx * lmcx + dcy * lmcy + dcz * lmcz;
				//cvval = cvval + sgnenv * dcmplx(tvecl, 0.d0)*cvval1
				tvecl = tvecl * (double)sgnenv;
				cvval = Z_ADD(cvval, Z_MUL(cvval1, Z_MAKE(tvecl, (double)0.e0)));
			}

			ifun = 2;
			cwire2_re(x, y, z, ipc, ips, ise, rad, ms, n, nsmax, ifun, rmcx, rmcy, rmcz, ka, &dcx, &dcy, &dcz, &dl, &cvval2);  // obliczanie calki na Sn-  pot.wek.
			tvecl = dcx * lmcx + dcy * lmcy + dcz * lmcz;       //obliczane  !lm - *1n -
			//	cvval = cvval + dcmplx(tvecl, 0.d0)*cvval2   //A(rm - )
			cvval = Z_ADD(cvval, Z_MUL(cvval2, Z_MAKE(tvecl, (double)0.e0)));
			
			if (ienv == 2)			//sk³adowa potencja³u A od obrazu dla cwire 2       PEC
			{
				//call cwire2_im(ms, n, cvval2)
				cwire2_im(x, y, z, ipc, ips, ise, rad, ms, n, nsmax, ifun, rmcx, rmcy, rmcz, ka, &dcx, &dcy, &dcz, &dl, sgnx, sgny, sgnz, &cvval2);
				tvecl = dcx * lmcx + dcy * lmcy + dcz * lmcz;
				//cvval = cvval + sgnenv * dcmplx(tvecl, 0.d0)*cvval2
				tvecl = tvecl * (double)sgnenv;
				cvval = Z_ADD(cvval, Z_MUL(cvval2, Z_MAKE(tvecl, (double)0.e0)));
			}

			ifun = 3;
			cwire1_re(x, y, z, ipc, ips, ise, rad, ms, n, nsmax, ifun, rmcx, rmcy, rmcz, ka, &dcx, &dcy, &dcz, &dl, &csval1);
			//csval = csval1 / dl;
			csval = Z_DIV(csval1, Z_MAKE(dl, (double)0.e0));
			
			if (ienv == 2)		//sk³adowa potencja³u skalarnego od cwire 1     PEC
			{
				cwire1_im(x, y, z, ipc, ips, ise, rad, ms, n, nsmax, ifun, rmcx, rmcy, rmcz, ka, &dcx, &dcy, &dcz, &dl, sgnx, sgny, sgnz, &csval1);
				//csval = csval + sgnenv * csval1 / dl
				dl = (double)sgnenv / dl;
				csval = Z_ADD(csval, Z_MUL(csval1, Z_MAKE(dl, (double)0.e0)));
			}

			cwire2_re(x, y, z, ipc, ips, ise, rad, ms, n, nsmax, ifun, rmcx, rmcy, rmcz, ka, &dcx, &dcy, &dcz, &dl, &csval2);
			//csval = csval - csval2 / dl       !O(rm - )
			csval = Z_SUB(csval, Z_DIV(csval2, Z_MAKE(dl, (double)0.e0)));

			//sk³adowa potencja³u skalarnego od cwire 2     PEC
			if (ienv == 2)
			{
				cwire2_im(x, y, z, ipc, ips, ise, rad, ms, n, nsmax, ifun, rmcx, rmcy, rmcz, ka, &dcx, &dcy, &dcz, &dl, sgnx, sgny, sgnz, &csval2);
				//csval = csval - sgnenv * csval2 / dl
				dl = (double)sgnenv / dl;
				csval = Z_SUB(csval, Z_MUL(csval2, Z_MAKE(dl, (double)0.e0)));
			}

			//cz(mm, nn) = cz(mm, nn) + ka * ka*cvval + csval			
			cpom = Z_MUL(cvval, Z_MAKE(ka*ka, (double)0.e0));
			cpom = Z_ADD(cpom, csval);
			cz[m + n * npls] = Z_ADD(cz[m + n * npls], cpom);

		}
	}
}

//----------------------------------------------------------------------

void cwire1_re(float *x, float *y, float *z, int *ipc, int *ips, int *ise, float *rad, 
	int ms, int n, int nsmax, int ifun, double rmcx, double rmcy, double rmcz, double ka,
	double *dcx1, double *dcy1, double *dcz1, double *dl1, clDoubleComplex *cval)
{
	int np, ns, nl, nr, ml, mr;
	double rnpx, rnpy, rnpz, an;
	double rnx, rny, rnz, rdx,rdy,rdz;
	double dl, dlh;
	double res;
	double xl, xu, a, b, c, s1, s2;
	clDoubleComplex cv, cv1, cv2, cvpom;
	double pi = acos((double)-1.0);
	double dcx, dcy, dcz;
	*cval = Z_ZERO;

	np = abs(ipc[n]);                 //na  jakim wezle n - ta funkcja bazowa
	rnpx = (double)x[np - 1];                 //wspolrzedne tego wezla
	rnpy = (double)y[np - 1];
	rnpz = (double)z[np - 1];
	ns = ips[n];
	nl = ise[ns - 1];                  //pierwszy wezel pierwszego sagmentu
	nr = ise[ns + nsmax - 1];
	ml = ise[ms - 1];
	mr = ise[ms + nsmax - 1];
	an = (double)rad[ns - 1];              //na ktorym n - ta f.bazowa
	rnx = (double)x[nl - 1];                //wspolrzedne pierwszego wezela pierwszego
	rny = (double)y[nl - 1];                //sagmentu na ktorym n - ta f.bazowa
	rnz = (double)z[nl - 1];
	rdx = rnpx - rnx;
	rdy = rnpy - rny;
	rdz = rnpz - rnz;
	dl = sqrt(rdx*rdx + rdy * rdy + rdz * rdz);
	*dl1 = dl;
	dlh = 0.5*dl;
	dcx = rdx / dl;
	dcy = rdy / dl;
	dcz = rdz / dl;
	*dcx1 = dcx;
	*dcy1 = dcy;
	*dcz1 = dcz;
	if (((nl == ml) && (nr == mr)) || ((nl == mr) && (nr == ml)))
	{
		// osobliwoœæ - dwa ca³kowania po kwadraturze 4 punktowej 
		//dcqg(0.0, dlh, cknn, cv1, 4);		
		xl = (double)0.e0;
		xu = dlh;
		a = (double)0.5e0*(xl + xu);
		b = xu - xl;
		c = (double)0.43056815579702629e0*b;
		s1 = a + c;
		s2 = a - c;
		//cout << "s1: " << s1 << "   s2:" << s2 << endl;
		cv1 = Z_ADD(cknn(s1, ifun, ka, dl, an), cknn(s2, ifun, ka, dl, an));
		cv1 = Z_MUL(cv1, Z_MAKE((double)0.17392742256872693e0, (double)0.e0));
		c = (double)0.16999052179242813e0*b;
		s1 = a + c;
		s2 = a - c;
		cvpom= Z_ADD(cknn(s1, ifun, ka, dl, an), cknn(s2, ifun, ka, dl, an));
		//y = b * (y + .32607257743127307d0*(f(a + c) + f(a - c)))
		cvpom = Z_MUL(cvpom, Z_MAKE((double)0.32607257743127307e0, (double)0.e0));
		cvpom = Z_ADD(cv1, cvpom);
		cv1 = Z_MUL(cvpom, Z_MAKE(b, (double)0.e0));		
		
		//dcqg(dlh, dl, cknn, cv2, 4);
		xl = dlh;
		xu = dl;
		a = (double)0.5e0*(xl + xu);
		b = xu - xl;
		c = (double)0.43056815579702629e0*b;
		s1 = a + c;
		s2 = a - c;
		cv2 = Z_ADD(cknn(s1, ifun, ka, dl, an), cknn(s2, ifun, ka, dl, an));
		cv2 = Z_MUL(cv2, Z_MAKE((double)0.17392742256872693e0, (double)0.e0));
		c = (double)0.16999052179242813e0*b;
		s1 = a + c;
		s2 = a - c;
		cvpom = Z_ADD(cknn(s1, ifun, ka, dl, an), cknn(s2, ifun, ka, dl, an));
		//y = b * (y + .32607257743127307d0*(f(a + c) + f(a - c)))
		cvpom = Z_MUL(cvpom, Z_MAKE((double)0.32607257743127307e0, (double)0.e0));
		cvpom = Z_ADD(cv2, cvpom);
		cv2 = Z_MUL(cvpom, Z_MAKE(b, (double)0.e0));				

		if ((ifun == 1) || (ifun == 2)) res = (double)0.5e0*dl*((double)1.e0 + log((double)16.e0*an / dl)) / (pi*an);
		if (ifun == 3) res = dl * ((double)1.e0 + log((double)16.e0*an / dl)) / (pi*an);
		//cv = cv1 + cv2 + Z_MAKE(res, (double)0.e0);
		cv = Z_ADD(cv1, cv2);
		cv = Z_ADD(cv, Z_MAKE(res, (double)0.e0));		
	}
	else
	{
		// ca³kowanie - kwadratura dwupunktowa   //dcqg(0.0, dl, ckmn, cv, 2);
		xl = (double)0.e0;
		xu = dl;
		a = (double)0.5e0*(xl + xu);
		b = xu - xl;
		c = (double)0.288675134594812882e0*b;
		s1 = a + c;
		s2 = a - c;
		cv = Z_ADD(ckmn(s1, ifun, ka, rmcx, rmcy, rmcz, rnx, rny, rnz, dcx, dcy, dcz, an, dl), ckmn(s2, ifun, ka, rmcx, rmcy, rmcz, rnx, rny, rnz, dcx, dcy, dcz, an, dl));
		cv = Z_MUL(cv, Z_MAKE(b * (double)0.5e0, (double)0.e0));
	}
	*cval = cv;
}
//----------------------------------------------------------------------

void cwire1_im(float *x, float *y, float *z, int *ipc, int *ips, int *ise, float *rad,
	int ms, int n, int nsmax, int ifun, double rmcx, double rmcy, double rmcz, double ka,
	double *dcx1, double *dcy1, double *dcz1, double *dl1, float sgnx, float sgny, float sgnz, clDoubleComplex *cval)
{
	//Obliczanie ca³ek potencja³owych dla segmentów Sn +
	//(1 segment na którym rozci¹gniêta n - ta f.bazowa)
	double dl, dlh;
	double res;
	int np, ns, nl, nr, ml, mr;
	double rnpx, rnpy, rnpz, an;
	double rnx, rny, rnz, rdx, rdy, rdz;
	double xl, xu, a, b, c, s1, s2;
	clDoubleComplex cv, cv1, cv2, cvpom;
	double pi = acos((double)-1.0);
	double dcx, dcy, dcz;

	*cval = Z_ZERO;	
	np = abs(ipc[n]);                 //na  jakim wezle n - ta funkcja bazowa
	rnpx = (double)sgnx * (double)x[np - 1];                 //wspolrzedne tego wezla
	rnpy = (double)sgny * (double)y[np - 1];
	rnpz = (double)sgnz * (double)z[np - 1];
	ns = ips[n];
	nl = ise[ns - 1];                  //pierwszy wezel pierwszego sagmentu
	nr = ise[ns + nsmax - 1];                  //drugi wezel pierwszego sagmentu
	ml = ise[ms - 1];
	mr = ise[ms + nsmax - 1];
	an = (double)rad[ns - 1];              //na ktorym n - ta f.bazowa
	rnx = (double)sgnx * (double)x[nl - 1];                //wspolrzedne pierwszego wezela pierwszego
	rny = (double)sgny * (double)y[nl - 1];                //sagmentu na ktorym n - ta f.bazowa
	rnz = (double)sgnz * (double)z[nl - 1];
	rdx = rnpx - rnx;
	rdy = rnpy - rny;
	rdz = rnpz - rnz;
	dl = sqrt(rdx*rdx + rdy * rdy + rdz * rdz);
	*dl1 = dl;
	dlh = (double)0.5e0*dl;
	dcx = rdx / dl;
	dcy = rdy / dl;
	dcz = rdz / dl;
	*dcx1 = dcx;
	*dcy1 = dcy;
	*dcz1 = dcz;
	if ((((nl == ml) && (nr == mr)) && (z[ml-1] < 1.e-6 && z[mr-1] < 1.e-6)) || (((nl == mr) && (nr == ml)) && (z[ml-1] < 1.e-6 && z[mr-1] < 1.e-6)))
		//if ((nl.eq.ml.and.nr.eq.mr).and.z(ml).eq.0..and.z(mr).eq.0.) goto 20
		//if ((nl.eq.mr.and.nr.eq.ml).and.z(ml).eq.0..and.z(mr).eq.0.) goto 20
	{
		//call dcqg(0.d0, dlh, cknn, cv1, 4)
		xl = (double)0.e0;
		xu = dlh;
		a = (double)0.5e0*(xl + xu);
		b = xu - xl;
		c = (double)0.43056815579702629e0*b;
		s1 = a + c;
		s2 = a - c;
		cv1 = Z_ADD(cknn(s1, ifun, ka, dl, an), cknn(s2, ifun, ka, dl, an));
		cv1 = Z_MUL(cv1, Z_MAKE((double)0.17392742256872693e0, (double)0.e0));
		c = (double)0.16999052179242813e0*b;
		s1 = a + c;
		s2 = a - c;
		cvpom = Z_ADD(cknn(s1, ifun, ka, dl, an), cknn(s2, ifun, ka, dl, an));
		//y = b * (y + .32607257743127307d0*(f(a + c) + f(a - c)))
		cvpom = Z_MUL(cvpom, Z_MAKE((double)0.32607257743127307e0, (double)0.e0));
		cvpom = Z_ADD(cv1, cvpom);
		cv1 = Z_MUL(cvpom, Z_MAKE(b, (double)0.e0));		

		//call dcqg(dlh, dl, cknn, cv2, 4)
		xl = dlh;
		xu = dl;
		a = (double)0.5e0*(xl + xu);
		b = xu - xl;
		c = (double)0.43056815579702629e0*b;
		s1 = a + c;
		s2 = a - c;
		cv2 = Z_ADD(cknn(s1, ifun, ka, dl, an), cknn(s2, ifun, ka, dl, an));
		cv2 = Z_MUL(cv2, Z_MAKE((double)0.17392742256872693e0, (double)0.e0));
		c = (double)0.16999052179242813e0*b;
		s1 = a + c;
		s2 = a - c;
		cvpom = Z_ADD(cknn(s1, ifun, ka, dl, an), cknn(s2, ifun, ka, dl, an));
		//y = b * (y + .32607257743127307d0*(f(a + c) + f(a - c)))
		cvpom = Z_MUL(cvpom, Z_MAKE((double)0.32607257743127307e0, (double)0.e0));
		cvpom = Z_ADD(cv2, cvpom);
		cv2 = Z_MUL(cvpom, Z_MAKE(b, (double)0.e0));		

		if ((ifun == 1) || (ifun == 2)) res = (double)0.5e0*dl*((double)1.e0 + log((double)16.e0*an / dl)) / (pi*an);
		if (ifun == 3) res = dl * ((double)1.e0 + log((double)16.e0*an / dl)) / (pi*an);
		//cv = cv1 + cv2 + Z_MAKE(res, (double)0.e0);
		cv = Z_ADD(cv1, cv2);
		cv = Z_ADD(cv, Z_MAKE(res, (double)0.e0));
	}
	else
	{
		// ca³kowanie - kwadratura dwupunktowa   //dcqg(0.0, dl, ckmn, cv, 2);
		xl = (double)0.e0;
		xu = dl;
		a = (double)0.5e0*(xl + xu);
		b = xu - xl;
		c = (double)0.288675134594812882e0*b;
		s1 = a + c;
		s2 = a - c;
		cv = Z_ADD(ckmn(s1, ifun, ka, rmcx, rmcy, rmcz, rnx, rny, rnz, dcx, dcy, dcz, an, dl), ckmn(s2, ifun, ka, rmcx, rmcy, rmcz, rnx, rny, rnz, dcx, dcy, dcz, an, dl));
		cv = Z_MUL(cv, Z_MAKE(b * (double)0.5e0, (double)0.e0));
	}
	*cval = cv;
}

//----------------------------------------------------------------------

void cwire2_re(float *x, float *y, float *z, int *ipc, int *ips, int *ise, float *rad,
	int ms, int n, int nsmax, int ifun, double rmcx, double rmcy, double rmcz, double ka,
	double *dcx1, double *dcy1, double *dcz1, double *dl1, clDoubleComplex *cval)
{
	int ml, mr, np, ns, nl, nr;

	double rnrx, rnry, rnrz, an;
	double rnx, rny, rnz, rdx, rdy, rdz;
	double dl, dlh;
	double res;
	double xl, xu, a, b, c, s1, s2;
	clDoubleComplex cv, cv1, cv2, cvpom;
	double pi = acos((double)-1.0);
	double dcx, dcy, dcz;
	
	*cval = Z_ZERO;
	ml = ise[ms - 1];
	mr = ise[ms + nsmax - 1];
	np = abs(ipc[n]);                //na  jakim wezle n - ta funkcja bazowa(wire)
	ns = ips[n + nsmax];
	nl = ise[ns - 1];
	nr = ise[ns + nsmax - 1];				//drugi wezel drugiego sagmentu
	rnx = (double)x[np - 1];                 //wspolrzedne tego wezla
	rny = (double)y[np - 1];                 //na ktorym n - ta f.bazowa
	rnz = (double)z[np - 1];
	an = (double)rad[ns - 1];
	rnrx = (double)x[nr - 1];                //wspolrzedne drugiego wezela drugiego sagmentu
	rnry = (double)y[nr - 1];                //na ktorym n - ta f.bazowa
	rnrz = (double)z[nr - 1];
	rdx = rnrx - rnx;
	rdy = rnry - rny;
	rdz = rnrz - rnz;
	dl = sqrt(rdx*rdx + rdy * rdy + rdz * rdz);
	*dl1 = dl;
	dlh = (double)0.5e0*dl;
	dcx = rdx / dl;
	dcy = rdy / dl;
	dcz = rdz / dl;
	*dcx1 = dcx;
	*dcy1 = dcy;
	*dcz1 = dcz;

	if (((nl == ml) && (nr == mr)) || ((nl == mr) && (nr == ml)))
		// if (nl.eq.ml.and.nr.eq.mr) goto 20
		// if (nl.eq.mr.and.nr.eq.ml) goto 20
	{
		// osobliwoœæ - dwa ca³kowania po kwadraturze 4 punktowej 		
		//dcqg(0.0, dlh, cknn, cv1, 4);		
		xl = (double)0.e0;
		xu = dlh;
		a = (double)0.5e0*(xl + xu);
		b = xu - xl;
		c = (double)0.43056815579702629e0*b;
		s1 = a + c;
		s2 = a - c;
		cv1 = Z_ADD(cknn(s1, ifun, ka, dl, an), cknn(s2, ifun, ka, dl, an));
		cv1 = Z_MUL(cv1, Z_MAKE((double)0.17392742256872693e0, (double)0.e0));
		c = (double)0.16999052179242813e0*b;
		s1 = a + c;
		s2 = a - c;
		cvpom = Z_ADD(cknn(s1, ifun, ka, dl, an), cknn(s2, ifun, ka, dl, an));
		//y = b * (y + .32607257743127307d0*(f(a + c) + f(a - c)))
		cvpom = Z_MUL(cvpom, Z_MAKE((double)0.32607257743127307e0, (double)0.e0));
		cvpom = Z_ADD(cv1, cvpom);
		cv1 = Z_MUL(cvpom, Z_MAKE(b, (double)0.e0));		

		//dcqg(dlh, dl, cknn, cv2, 4);
		xl = dlh;
		xu = dl;
		a = (double)0.5e0*(xl + xu);
		b = xu - xl;
		c = (double)0.43056815579702629e0*b;
		s1 = a + c;
		s2 = a - c;
		cv2 = Z_ADD(cknn(s1, ifun, ka, dl, an), cknn(s2, ifun, ka, dl, an));
		cv2 = Z_MUL(cv2, Z_MAKE((double)0.17392742256872693e0, (double)0.e0));
		c = (double)0.16999052179242813e0*b;
		s1 = a + c;
		s2 = a - c;
		cvpom = Z_ADD(cknn(s1, ifun, ka, dl, an), cknn(s2, ifun, ka, dl, an));
		//y = b * (y + .32607257743127307d0*(f(a + c) + f(a - c)))
		cvpom = Z_MUL(cvpom, Z_MAKE((double)0.32607257743127307e0, (double)0.e0));
		cvpom = Z_ADD(cv2, cvpom);
		cv2 = Z_MUL(cvpom, Z_MAKE(b, (double)0.e0));		

		if ((ifun == 1) || (ifun == 2)) res = (double)0.5e0*dl*((double)1.e0 + log((double)16.e0*an / dl)) / (pi*an);
		if (ifun == 3) res = dl * ((double)1.e0 + log((double)16.e0*an / dl)) / (pi*an);
		//cv = cv1 + cv2 + Z_MAKE(res, (double)0.e0);
		cv = Z_ADD(cv1, cv2);
		cv = Z_ADD(cv, Z_MAKE(res, (double)0.e0));
	}
	else
	{
		// ca³kowanie - kwadratura dwupunktowa   //dcqg(0.0, dl, ckmn, cv, 2); 
		xl = (double)0.e0;
		xu = dl;
		a = (double)0.5e0*(xl + xu);
		b = xu - xl;
		c = (double)0.288675134594812882e0*b;
		s1 = a + c;
		s2 = a - c;
		cv = Z_ADD(ckmn(s1, ifun, ka, rmcx, rmcy, rmcz, rnx, rny, rnz, dcx, dcy, dcz, an, dl), ckmn(s2, ifun, ka, rmcx, rmcy, rmcz, rnx, rny, rnz, dcx, dcy, dcz, an, dl));
		cv = Z_MUL(cv, Z_MAKE(b * (double)0.5e0, (double)0.e0));
	}
	*cval = cv;
}

//----------------------------------------------------------------------

void cwire2_im(float *x, float *y, float *z, int *ipc, int *ips, int *ise, float *rad,
	int ms, int n, int nsmax, int ifun, double rmcx, double rmcy, double rmcz, double ka,
	double *dcx1, double *dcy1, double *dcz1, double *dl1, float sgnx, float sgny, float sgnz, clDoubleComplex *cval)
{
	//Obliczanie ca³ek potencja³owych dla segmentów Sn +
	//(1 segment na którym rozci¹gniêta n - ta f.bazowa)
	double dl, dlh;
	double res;
	int np, ns, nl, nr, ml, mr;
	double rnrx, rnry, rnrz, an;
	double rnx, rny, rnz, rdx, rdy, rdz;
	double xl, xu, a, b, c, s1, s2;
	clDoubleComplex cv, cv1, cv2, cvpom;
	double pi = acos((double)-1.0);
	double dcx, dcy, dcz;

	*cval = Z_ZERO;
	ml = ise[ms - 1];
	mr = ise[ms + nsmax - 1];
	np = abs(ipc[n]);                 //na  jakim wezle n - ta funkcja bazowa
	ns = ips[n + nsmax]; //ns = ips(n, 2)
	nl = ise[ns - 1];                  //pierwszy wezel drugiegoo sagmentu
	nr = ise[ns + nsmax - 1];                  //drugi wezel drugiego sagmentu
	rnx = (double)sgnx * (double)x[np - 1];                //wspolrzedne tego wezela
	rny = (double)sgny * (double)y[np - 1];                //na ktorym n - ta f.bazowa
	rnz = (double)sgnz * (double)z[np - 1];
	an = (double)rad[ns - 1];
	rnrx = (double)sgnx * (double)x[nr - 1];                 //! wspolrzedne drugiego wezela drugiego sagmentu
	rnry = (double)sgny * (double)y[nr - 1];				//na ktorym n-ta f.bazowa
	rnrz = (double)sgnz * (double)z[nr - 1];
	rdx = rnrx - rnx;
	rdy = rnry - rny;
	rdz = rnrz - rnz;
	dl = sqrt(rdx*rdx + rdy * rdy + rdz * rdz);
	*dl1 = dl;
	dlh = (double)0.5e0*dl;
	dcx = rdx / dl;
	dcy = rdy / dl;
	dcz = rdz / dl;
	*dcx1 = dcx;
	*dcy1 = dcy;
	*dcz1 = dcz;
	if ((((nl == ml) && (nr == mr)) && (z[ml-1] < 1.e-6 && z[mr-1] < 1.e-6)) || (((nl == mr) && (nr == ml)) && (z[ml-1] < 1.e-6 && z[mr-1] < 1.e-6)))
		//if ((nl.eq.ml.and.nr.eq.mr).and.z(ml).eq.0..and.z(mr).eq.0.) goto 20
		//if ((nl.eq.mr.and.nr.eq.ml).and.z(ml).eq.0..and.z(mr).eq.0.) goto 20
	{
		//call dcqg(0.d0, dlh, cknn, cv1, 4)
		xl = (double)0.e0;
		xu = dlh;
		a = (double)0.5e0*(xl + xu);
		b = xu - xl;
		c = (double)0.43056815579702629e0*b;
		s1 = a + c;
		s2 = a - c;
		cv1 = Z_ADD(cknn(s1, ifun, ka, dl, an), cknn(s2, ifun, ka, dl, an));
		cv1 = Z_MUL(cv1, Z_MAKE((double)0.17392742256872693e0, (double)0.e0));
		c = (double)0.16999052179242813e0*b;
		s1 = a + c;
		s2 = a - c;
		cvpom = Z_ADD(cknn(s1, ifun, ka, dl, an), cknn(s2, ifun, ka, dl, an));
		//y = b * (y + .32607257743127307d0*(f(a + c) + f(a - c)))
		cvpom = Z_MUL(cvpom, Z_MAKE((double)0.32607257743127307e0, (double)0.e0));
		cvpom = Z_ADD(cv1, cvpom);
		cv1 = Z_MUL(cvpom, Z_MAKE(b, (double)0.e0));		

		//call dcqg(dlh, dl, cknn, cv2, 4)
		xl = dlh;
		xu = dl;
		a = (double)0.5e0*(xl + xu);
		b = xu - xl;
		c = (double)0.43056815579702629e0*b;
		s1 = a + c;
		s2 = a - c;
		cv2 = Z_ADD(cknn(s1, ifun, ka, dl, an), cknn(s2, ifun, ka, dl, an));
		cv2 = Z_MUL(cv2, Z_MAKE((double)0.17392742256872693e0, (double)0.e0));
		c = (double)0.16999052179242813e0*b;
		s1 = a + c;
		s2 = a - c;
		cvpom = Z_ADD(cknn(s1, ifun, ka, dl, an), cknn(s2, ifun, ka, dl, an));
		//y = b * (y + .32607257743127307d0*(f(a + c) + f(a - c)))
		cvpom = Z_MUL(cvpom, Z_MAKE((double)0.32607257743127307e0, (double)0.e0));
		cvpom = Z_ADD(cv2, cvpom);
		cv2 = Z_MUL(cvpom, Z_MAKE(b, (double)0.e0));		

		if ((ifun == 1) || (ifun == 2)) res = (double)0.5e0*dl*((double)1.e0 + log((double)16.e0*an / dl)) / (pi*an);
		if (ifun == 3) res = dl * ((double)1.e0 + log((double)16.e0*an / dl)) / (pi*an);
		//cv = cv1 + cv2 + Z_MAKE(res, (double)0.e0);
		cv = Z_ADD(cv1, cv2);
		cv = Z_ADD(cv, Z_MAKE(res, (double)0.e0));
	}
	else
	{
		// ca³kowanie - kwadratura dwupunktowa   //dcqg(0.0, dl, ckmn, cv, 2);
		xl = (double)0.e0;
		xu = dl;
		a = (double)0.5e0*(xl + xu);
		b = xu - xl;
		c = (double)0.288675134594812882e0*b;
		s1 = a + c;
		s2 = a - c;
		cv = Z_ADD(ckmn(s1, ifun, ka, rmcx, rmcy, rmcz, rnx, rny, rnz, dcx, dcy, dcz, an, dl), ckmn(s2, ifun, ka, rmcx, rmcy, rmcz, rnx, rny, rnz, dcx, dcy, dcz, an, dl));
		cv = Z_MUL(cv, Z_MAKE(b * (double)0.5e0, (double)0.e0));
	}
	*cval = cv;
}


//----------------------------------------------------------------------
clDoubleComplex ckmn(double s, int ifun, double ka, double rmcx, double rmcy, double rmcz,
	double rnx, double rny, double rnz, double dcx, double dcy, double dcz, double an, double dl)
{
	//Jadro zredukowane(ca³ki potencja³owe dla segmentów)
	double r, drx, dry, drz;
	clDoubleComplex cj = Z_MAKE(0.0, 1.0);
	clDoubleComplex cval;

	drx = rmcx - rnx - s * dcx;
	dry = rmcy - rny - s * dcy;
	drz = rmcz - rnz - s * dcz;
	r = sqrt(an*an + drx * drx + dry * dry + drz * drz);
	if (ifun == 1)
	{
		cval = Z_MUL(Z_MAKE(s / dl, (double)0.e0), Z_MAKE(cos(ka*r) / r, -sin(ka*r) / r));      //pot.wektorowy Smn + 
	}
	else if (ifun == 2)
	{
		cval = Z_MUL(Z_MAKE((dl - s) / dl, (double)0.e0), Z_MAKE(cos(ka*r) / r, -sin(ka*r) / r));  //pot.wektorowy Smn -
	}
	else if (ifun == 3)
	{
		cval = Z_MAKE(cos(ka*r)/r,-sin(ka*r)/r);                      //pot.skalarny Smn + i Smn -
	}
	return cval;
}

//----------------------------------------------------------------------
clDoubleComplex cknn(double s, int ifun, double ka, double dl, double an)
{
	double bet, res, bnd2, r, rr;
	double pi = acos((double)-1.0);
	clDoubleComplex cj = Z_MAKE(0.0, 1.0);
	clDoubleComplex cval, cbnd1;
	
	r = 0.5*dl - s;
	rr = sqrt(an*an + r * r);	
	cbnd1 = Z_MAKE(cos(ka*rr), -sin(ka*rr));
	cbnd1 = Z_SUB(cbnd1, Z_ONE);
	cbnd1 = Z_DIV(cbnd1, Z_MAKE(rr, 0.));	
	bet = (an + an) / sqrt((double)4.0*an*an + r * r);	
	delick(bet, &res);	
	bnd2 = (bet*res + log(abs(r / ((double)8.0*an)))) / (pi*an);	
	if (ifun == 1)
	{
		//cknn = dcmplx(s / dl, 0.d0)*(cbnd1 + DCMPLX(bnd2, 0.d0)) //1 segment n - tej f.bazowej
		cval = Z_MAKE(s / dl, (double)0.e0);
		cval = Z_MUL(cval, Z_ADD(cbnd1, Z_MAKE(bnd2, (double)0.e0)));		
	}
	if (ifun == 2)
	{
		//cknn = dcmplx((dl - s) / dl, 0.d0)*(cbnd1 + DCMPLX(bnd2, 0.d0))!2 segment n - tej f.bazowej
		cval = Z_MAKE((dl - s) / dl, (double)0.e0);
		cval = Z_MUL(cval, Z_ADD(cbnd1, Z_MAKE(bnd2, (double)0.e0)));		
	}
	if (ifun == 3)
	{
		//cknn = cbnd1 + DCMPLX(bnd2, 0.d0)!pot.scalarny
		cval = Z_ADD(cbnd1, Z_MAKE(bnd2, (double)0.e0));		
	}		
	return cval;
}

//----------------------------------------------------------------------
void delick(double bet, double *res)
{
	//Complete elliptic integral of the first kind
	double a0, a1, a2, a3, a4;
	double b0, b1, b2, b3, b4;
	double am1, a, b, am12, am13, am14;

	a0 = (double)1.38629436112e0;
	a1 = (double)0.09666344259e0;
	a2 = (double)0.03590092383e0;
	a3 = (double)0.03742563713e0;
	a4 = (double)0.01451196212e0;
	b0 = (double)0.5e0;
	b1 = (double)0.12498593597e0;
	b2 = (double)0.06880248576e0;
	b3 = (double)0.03328355346e0;
	b4 = (double)0.00441787012e0;
	am1 = (double)1.e0 - bet * bet;
	a = a0 + a1 * am1;
	b = b0 + b1 * am1;
	if (am1 >= (double)1.0e-18)
	{
		am12 = am1 * am1;
		a = a + a2 * am12;
		b = b + b2 * am12;
	}
	if (am1 >= (double)1.0e-12)
	{
		am13 = am12 * am1;
		a = a + a3 * am13;
		b = b + b3 * am13;
	}
	if (am1 >= (double)1.e-9)
	{
		am14 = am13 * am1;
		a = a + a4 * am14;
		b = b + b4 * am14;
	}
	*res = a - b * log(am1);
}


