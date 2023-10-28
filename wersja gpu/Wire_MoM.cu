#include "Wire_MoM.h"

//using namespace std;

double cspeed = 2.99792448E+8;;
double pi = std::acos(-1.);

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
		std::cout << std::endl << " Process - " << process << " Timer took " << std::setprecision(4) << s << "s for " << info << std::endl;
	}
};

int main(int argc, char** argv)
{
	Timer timer;
	int rank;
	int num_procs;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	timer.process = rank;
	timer.info = "Whole calculations";

	if (rank == BROADCASTER_RANK)
	{
		std::cout << std::endl;
		std::cout << " -------------------------" << std::endl;
		std::cout << " Implementacja MoM w C++" << std::endl;
		std::cout << " Wersja dla drucikow z PEC" << std::endl;
		std::cout << " Analiza instalacji odgromowych" << std::endl;
		std::cout << " -------------------------" << std::endl;
		std::cout << std::endl;
	}

	//inicjalizacja tablic
	const int nnmax = 3000, ntmax = 5000;
	const int nwmax = 10000, npwmax = 20000, nsmax = nwmax + npwmax;
	int nnod = 0, nw = 0, nt = 0, npnt = 0;
	float *xnod = 0, *ynod = 0, *znod = 0; //wrad;
	int *iwn, *nsw, *iwt;
	float *wrad;
	int nseg;

	//zmienne
	bool wyb;
	int ienv, iexct, iload, ngap, ispc;
	float ampl, phas;


	xnod = new float[nnmax];
	ynod = new float[nnmax];
	znod = new float[nnmax];
	iwn = new int[2 * nwmax];
	wrad = new float[nwmax];
	nsw = new int[nwmax];
	iwt = new int[3 * ntmax];

	//zerowanie tablic
	if (rank == BROADCASTER_RANK)
	{
		for (int i = 0; i < nnmax; i++)
		{
			xnod[i] = 0.;
			ynod[i] = 0.;
			znod[i] = 0.;
		}
		for (int i = 0; i < nwmax; i++)
		{
			iwn[i] = 0;
			iwn[i + nwmax] = 0;
			wrad[i] = 0.;
			nsw[i] = 0;
		}
		for (int i = 0; i < ntmax; i++)
		{
			iwt[i] = 0;
			iwt[i + ntmax] = 0;
			iwt[i + 2 * ntmax] = 0;
		}
	}
	// wczytanie pliku geo
	if(rank == BROADCASTER_RANK)
		geoinp(nnmax, nwmax, ntmax, nnod, xnod, ynod, znod, nw, iwn, wrad, nsw, nt, iwt);

	wyb = true;
	// Wersja umo�liwiaj�ca uwzgl�dnienie PEC dla z=0 (ienv=0)
	//ienv = 2;
	ienv = 1;
	if (rank == BROADCASTER_RANK)
	{
		while (wyb)
		{
			std::cout << " Environment: free space       --> 1" << std::endl;
			std::cout << "              pec ground (z=0) --> 2" << std::endl;
			std::cout << " Select indicator (1 or 2):  " << ienv << std::endl;
			//cin >> ienv;
			if (ienv == 1 || ienv == 2) wyb = false;
		}
	}
	int isym;
	float sgnenv, sgnx, sgny, sgnz;
	if (ienv == 1)
	{
		sgnenv = 1.;
		isym = 0;
		sgnx = 1.;
		sgny = 1.;
		sgnz = 1.;
	}
	else
	{
		sgnenv = -1.;
		isym = 3;
		sgnx = 1.;
		sgny = 1.;
		sgnz = -1.;
	}

	//w�z�y
	float *x, *y, *z;
	// wi�ksza liczba punkt�w ni� nnod poniewa� w procedurze config tworzone s� nadmiarowe w�z�y, kt�re potem s� eliminowane
	x = new float[nsmax];
	y = new float[nsmax];
	z = new float[nsmax];

	//tablice zwi�zane z drucikami
	int *ise, *isp, *ips, *ipc, *iwe;
	float *rad;
	ise = new int[2 * nsmax];
	isp = new int[2 * nsmax];
	ips = new int[2 * nsmax];
	ipc = new int[nsmax];
	iwe = new int[2 * nsmax];
	rad = new float[nsmax];

	int npls = 0, nplsw = 0;

	if (rank == BROADCASTER_RANK)
	{
		//zerowanie tablic
		for (int i = 0; i < nsmax; i++)
		{
			x[i] = 0.;
			y[i] = 0.;
			z[i] = 0.;
			ise[i] = 0;
			ise[i + nsmax] = 0;
			isp[i] = 0;
			isp[i + nsmax] = 0;
			iwe[i] = 0;
			iwe[i + nsmax] = 0;
			rad[i] = 0.;
			ips[i] = 0;
			ips[i + nsmax] = 0;
			ipc[i] = 0;
		}
	}

	// konfiguracja modelu - wype�nienie wszystkich tablic zwi�zanych z przewodami
	if(rank == BROADCASTER_RANK)
		config(npnt, nnod, xnod, ynod, znod, x, y, z, nw, iwn, wrad, nsw, ise, isp, ips, ipc, iwe, rad, npls, nplsw, nseg, nsmax, nwmax, ienv);

	// Pobudzenie -> W tej chwili domy�lnie ustawione na 1 -> generator
	wyb = true;
	iexct = 1;
	if(rank == BROADCASTER_RANK)
	{
		while (wyb)
		{
			std::cout << " -------------------" << std::endl;
			std::cout << " EXCITATION: " << std::endl;
			std::cout << " 1 - delta gap source" << std::endl;
			std::cout << " 2 - incident wave" << std::endl;
			std::cout << " Select 1 or 2: " << iexct << std::endl;
			//cin >> iexct;
			//std::cout << std::endl;
			if (iexct == 1 || iexct == 2) wyb = false;
		}
	}
	//if (iexct == 2) std::cout << " ERROR - ten typ pobudzenia nie zaimplementowany!" << std::endl;

	//definiowanie pobudzenia - uproszczone - zak�adamy tylko jeden generator na przewodach -> do analizy LPS
	if (rank == BROADCASTER_RANK)
	{
		if (iexct == 1) // pobudzenie w postaci generatora nexct=1, nexct - liczba generator�w
		{
			wyb = true;
			while (wyb)
			{
				//ngap = 592;
				ngap = 100;
				std::cout << " ------------------------------" << std::endl;
				std::cout << " Source (source zone on wire): " << ngap << std::endl;
				//cin >> ngap; // numer funkcji bazowej na przewodach definiuj�cy miejsce generatora
				if (ngap > 0 && ngap <= npls)
				{
					wyb = false;
				}
				else
				{
					std::cout << "    ERROR (wrong number of source zone)" << std::endl;
				}
			}
			ampl = 1.0;
			std::cout << "       Specify source attributes" << std::endl;
			std::cout << "                 magnitude [V]: " << ampl << std::endl;
			//cin >> ampl;
			ampl = abs(ampl);
			phas = 0.0;
			std::cout << "               phase [degrees]: " << phas << std::endl;
			//cin >> phas;				
		}
	}
	// Obci��enia
	wyb = true;
	//iload = 1;
	iload = 2;
	if (rank == BROADCASTER_RANK)
	{
		while (wyb)
		{
			std::cout << " -----------------------------------" << std::endl;
			std::cout << " Lumped loadings (1-> yes, 2-> no): " << iload << std::endl;
			//cin >> iload;
			if (iload == 1 || iload == 2)
			{
				wyb = false;
			}
			else
			{
				std::cout << "     ERROR!" << std::endl;
			}
		}
	}
	//skalowanie widma -> W tej chwili domy�lnie ustawione na NIE -> 2
	ispc = 2;
	wyb = true;
	if (rank == BROADCASTER_RANK)
	{
		while (wyb)
		{
			std::cout << " -----------------------------------" << std::endl;
			std::cout << " spectrum weight (1-> yes, 2-> no): " << ispc << std::endl;
			//cin >> ispc;
			if (ispc == 1 || ispc == 2)
			{
				wyb = false;
			}
			else
			{
				std::cout << "    ERROR!" << std::endl;
			}
		}
	}
	//wczytanie cz�stotliwo�ci z pliku
	float *freqtab;
	int nfreq;
	if (rank == BROADCASTER_RANK)
	{
		// plik z cz�stotliwo�ciami dla kt�rych wykonywane obliczenia -> pierwszy wiersz zawiera liczb� punkt�w
		fstream freqfile;
		freqfile.open("freq.txt", ios::in);
		freqfile >> nfreq;
		freqtab = new float[nfreq];

		std::cout << std::endl << std::endl;
		std::cout << "  -------------------------------------" << std::endl;
		std::cout << "  Number of frequency points: " << nfreq << std::endl;
		std::cout << "  -------------------------------------" << std::endl;

		for (int i = 0; i < nfreq; i++)
		{
			freqfile >> freqtab[i];
		}

		//resfile.close();
		freqfile.close();
	}
	MPI_Bcast(&npls, 1, MPI_INT, BROADCASTER_RANK, MPI_COMM_WORLD);
	MPI_Bcast(&nfreq, 1, MPI_INT, BROADCASTER_RANK, MPI_COMM_WORLD);

	if (rank != BROADCASTER_RANK) //only BROADCASTER_RANK allocated memory for freqtab
		freqtab = new float[nfreq];

	MPI_Bcast(freqtab, nfreq, MPI_FLOAT, BROADCASTER_RANK, MPI_COMM_WORLD);
	MPI_Bcast(x, nsmax, MPI_FLOAT, BROADCASTER_RANK, MPI_COMM_WORLD);
	MPI_Bcast(y, nsmax, MPI_FLOAT, BROADCASTER_RANK, MPI_COMM_WORLD);
	MPI_Bcast(z, nsmax, MPI_FLOAT, BROADCASTER_RANK, MPI_COMM_WORLD);
	MPI_Bcast(rad, nsmax, MPI_FLOAT, BROADCASTER_RANK, MPI_COMM_WORLD);
	MPI_Bcast(ipc, nsmax, MPI_INT, BROADCASTER_RANK, MPI_COMM_WORLD);
	MPI_Bcast(ips, 2 * nsmax, MPI_INT, BROADCASTER_RANK, MPI_COMM_WORLD);
	MPI_Bcast(ise, 2 * nsmax, MPI_INT, BROADCASTER_RANK, MPI_COMM_WORLD);
	MPI_Bcast((void*)&ienv, 1, MPI_INT, BROADCASTER_RANK, MPI_COMM_WORLD);
	MPI_Bcast((void*)&sgnx, 1, MPI_FLOAT, BROADCASTER_RANK, MPI_COMM_WORLD);
	MPI_Bcast((void*)&sgny, 1, MPI_FLOAT, BROADCASTER_RANK, MPI_COMM_WORLD);
	MPI_Bcast((void*)&sgnz, 1, MPI_FLOAT, BROADCASTER_RANK, MPI_COMM_WORLD);
	MPI_Bcast((void*)&sgnenv, 1, MPI_FLOAT, BROADCASTER_RANK, MPI_COMM_WORLD);
	MPI_Bcast((void*)&iload, 1, MPI_INT, BROADCASTER_RANK, MPI_COMM_WORLD);
	MPI_Bcast((void*)&phas, 1, MPI_FLOAT, BROADCASTER_RANK, MPI_COMM_WORLD);
	MPI_Bcast((void*)&ngap, 1, MPI_INT, BROADCASTER_RANK, MPI_COMM_WORLD);
	MPI_Bcast((void*)&ampl, 1, MPI_FLOAT, BROADCASTER_RANK, MPI_COMM_WORLD);

	//alokacja macierzy z
	clDoubleComplex  *cz;
	cz = new clDoubleComplex[npls*npls];

	cuDoubleComplex* cu_cu;
	cuDoubleComplex* ci_cu;
	cu_cu = new clDoubleComplex[npls];
	ci_cu = new clDoubleComplex[npls];

	// zaczynamy p�tl� po cz�stotliwo�ci
	
	for (int i = rank; i < nfreq; i = i+num_procs)
	{
		CalculateByFreq
		(
			i, npls, nfreq, freqtab[i], x, y, z, rad, ipc, ips, ise, nsmax, ienv, 
			sgnx, sgny, sgnz,sgnenv, iload, phas, ngap, ampl, cz, cu_cu, ci_cu, rank
		);
		
	}

	delete[] xnod, ynod, znod;
	delete[] iwn, wrad, nsw, iwt;
	delete[] x, y, z;
	delete[] ise, isp, ips, ipc, iwe, rad;
	delete[] cz;
	delete[] freqtab;

	//DS - remove if changed to unique_ptr
	delete[] cu_cu, ci_cu;
	MPI_Finalize();
	return (0);

}
//---------------------------------------------------------------------------------------------------------------------------
void CalculateByFreq(int i, int npls, int nfreq, float freq,  
	float *x, float* y, float *z, float* rad, int* ipc, int* ips, int* ise,  const int nsmax, 
	int ienv, float sgnx, float sgny, float sgnz, float sgnenv, int iload, 
	float phas, int ngap, float ampl, cuDoubleComplex* cz, cuDoubleComplex* cu_cu, cuDoubleComplex* ci_cu, int rank)
{
	double ambda, ka, arg;
	clDoubleComplex cload;

	ambda = cspeed / freq * 1E-6;
	ka = (pi + pi) / ambda;
	
	// zmatrix
	zmatrix_gpu(cz, x, y, z, rad, ipc, ips, ise, nsmax, ka, npls, ienv, sgnx, sgny, sgnz, sgnenv, rank);

	// addloads
	// w tej chwili bez osobnej procedury
	if (iload == 1)
	{
		// dodanie obi��enia do kana�u
		cload = Z_MAKE(1.e0, 2. * pi * freq * 4.5);		// 4.5 uH ale freq w MHz	!e - 6
		int nl;
		for (nl = 592 - 1; nl <= 2581 - 1; nl++)  //anis1m_loop.geo
		{
			cz[nl + nl * npls] = Z_ADD(cz[nl + nl * npls], cload);
		}
		//dodanie obci��enia do p�tli
		nl = 2613 - 1;	//anis1m_loop.geo
		cload = Z_MAKE(30., 0.);
		cz[nl + nl * npls] = Z_ADD(cz[nl + nl * npls], cload);
		//dodanie obci��enia do kr�tkiego dipola w �rodku p�tli -> do wyznaczania pola E
		nl = 2582 - 1;	//anis1m_loop.geo
		cload = Z_MAKE(1.e12, 0.);
		cz[nl + nl * npls] = Z_ADD(cz[nl + nl * npls], cload);
	}

	int batchSize = 1;
	int nrhs = 1;

	for (int i = 0; i < npls; i++)
	{
		cu_cu[i] = Z_ZERO;
	}
	arg = pi * phas / 180.;
	cu_cu[ngap - 1] = Z_MUL(Z_MAKE(ampl, 0.), Z_MAKE(cos(arg), sin(arg)));

	// Substitution and solution equation
	for (int i = 0; i < npls; i++) ci_cu[i] = cu_cu[i];

	LU_cuBLAS(npls, nrhs, ci_cu, cu_cu, cz, batchSize, rank);

	clDoubleComplex csample_cu = Z_DIV(cu_cu[ngap - 1], ci_cu[ngap - 1]);
	std::cout << std::endl << " Process - " << rank << " Freq = " << freq << "\t" << Z_REAL(csample_cu) << "\t" << Z_IMAG(csample_cu) << std::endl;
	//resfile << freq << "\t" << Z_REAL(csample_cu) << "\t" << Z_IMAG(csample_cu) << std::endl;

	// Pole E unormowane do pr�du u podstawy kana�u - anis1m_loop.geo
	/*clDoubleComplex csample = Z_DIV(ci[2582 - 1], ci[592 - 1]);
	csample = Z_MUL(csample, Z_MAKE((double)0.4e15, (double)0.e0));
	resfile << freq << "\t" << Z_REAL(csample) << "\t" << Z_IMAG(csample) << "\t" << Z_ABS(csample) << std::endl; */
}

//---------------------------------------------------------------------------------------------------------------------------

void LU_cuBLAS(int n, int nrhs, cuDoubleComplex* ci_cu, cuDoubleComplex* cu_cu, cuDoubleComplex* cz, int batchSize, int rank)
{
	//Timer timer;
	//timer.rank = rank;
	//timer.info = "LU_cuBLAS";

	int* pivotArray = new int;
	int* infoArray = new int;
	int lda = n;
	int infoStatus;
	
	float elapsedTime;
	float elapsedTimeZgetrf;
	float elapsedTimeZgetrs;

	cudaEvent_t start, stop;
	cudaEvent_t startZgetrf, stopZgetrf;
	cudaEvent_t startZgetrs, stopZgetrs;

	cudaStream_t stream;
	cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking);

	cublasHandle_t handle;
	cublasCreate_v2(&handle);
	cublasSetStream_v2(handle, stream);

	cuDoubleComplex* cz_cu = nullptr;
	cuDoubleComplex* results_cu = nullptr;
	cuDoubleComplex** cz_d;
	cuDoubleComplex** results_d;

	CUDA_CALL(cudaMalloc<cuDoubleComplex>(&cz_cu, sizeof(cuDoubleComplex) * n * n),"cudaMalloc Failed!\n");
	cudaMalloc<cuDoubleComplex>(&results_cu, sizeof(cuDoubleComplex) * n);
	cudaMalloc<cuDoubleComplex*>(&cz_d, sizeof(cuDoubleComplex*));
	cudaMalloc<cuDoubleComplex*>(&results_d, sizeof(cuDoubleComplex*));
	cudaMalloc<int>(&pivotArray, sizeof(int) * n * batchSize);
	cudaMalloc<int>(&infoArray, sizeof(int) * n);

	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventCreate(&startZgetrf);
	cudaEventCreate(&stopZgetrf);
	cudaEventCreate(&startZgetrs);
	cudaEventCreate(&stopZgetrs);

	cudaEventRecord(start, 0);
	cudaMemcpy(cz_cu, cz, n * n * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
	cudaMemcpy(results_cu, ci_cu, n * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
	cudaMemcpy(cz_d, &cz_cu, sizeof(cuDoubleComplex*), cudaMemcpyHostToDevice);

	cudaEventRecord(startZgetrf, 0);
	cublasZgetrfBatched(handle, n, cz_d, lda, pivotArray, infoArray, batchSize);
	cudaEventRecord(stopZgetrf, 0);
	cudaEventSynchronize(stopZgetrf);
	cudaEventElapsedTime(&elapsedTimeZgetrf, startZgetrf, stopZgetrf);
	std::cout << std::endl << " Process - " << rank << " cublasZgetrf execution time = " << std::setprecision(4) << elapsedTimeZgetrf / 1.0e3 << " s" << std::endl;

	cudaMemcpy(results_d, &results_cu, sizeof(cuDoubleComplex*), cudaMemcpyHostToDevice);

	cudaEventRecord(startZgetrs, 0);
	cublasZgetrsBatched(handle, CUBLAS_OP_N, n, nrhs, cz_d, lda, pivotArray, results_d, lda, &infoStatus, batchSize);
	cudaEventRecord(stopZgetrs, 0);
	cudaEventSynchronize(stopZgetrs);
	cudaEventElapsedTime(&elapsedTimeZgetrs, startZgetrs, stopZgetrs);
	std::cout << std::endl << " Process - " << rank << " cublasZgetrs execution time = " << std::setprecision(4) << elapsedTimeZgetrs / 1.0e3 << " s" << std::endl;
	
	cudaFree(pivotArray), cudaFree(infoArray), cublasDestroy_v2(handle);

	cudaMemcpy(ci_cu, results_cu, sizeof(clDoubleComplex) * n, cudaMemcpyDeviceToHost);

	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&elapsedTime, start, stop);

	std::cout << std::endl << " Process - " << rank << " LU_cuBLAS execution time including data transfer = " << std::setprecision(4) << elapsedTime / 1.0e3 << " s" << std::endl;
	std::cout << std::endl << " Process - " << rank << " LU_cuBLAS Data transfer time = " << std::setprecision(4) << (elapsedTime - elapsedTimeZgetrf - elapsedTimeZgetrs) / 1.0e3 << " s" << std::endl;
	
	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	cudaEventDestroy(startZgetrf);
	cudaEventDestroy(stopZgetrf);
	cudaEventDestroy(startZgetrs);
	cudaEventDestroy(stopZgetrs);

	cudaFree(cz_cu), cudaFree(results_cu);
	cudaFree(cz_d), cudaFree(results_d);
} 

//----------------------------------------------------------------------

__device__ double delick_gpu(double bet)
{
	double a0 = (double)1.38629436112e0;
	double a1 = (double)0.09666344259e0;
	double a2 = (double)0.03590092383e0;
	double a3 = (double)0.03742563713e0;
	double a4 = (double)0.01451196212e0;
	double b0 = (double)0.5e0;
	double b1 = (double)0.12498593597e0;
	double b2 = (double)0.06880248576e0;
	double b3 = (double)0.03328355346e0;
	double b4 = (double)0.00441787012e0;
	double res, a, b;
	double am1 = (double)1.0e0 - bet*bet;
	double am12, am13, am14;

	a = a0 + a1*am1;
	b = b0 + b1*am1;

	if (am1 >= (double)1.0e-18){
		am12 = am1*am1;
		a = a + a2*am12;
		b = b + b2*am12;
		if (am1 >= (double)1.0e-12){
			am13 = am12*am1;
			a = a + a3*am13;
			b = b + b3*am13;
			if (am1 >= (double)1.0e-9){
				am14 = am13*am1;
				a = a + a4*am14;
				b = b + b4*am14;
			}
		}
	}

	res = a - b*log(am1);
	return res;

}

//------------------------------------------------

__device__ bool war_gpu_pec_re(int ise_ms1, int ise_ms2, int ise_ns1, int ise_ns2)
{
	if ((ise_ns1 == ise_ms1  && ise_ns2 == ise_ms2) || (ise_ns1 == ise_ms2  && ise_ns2 == ise_ms1))
	{
		return true;
	}
	else {
		return false;
	}
}
//------------------------------------------------

__device__ bool war_gpu_pec_im(int ise_ms1, int ise_ms2, int ise_ns1, int ise_ns2, float z_ise_ms1, float z_ise_ms2)
{
	if (((ise_ns1 == ise_ms1 && ise_ns2 == ise_ms2) && (z_ise_ms1 == 0.0f && z_ise_ms2 == 0.0f))
		|| ((ise_ns1 == ise_ms2 && ise_ns2 == ise_ms1) && (z_ise_ms1 == 0.0f && z_ise_ms2 == 0.0f)))
	{
		return true;
	}
	else {
		return false;
	}
}
//----------------------------------------------------------------------
__global__ void wim_gpu_cknn_cwire1_re1(cuDoubleComplex *cz, float *x, float *y, float *z,
	float *rad, int *ipc, int *ips, int *ise, int nsmax,
	int n, double ka, double eta0,
	int realSubSliceSize, int k, int kk)
{

	double pi = acos(double(-1.e0));
	cuDoubleComplex cvval1, csval1;

	// Block index
	int bx = blockIdx.x;
	int by = blockIdx.y;

	// Thread index
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	// Accumulate row i of A and column j of B
	int i = tx + bx * blockDim.x;
	int j = ty + by * blockDim.y;

	//c16vec c1 = make_c16vec(Z_ONE,Z_MAKE(0.0f, 1.0f),Z_MUL(Z_MAKE(2.0f, 0.0f),Z_ONE));

	if (i < n  && j < realSubSliceSize) {
		int mp = abs(ipc[i]);
		int ms = ips[i];
		int mlr = ise[ms - 1];
		//double am = rad[ms-1];
		double3 rmp = make_double3((double)x[mp - 1], (double)y[mp - 1], (double)z[mp - 1]);
		double3 rmlr = make_double3((double)x[mlr - 1], (double)y[mlr - 1], (double)z[mlr - 1]);
		double3 rmc = make_double3((double)0.5e0*(rmlr.x + rmp.x), (double)0.5e0*(rmlr.y + rmp.y), (double)0.5e0*(rmlr.z + rmp.z));
		float znak = -1.0f;
		double3 lmc = make_double3((double)znak*(rmc.x - rmp.x), (double)znak*(rmc.y - rmp.y), (double)znak*(rmc.z - rmp.z));

		// ------------ Zamiast prcodeury cwire1 -------------------
		//ifun=1 i ifun=2=> res_cv
		//fun=3 => res_cs		  
		//call cwire1_gpu(ms,n,cvval1)             ! obliczanie calki na Sn+  pot.wek.
		//call cwire1(ms,n,csval1)

		//int np = abs(ipc[k-1+kk-1+j]);
		//int ns = ips[k-1+kk-1+j];
		int np = abs(ipc[k + kk + j]);
		int ns = ips[k + kk + j];
		double3 rnp = make_double3((double)x[np - 1], (double)y[np - 1], (double)z[np - 1]);
		int nlr = ise[ns - 1];
		double an = rad[ns - 1];
		double3 rn = make_double3((double)x[nlr - 1], (double)y[nlr - 1], (double)z[nlr - 1]);
		double3 rd = make_double3(rnp.x - rn.x, rnp.y - rn.y, rnp.z - rn.z);
		double dl = sqrt(rd.x*rd.x + rd.y*rd.y + rd.z*rd.z);
		double dlh = (double)0.5e0*dl;
		double3 dc = make_double3(rd.x / dl, rd.y / dl, rd.z / dl);

		int ise_ms1 = ise[ms - 1];
		int ise_ms2 = ise[ms - 1 + nsmax];
		int ise_ns1 = ise[ns - 1];
		int ise_ns2 = ise[ns - 1 + nsmax];


		if (war_gpu_pec_re(ise_ms1, ise_ms2, ise_ns1, ise_ns2)) {
			//call dcqg(0.d0,dlh,cknn,cv1,4)
			//nord=4
			//  res_cv = 0.5d0*dl*( 1.d0+DLOG(16.d0*an/dl) )/(pi*an)	! dla ifun=1
			//  res_cs = dl*( 1.d0+DLOG(16.d0*an/dl) )/(pi*an)	        ! dla ifun=3		    
			// --------------- cv1 ------------------
			//dcqg_a = 0.5d0*(0.0d0+dlh)
			//dcqg_b = dlh-0.0d0
			double dcqg_a = (double)0.5e0*dlh;
			double dcqg_b = dlh;
			//
			double dcqg_c = (double)0.43056815579702629e0 * dcqg_b;
			//y=.17392742256872693d0*(f(a+c)+f(a-c))
			//s=a+c
			double cknn_s = dcqg_a + dcqg_c;
			// ---- cknn ----
			double cknn_r = (double)0.5e0*dl - cknn_s;
			double cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			///cknn_cbnd1 = ( dcos(ka*cknn_rr) - cj*dsin(ka*cknn_rr) - cone )/cknn_rr
			///cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cuDoubleComplex cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			double cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			double cknn_res = delick_gpu(cknn_bet);
			double cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))        !1 segment n-tej f.bazowej
			cuDoubleComplex cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//cv = Z_MUL( Z_MAKE(cknn_s/dl,(double)0.0e0),cv );
			cuDoubleComplex cvv1a = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)          ! pot. scalarny 
			//cv = Z_ADD( cknn_cbnd1,Z_MAKE(cknn_bnd2,(double)0.0e0) );
			cuDoubleComplex csv1a = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			// ---- cknn ----
			//s=a-c
			cknn_s = dcqg_a - dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			//cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//!cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))        !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			cv = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			cvv1a = Z_MUL(Z_MAKE((double)0.17392742256872693e0, (double)0.0e0), Z_ADD(cv, cvv1a));
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)            ! pot. scalarny 
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//csv1a = .17392742256872693d0*( cv + csv1a )
			csv1a = Z_MUL(Z_MAKE((double)0.17392742256872693e0, (double)0.0e0), Z_ADD(cv, csv1a));
			// ---- cknn ---- 

			dcqg_c = (double)0.16999052179242813e0 * dcqg_b;
			//y=b*(y+.32607257743127307d0*(f(a+c)+f(a-c))) 
			//s=a+c
			cknn_s = dcqg_a + dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			//cknn_cbnd1 = ( dcos(ka*cknn_rr) - cj*dsin(ka*cknn_rr) - cone )/cknn_rr
			//cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))         !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//cv = Z_MUL( Z_MAKE(cknn_s/dl,(double)0.0e0),cv );
			cuDoubleComplex cvv1b = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)               ! pot. scalarny 
			//cv = Z_ADD( cknn_cbnd1,Z_MAKE(cknn_bnd2,(double)0.0e0) );
			cuDoubleComplex csv1b = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			// ---- cknn ----
			//s=a-c
			cknn_s = dcqg_a - dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			//cknn_cbnd1 = ( dcos(ka*cknn_rr) - cj*dsin(ka*cknn_rr) - cone )/cknn_rr
			//cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))            !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			cv = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			cvv1b = Z_MUL(Z_MAKE((double)0.32607257743127307e0, (double)0.0e0), Z_ADD(cv, cvv1b));
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)                      ! pot. scalarny
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//csv1b = .32607257743127307d0*( cv + csv1b ) 
			csv1b = Z_MUL(Z_MAKE((double)0.32607257743127307e0, (double)0.0e0), Z_ADD(cv, csv1b));
			// ---- cknn ---- 
			//cv1v = dcqg_b*( cvv1a+cvv1b )
			//cv1s = dcqg_b*( csv1a+csv1b ) 
			cuDoubleComplex cv1v = Z_MUL(Z_MAKE(dcqg_b, (double)0.0e0), Z_ADD(cvv1a, cvv1b));
			cuDoubleComplex cv1s = Z_MUL(Z_MAKE(dcqg_b, (double)0.0e0), Z_ADD(csv1a, csv1b));
			// --------------- cv1 ------------------


			//call dcqg(dlh,dl,cknn,cv2,4)	      
			//nord=4
			// --------------- cv2 ------------------
			dcqg_a = (double)0.5e0*(dlh + dl);
			dcqg_b = dl - dlh;
			//
			dcqg_c = (double)0.43056815579702629e0 * dcqg_b;
			//y=.17392742256872693d0*(f(a+c)+f(a-c))
			//s=a+c
			cknn_s = dcqg_a + dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			///cknn_cbnd1 = ( dcos(ka*cknn_rr) - cj*dsin(ka*cknn_rr) - cone )/cknn_rr
			///cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))        !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			cv = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			cvv1a = cv;
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)          ! pot. scalarny 
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			csv1a = cv;
			// ---- cknn ----
			//s=a-c
			cknn_s = dcqg_a - dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			//cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//!cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))        !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			cv = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			cvv1a = Z_MUL(Z_MAKE((double)0.17392742256872693e0, (double)0.0e0), Z_ADD(cv, cvv1a));
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)            ! pot. scalarny 
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//csv1a = .17392742256872693d0*( cv + csv1a )
			csv1a = Z_MUL(Z_MAKE((double)0.17392742256872693e0, (double)0.0e0), Z_ADD(cv, csv1a));
			// ---- cknn ---- 

			dcqg_c = (double)0.16999052179242813e0 * dcqg_b;
			//y=b*(y+.32607257743127307d0*(f(a+c)+f(a-c))) 
			//s=a+c
			cknn_s = dcqg_a + dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			//cknn_cbnd1 = ( dcos(ka*cknn_rr) - cj*dsin(ka*cknn_rr) - cone )/cknn_rr
			//cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))         !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//cv = Z_MUL( Z_MAKE(cknn_s/dl,(double)0.0e0),cv );
			cvv1b = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)               ! pot. scalarny 
			//cv = Z_ADD( cknn_cbnd1,Z_MAKE(cknn_bnd2,(double)0.0e0) );
			csv1b = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			// ---- cknn ----
			//s=a-c
			cknn_s = dcqg_a - dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			//cknn_cbnd1 = ( dcos(ka*cknn_rr) - cj*dsin(ka*cknn_rr) - cone )/cknn_rr
			//cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))            !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			cv = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			cvv1b = Z_MUL(Z_MAKE((double)0.32607257743127307e0, (double)0.0e0), Z_ADD(cv, cvv1b));
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)                      ! pot. scalarny
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//csv1b = .32607257743127307d0*( cv + csv1b ) 
			csv1b = Z_MUL(Z_MAKE((double)0.32607257743127307e0, (double)0.0e0), Z_ADD(cv, csv1b));
			// ---- cknn ---- 
			//cv2v = dcqg_b*( cvv1a+cvv1b )
			//cv2s = dcqg_b*( csv1a+csv1b  
			cuDoubleComplex cv2v = Z_MUL(Z_MAKE(dcqg_b, (double)0.0e0), Z_ADD(cvv1a, cvv1b));
			cuDoubleComplex cv2s = Z_MUL(Z_MAKE(dcqg_b, (double)0.0e0), Z_ADD(csv1a, csv1b));
			// --------------- cv2 ------------------

			//res_cv = 0.5d0*dl*( 1.d0+DLOG(16.d0*an/dl) )/(pi*an)	! dla ifun=1 i ifun=2
			//res_cs = dl*( 1.d0+DLOG(16.d0*an/dl) )/(pi*an)	        ! dla ifun=3
			double res_cs = dl*((double)1.0e0 + log((double)16.0e0*an / dl)) / (pi*an);
			double res_cv = (double)0.5e0*res_cs;
			//cvval1 = cv1v + cv2v + DCMPLX(res_cv,0.d0)     
			//csval1 = cv1s + cv2s + DCMPLX(res_cs,0.d0)
			cvval1 = Z_ADD(cv2v, Z_MAKE(res_cv, (double)0.0e0));
			cvval1 = Z_ADD(cv1v, cvval1);
			csval1 = Z_ADD(cv2s, Z_MAKE(res_cs, (double)0.0e0));
			csval1 = Z_ADD(cv1s, csval1);
		}
		else {
			//call dcqg(0.d0,dl,ckmn,cv,2)
			//cvval1 = dcmplx(0.0d0,0.0d0)
			//csval1 = dcmplx(0.0d0,0.0d0)
			cvval1 = Z_MAKE((double)0.0e0, (double)0.0e0);
			csval1 = Z_MAKE((double)0.0e0, (double)0.0e0);
		}
		// ------------ Zamiast prcodeury cwire1 -------------------

		//tvecl1 = dc%x*lmc%x+dc%y*lmc%y+dc%z*lmc%z       ! obliczane !lm+ * 1n+ dla k=1, obliczane !lm- * 1n+ dla k=1
		double tvecl1 = dc.x*lmc.x + dc.y*lmc.y + dc.z*lmc.z;
		//csval1 = csval1/dl
		csval1 = Z_DIV(csval1, Z_MAKE(dl, (double)0.0e0));

		//cvval = dcmplx(tvecl1,0.d0)*cvval1
		cuDoubleComplex cvval = Z_MUL(cvval1, Z_MAKE(tvecl1, (double)0.0e0));
		//csval = csval1
		cuDoubleComplex csval = csval1;

		//cz_ij = cz_ij+ka*ka*cvval + znak * csval
		cuDoubleComplex cz_ij = Z_MUL(Z_MAKE((double)znak, (double)0.0e0), csval);
		cz_ij = Z_ADD(Z_MUL(Z_MAKE(ka*ka, (double)0.0e0), cvval), cz_ij);

		//cz_ij = cz_ij*cj*eta0/(4.0e0*pi*ka)
		cz_ij = Z_MUL(cz_ij, Z_MAKE(eta0 / ((double)4.0e0*pi*ka), (double)0.0e0));
		cz_ij = Z_MUL(cz_ij, Z_ONE_IMG);


		//cz_d[i_m+j_n*n] = Z_MAKE(10.0f,0.0f);
		//cz_d[i_m+j_n*n] = Z_MUL ( cz_d[i_m+j_n*n] , (Z_MAKE(10.0f,0.0f));
		//cz_d[i_m+j_n*n] = Z_MAKE((double)i_m,(double)j_n);
		//cz_d[i_m+j_n*n]=make_cuDoubleComplex( (i_m),(j_n) );
		//cz_d[ i+j*n ] =  Z_MAKE(float(i*n+j),0.0f) ;
		//cz[ i+j*n ] = Z_MAKE(x[mp-1],z[mp-1]);
		//cz[ i+j*n ] = Z_MAKE(dc.x,dc.z);
		//cz[ i+j*n ] = cknn_cbnd1;
		//cz[ i+j*n ] = Z_MAKE(cknn_res,cknn_bnd2);
		cz[i + (j + kk)*n] = cz_ij;
		//cz[ i+j*n ] =  Z_ADD(c1.z,Z_MAKE(float(i*n+j),0.0f)) ;
	}

}
//----------------------------------------------------------------------
__global__ void wim_gpu_cknn_cwire1_re2(cuDoubleComplex *cz, float *x, float *y, float *z,
	float *rad, int *ipc, int *ips, int *ise, int nsmax,
	int n, double ka, double eta0,
	int realSubSliceSize, int k, int kk)
{

	double pi = acos(double(-1.e0));
	cuDoubleComplex cvval1, csval1;

	// Block index
	int bx = blockIdx.x;
	int by = blockIdx.y;

	// Thread index
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	// Accumulate row i of A and column j of B
	int i = tx + bx * blockDim.x;
	int j = ty + by * blockDim.y;

	//c16vec c1 = make_c16vec(Z_ONE,Z_MAKE(0.0f, 1.0f),Z_MUL(Z_MAKE(2.0f, 0.0f),Z_ONE));

	if (i < n  && j < realSubSliceSize) {
		int mp = abs(ipc[i]);
		//int ms = ips[i];
		//int mlr = ise[ms-1];
		int ms = ips[i + nsmax];
		int mlr = ise[ms - 1 + nsmax];
		//double am = rad[ms-1];
		double3 rmp = make_double3((double)x[mp - 1], (double)y[mp - 1], (double)z[mp - 1]);
		double3 rmlr = make_double3((double)x[mlr - 1], (double)y[mlr - 1], (double)z[mlr - 1]);
		double3 rmc = make_double3((double)0.5e0*(rmlr.x + rmp.x), (double)0.5e0*(rmlr.y + rmp.y), (double)0.5e0*(rmlr.z + rmp.z));
		float znak = 1.0f;
		double3 lmc = make_double3((double)znak*(rmc.x - rmp.x), (double)znak*(rmc.y - rmp.y), (double)znak*(rmc.z - rmp.z));

		// ------------ Zamiast prcodeury cwire1 -------------------
		//ifun=1 i ifun=2=> res_cv
		//fun=3 => res_cs		  
		//call cwire1_gpu(ms,n,cvval1)             ! obliczanie calki na Sn+  pot.wek.
		//call cwire1(ms,n,csval1)

		//int np = abs(ipc[k-1+kk-1+j]);
		//int ns = ips[k-1+kk-1+j];
		int np = abs(ipc[k + kk + j]);
		int ns = ips[k + kk + j];
		double3 rnp = make_double3((double)x[np - 1], (double)y[np - 1], (double)z[np - 1]);
		int nlr = ise[ns - 1];
		double an = rad[ns - 1];
		double3 rn = make_double3((double)x[nlr - 1], (double)y[nlr - 1], (double)z[nlr - 1]);
		double3 rd = make_double3(rnp.x - rn.x, rnp.y - rn.y, rnp.z - rn.z);
		double dl = sqrt(rd.x*rd.x + rd.y*rd.y + rd.z*rd.z);
		double dlh = (double)0.5e0*dl;
		double3 dc = make_double3(rd.x / dl, rd.y / dl, rd.z / dl);

		int ise_ms1 = ise[ms - 1];
		int ise_ms2 = ise[ms - 1 + nsmax];
		int ise_ns1 = ise[ns - 1];
		int ise_ns2 = ise[ns - 1 + nsmax];


		if (war_gpu_pec_re(ise_ms1, ise_ms2, ise_ns1, ise_ns2)) {
			//call dcqg(0.d0,dlh,cknn,cv1,4)
			//nord=4
			//  res_cv = 0.5d0*dl*( 1.d0+DLOG(16.d0*an/dl) )/(pi*an)	! dla ifun=1
			//  res_cs = dl*( 1.d0+DLOG(16.d0*an/dl) )/(pi*an)	        ! dla ifun=3		    
			// --------------- cv1 ------------------
			//dcqg_a = 0.5d0*(0.0d0+dlh)
			//dcqg_b = dlh-0.0d0
			double dcqg_a = (double)0.5e0*dlh;
			double dcqg_b = dlh;
			//
			double dcqg_c = (double)0.43056815579702629e0 * dcqg_b;
			//y=.17392742256872693d0*(f(a+c)+f(a-c))
			//s=a+c
			double cknn_s = dcqg_a + dcqg_c;
			// ---- cknn ----
			double cknn_r = (double)0.5e0*dl - cknn_s;
			double cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			///cknn_cbnd1 = ( dcos(ka*cknn_rr) - cj*dsin(ka*cknn_rr) - cone )/cknn_rr
			///cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cuDoubleComplex cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			double cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			double cknn_res = delick_gpu(cknn_bet);
			double cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))        !1 segment n-tej f.bazowej
			cuDoubleComplex cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//cv = Z_MUL( Z_MAKE(cknn_s/dl,(double)0.0e0),cv );
			cuDoubleComplex cvv1a = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)          ! pot. scalarny 
			//cv = Z_ADD( cknn_cbnd1,Z_MAKE(cknn_bnd2,(double)0.0e0) );
			cuDoubleComplex csv1a = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			// ---- cknn ----
			//s=a-c
			cknn_s = dcqg_a - dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			//cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//!cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))        !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			cv = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			cvv1a = Z_MUL(Z_MAKE((double)0.17392742256872693e0, (double)0.0e0), Z_ADD(cv, cvv1a));
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)            ! pot. scalarny 
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//csv1a = .17392742256872693d0*( cv + csv1a )
			csv1a = Z_MUL(Z_MAKE((double)0.17392742256872693e0, (double)0.0e0), Z_ADD(cv, csv1a));
			// ---- cknn ---- 

			dcqg_c = (double)0.16999052179242813e0 * dcqg_b;
			//y=b*(y+.32607257743127307d0*(f(a+c)+f(a-c))) 
			//s=a+c
			cknn_s = dcqg_a + dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			//cknn_cbnd1 = ( dcos(ka*cknn_rr) - cj*dsin(ka*cknn_rr) - cone )/cknn_rr
			//cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))         !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//cv = Z_MUL( Z_MAKE(cknn_s/dl,(double)0.0e0),cv );
			cuDoubleComplex cvv1b = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)               ! pot. scalarny 
			//cv = Z_ADD( cknn_cbnd1,Z_MAKE(cknn_bnd2,(double)0.0e0) );
			cuDoubleComplex csv1b = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			// ---- cknn ----
			//s=a-c
			cknn_s = dcqg_a - dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			//cknn_cbnd1 = ( dcos(ka*cknn_rr) - cj*dsin(ka*cknn_rr) - cone )/cknn_rr
			//cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))            !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			cv = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			cvv1b = Z_MUL(Z_MAKE((double)0.32607257743127307e0, (double)0.0e0), Z_ADD(cv, cvv1b));
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)                      ! pot. scalarny
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//csv1b = .32607257743127307d0*( cv + csv1b ) 
			csv1b = Z_MUL(Z_MAKE((double)0.32607257743127307e0, (double)0.0e0), Z_ADD(cv, csv1b));
			// ---- cknn ---- 
			//cv1v = dcqg_b*( cvv1a+cvv1b )
			//cv1s = dcqg_b*( csv1a+csv1b ) 
			cuDoubleComplex cv1v = Z_MUL(Z_MAKE(dcqg_b, (double)0.0e0), Z_ADD(cvv1a, cvv1b));
			cuDoubleComplex cv1s = Z_MUL(Z_MAKE(dcqg_b, (double)0.0e0), Z_ADD(csv1a, csv1b));
			// --------------- cv1 ------------------


			//call dcqg(dlh,dl,cknn,cv2,4)	      
			//nord=4
			// --------------- cv2 ------------------
			dcqg_a = (double)0.5e0*(dlh + dl);
			dcqg_b = dl - dlh;
			//
			dcqg_c = (double)0.43056815579702629e0 * dcqg_b;
			//y=.17392742256872693d0*(f(a+c)+f(a-c))
			//s=a+c
			cknn_s = dcqg_a + dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			///cknn_cbnd1 = ( dcos(ka*cknn_rr) - cj*dsin(ka*cknn_rr) - cone )/cknn_rr
			///cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))        !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			cv = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			cvv1a = cv;
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)          ! pot. scalarny 
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			csv1a = cv;
			// ---- cknn ----
			//s=a-c
			cknn_s = dcqg_a - dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			//cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//!cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))        !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			cv = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			cvv1a = Z_MUL(Z_MAKE((double)0.17392742256872693e0, (double)0.0e0), Z_ADD(cv, cvv1a));
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)            ! pot. scalarny 
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//csv1a = .17392742256872693d0*( cv + csv1a )
			csv1a = Z_MUL(Z_MAKE((double)0.17392742256872693e0, (double)0.0e0), Z_ADD(cv, csv1a));
			// ---- cknn ---- 

			dcqg_c = (double)0.16999052179242813e0 * dcqg_b;
			//y=b*(y+.32607257743127307d0*(f(a+c)+f(a-c))) 
			//s=a+c
			cknn_s = dcqg_a + dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			//cknn_cbnd1 = ( dcos(ka*cknn_rr) - cj*dsin(ka*cknn_rr) - cone )/cknn_rr
			//cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))         !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//cv = Z_MUL( Z_MAKE(cknn_s/dl,(double)0.0e0),cv );
			cvv1b = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)               ! pot. scalarny 
			//cv = Z_ADD( cknn_cbnd1,Z_MAKE(cknn_bnd2,(double)0.0e0) );
			csv1b = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			// ---- cknn ----
			//s=a-c
			cknn_s = dcqg_a - dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			//cknn_cbnd1 = ( dcos(ka*cknn_rr) - cj*dsin(ka*cknn_rr) - cone )/cknn_rr
			//cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))            !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			cv = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			cvv1b = Z_MUL(Z_MAKE((double)0.32607257743127307e0, (double)0.0e0), Z_ADD(cv, cvv1b));
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)                      ! pot. scalarny
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//csv1b = .32607257743127307d0*( cv + csv1b ) 
			csv1b = Z_MUL(Z_MAKE((double)0.32607257743127307e0, (double)0.0e0), Z_ADD(cv, csv1b));
			// ---- cknn ---- 
			//cv2v = dcqg_b*( cvv1a+cvv1b )
			//cv2s = dcqg_b*( csv1a+csv1b  
			cuDoubleComplex cv2v = Z_MUL(Z_MAKE(dcqg_b, (double)0.0e0), Z_ADD(cvv1a, cvv1b));
			cuDoubleComplex cv2s = Z_MUL(Z_MAKE(dcqg_b, (double)0.0e0), Z_ADD(csv1a, csv1b));
			// --------------- cv2 ------------------

			//res_cv = 0.5d0*dl*( 1.d0+DLOG(16.d0*an/dl) )/(pi*an)	! dla ifun=1 i ifun=2
			//res_cs = dl*( 1.d0+DLOG(16.d0*an/dl) )/(pi*an)	        ! dla ifun=3
			double res_cs = dl*((double)1.0e0 + log((double)16.0e0*an / dl)) / (pi*an);
			double res_cv = (double)0.5e0*res_cs;
			//cvval1 = cv1v + cv2v + DCMPLX(res_cv,0.d0)     
			//csval1 = cv1s + cv2s + DCMPLX(res_cs,0.d0)
			cvval1 = Z_ADD(cv2v, Z_MAKE(res_cv, (double)0.0e0));
			cvval1 = Z_ADD(cv1v, cvval1);
			csval1 = Z_ADD(cv2s, Z_MAKE(res_cs, (double)0.0e0));
			csval1 = Z_ADD(cv1s, csval1);
		}
		else {
			//call dcqg(0.d0,dl,ckmn,cv,2)
			//cvval1 = dcmplx(0.0d0,0.0d0)
			//csval1 = dcmplx(0.0d0,0.0d0)
			cvval1 = Z_MAKE((double)0.0e0, (double)0.0e0);
			csval1 = Z_MAKE((double)0.0e0, (double)0.0e0);
		}
		// ------------ Zamiast prcodeury cwire1 -------------------

		//tvecl1 = dc%x*lmc%x+dc%y*lmc%y+dc%z*lmc%z       ! obliczane !lm+ * 1n+ dla k=1, obliczane !lm- * 1n+ dla k=1
		double tvecl1 = dc.x*lmc.x + dc.y*lmc.y + dc.z*lmc.z;
		//csval1 = csval1/dl
		csval1 = Z_DIV(csval1, Z_MAKE(dl, (double)0.0e0));

		//cvval = dcmplx(tvecl1,0.d0)*cvval1
		cuDoubleComplex cvval = Z_MUL(cvval1, Z_MAKE(tvecl1, (double)0.0e0));
		//csval = csval1
		cuDoubleComplex csval = csval1;

		//cz_ij = cz_ij+ka*ka*cvval + znak * csval
		cuDoubleComplex cz_ij = Z_MUL(Z_MAKE((double)znak, (double)0.0e0), csval);
		cz_ij = Z_ADD(Z_MUL(Z_MAKE(ka*ka, (double)0.0e0), cvval), cz_ij);

		//cz_ij = cz_ij*cj*eta0/(4.0e0*pi*ka)
		cz_ij = Z_MUL(cz_ij, Z_MAKE(eta0 / ((double)4.0e0*pi*ka), (double)0.0e0));
		cz_ij = Z_MUL(cz_ij, Z_ONE_IMG);
		//cuDoubleComplex cz_ij1 = cz[ i+j*n ];

		//cz_d[i_m+j_n*n] = Z_MAKE(10.0f,0.0f);
		//cz_d[i_m+j_n*n] = Z_MUL ( cz_d[i_m+j_n*n] , (Z_MAKE(10.0f,0.0f));
		//cz_d[i_m+j_n*n] = Z_MAKE((double)i_m,(double)j_n);
		//cz_d[i_m+j_n*n]=make_cuDoubleComplex( (i_m),(j_n) );
		//cz_d[ i+j*n ] =  Z_MAKE(float(i*n+j),0.0f) ;
		//cz[ i+j*n ] = Z_MAKE(x[mp-1],z[mp-1]);
		//cz[ i+j*n ] = Z_MAKE(lmc.x,lmc.z);
		//cz[ i+j*n ] = cknn_cbnd1;
		//cz[ i+j*n ] = Z_MAKE(cknn_res,cknn_bnd2);
		cz[i + (j + kk)*n] = Z_ADD(cz[i + (j + kk)*n], cz_ij);
		//cz[ i+j*n ] = Z_ADD( cz[ i+j*n ],Z_MAKE((double)10.0f,(double)0.0f) );
		//cz[ i+j*n ] = Z_ADD(cz_ij1, cz_ij );
		//cz[ i+j*n ] =  Z_ADD(c1.z,Z_MAKE(float(i*n+j),0.0f)) ;
	}
}
//----------------------------------------------------------------------
__global__ void wim_gpu_cknn_cwire1_im1(cuDoubleComplex *cz, float *x, float *y, float *z,
	float *rad, int *ipc, int *ips, int *ise, int nsmax,
	int n, double ka, double eta0,
	float sgnx, float sgny, float sgnz, float sgnenv,
	int realSubSliceSize, int k, int kk)
{

	double pi = acos(double(-1.e0));
	cuDoubleComplex cvval1, csval1;

	// Block index
	int bx = blockIdx.x;
	int by = blockIdx.y;

	// Thread index
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	// Accumulate row i of A and column j of B
	int i = tx + bx * blockDim.x;
	int j = ty + by * blockDim.y;

	//c16vec c1 = make_c16vec(Z_ONE,Z_MAKE(0.0f, 1.0f),Z_MUL(Z_MAKE(2.0f, 0.0f),Z_ONE));

	if (i < n  && j < realSubSliceSize) {
		int mp = abs(ipc[i]);
		int ms = ips[i];
		int mlr = ise[ms - 1];
		//double am = rad[ms-1];
		double3 rmp = make_double3((double)x[mp - 1], (double)y[mp - 1], (double)z[mp - 1]);
		double3 rmlr = make_double3((double)x[mlr - 1], (double)y[mlr - 1], (double)z[mlr - 1]);
		double3 rmc = make_double3((double)0.5e0*(rmlr.x + rmp.x), (double)0.5e0*(rmlr.y + rmp.y), (double)0.5e0*(rmlr.z + rmp.z));
		float znak = -1.0f;
		double3 lmc = make_double3((double)znak*(rmc.x - rmp.x), (double)znak*(rmc.y - rmp.y), (double)znak*(rmc.z - rmp.z));

		// ------------ Zamiast prcodeury cwire1 -------------------
		//ifun=1 i ifun=2=> res_cv
		//fun=3 => res_cs		  
		//call cwire1_gpu(ms,n,cvval1)             ! obliczanie calki na Sn+  pot.wek.
		//call cwire1(ms,n,csval1)

		//int np = abs(ipc[k-1+kk-1+j]);
		//int ns = ips[k-1+kk-1+j];
		int np = abs(ipc[k + kk + j]);
		int ns = ips[k + kk + j];
		double3 rnp = make_double3((double)sgnx*(double)x[np - 1], (double)sgny*(double)y[np - 1], (double)sgnz*(double)z[np - 1]);
		int nlr = ise[ns - 1];
		double an = rad[ns - 1];
		double3 rn = make_double3((double)sgnx*(double)x[nlr - 1], (double)sgny*(double)y[nlr - 1], (double)sgnz*(double)z[nlr - 1]);
		double3 rd = make_double3(rnp.x - rn.x, rnp.y - rn.y, rnp.z - rn.z);
		double dl = sqrt(rd.x*rd.x + rd.y*rd.y + rd.z*rd.z);
		double dlh = (double)0.5e0*dl;
		double3 dc = make_double3(rd.x / dl, rd.y / dl, rd.z / dl);

		int ise_ms1 = ise[ms - 1];
		int ise_ms2 = ise[ms - 1 + nsmax];
		int ise_ns1 = ise[ns - 1];
		int ise_ns2 = ise[ns - 1 + nsmax];
		//float z_ise_ms1 = z[ise[ms-1]];
		//float z_ise_ms2 = z[ise[ms-1+nsmax]];
		float z_ise_ms1 = z[ise[ms - 1] - 1];
		float z_ise_ms2 = z[ise[ms - 1 + nsmax] - 1];

		if (war_gpu_pec_im(ise_ms1, ise_ms2, ise_ns1, ise_ns2, z_ise_ms1, z_ise_ms2)) {
			//call dcqg(0.d0,dlh,cknn,cv1,4)
			//nord=4
			//  res_cv = 0.5d0*dl*( 1.d0+DLOG(16.d0*an/dl) )/(pi*an)	! dla ifun=1
			//  res_cs = dl*( 1.d0+DLOG(16.d0*an/dl) )/(pi*an)	        ! dla ifun=3		    
			// --------------- cv1 ------------------
			//dcqg_a = 0.5d0*(0.0d0+dlh)
			//dcqg_b = dlh-0.0d0
			double dcqg_a = (double)0.5e0*dlh;
			double dcqg_b = dlh;
			//
			double dcqg_c = (double)0.43056815579702629e0 * dcqg_b;
			//y=.17392742256872693d0*(f(a+c)+f(a-c))
			//s=a+c
			double cknn_s = dcqg_a + dcqg_c;
			// ---- cknn ----
			double cknn_r = (double)0.5e0*dl - cknn_s;
			double cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			///cknn_cbnd1 = ( dcos(ka*cknn_rr) - cj*dsin(ka*cknn_rr) - cone )/cknn_rr
			///cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cuDoubleComplex cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			double cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			double cknn_res = delick_gpu(cknn_bet);
			double cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))        !1 segment n-tej f.bazowej
			cuDoubleComplex cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//cv = Z_MUL( Z_MAKE(cknn_s/dl,(double)0.0e0),cv );
			cuDoubleComplex cvv1a = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)          ! pot. scalarny 
			//cv = Z_ADD( cknn_cbnd1,Z_MAKE(cknn_bnd2,(double)0.0e0) );
			cuDoubleComplex csv1a = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			// ---- cknn ----
			//s=a-c
			cknn_s = dcqg_a - dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			//cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//!cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))        !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			cv = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			cvv1a = Z_MUL(Z_MAKE((double)0.17392742256872693e0, (double)0.0e0), Z_ADD(cv, cvv1a));
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)            ! pot. scalarny 
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//csv1a = .17392742256872693d0*( cv + csv1a )
			csv1a = Z_MUL(Z_MAKE((double)0.17392742256872693e0, (double)0.0e0), Z_ADD(cv, csv1a));
			// ---- cknn ---- 

			dcqg_c = (double)0.16999052179242813e0 * dcqg_b;
			//y=b*(y+.32607257743127307d0*(f(a+c)+f(a-c))) 
			//s=a+c
			cknn_s = dcqg_a + dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			//cknn_cbnd1 = ( dcos(ka*cknn_rr) - cj*dsin(ka*cknn_rr) - cone )/cknn_rr
			//cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))         !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//cv = Z_MUL( Z_MAKE(cknn_s/dl,(double)0.0e0),cv );
			cuDoubleComplex cvv1b = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)               ! pot. scalarny 
			//cv = Z_ADD( cknn_cbnd1,Z_MAKE(cknn_bnd2,(double)0.0e0) );
			cuDoubleComplex csv1b = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			// ---- cknn ----
			//s=a-c
			cknn_s = dcqg_a - dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			//cknn_cbnd1 = ( dcos(ka*cknn_rr) - cj*dsin(ka*cknn_rr) - cone )/cknn_rr
			//cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))            !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			cv = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			cvv1b = Z_MUL(Z_MAKE((double)0.32607257743127307e0, (double)0.0e0), Z_ADD(cv, cvv1b));
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)                      ! pot. scalarny
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//csv1b = .32607257743127307d0*( cv + csv1b ) 
			csv1b = Z_MUL(Z_MAKE((double)0.32607257743127307e0, (double)0.0e0), Z_ADD(cv, csv1b));
			// ---- cknn ---- 
			//cv1v = dcqg_b*( cvv1a+cvv1b )
			//cv1s = dcqg_b*( csv1a+csv1b ) 
			cuDoubleComplex cv1v = Z_MUL(Z_MAKE(dcqg_b, (double)0.0e0), Z_ADD(cvv1a, cvv1b));
			cuDoubleComplex cv1s = Z_MUL(Z_MAKE(dcqg_b, (double)0.0e0), Z_ADD(csv1a, csv1b));
			// --------------- cv1 ------------------


			//call dcqg(dlh,dl,cknn,cv2,4)	      
			//nord=4
			// --------------- cv2 ------------------
			dcqg_a = (double)0.5e0*(dlh + dl);
			dcqg_b = dl - dlh;
			//
			dcqg_c = (double)0.43056815579702629e0 * dcqg_b;
			//y=.17392742256872693d0*(f(a+c)+f(a-c))
			//s=a+c
			cknn_s = dcqg_a + dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			///cknn_cbnd1 = ( dcos(ka*cknn_rr) - cj*dsin(ka*cknn_rr) - cone )/cknn_rr
			///cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))        !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			cv = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			cvv1a = cv;
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)          ! pot. scalarny 
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			csv1a = cv;
			// ---- cknn ----
			//s=a-c
			cknn_s = dcqg_a - dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			//cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//!cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))        !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			cv = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			cvv1a = Z_MUL(Z_MAKE((double)0.17392742256872693e0, (double)0.0e0), Z_ADD(cv, cvv1a));
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)            ! pot. scalarny 
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//csv1a = .17392742256872693d0*( cv + csv1a )
			csv1a = Z_MUL(Z_MAKE((double)0.17392742256872693e0, (double)0.0e0), Z_ADD(cv, csv1a));
			// ---- cknn ---- 

			dcqg_c = (double)0.16999052179242813e0 * dcqg_b;
			//y=b*(y+.32607257743127307d0*(f(a+c)+f(a-c))) 
			//s=a+c
			cknn_s = dcqg_a + dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			//cknn_cbnd1 = ( dcos(ka*cknn_rr) - cj*dsin(ka*cknn_rr) - cone )/cknn_rr
			//cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))         !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//cv = Z_MUL( Z_MAKE(cknn_s/dl,(double)0.0e0),cv );
			cvv1b = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)               ! pot. scalarny 
			//cv = Z_ADD( cknn_cbnd1,Z_MAKE(cknn_bnd2,(double)0.0e0) );
			csv1b = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			// ---- cknn ----
			//s=a-c
			cknn_s = dcqg_a - dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			//cknn_cbnd1 = ( dcos(ka*cknn_rr) - cj*dsin(ka*cknn_rr) - cone )/cknn_rr
			//cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))            !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			cv = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			cvv1b = Z_MUL(Z_MAKE((double)0.32607257743127307e0, (double)0.0e0), Z_ADD(cv, cvv1b));
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)                      ! pot. scalarny
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//csv1b = .32607257743127307d0*( cv + csv1b ) 
			csv1b = Z_MUL(Z_MAKE((double)0.32607257743127307e0, (double)0.0e0), Z_ADD(cv, csv1b));
			// ---- cknn ---- 
			//cv2v = dcqg_b*( cvv1a+cvv1b )
			//cv2s = dcqg_b*( csv1a+csv1b  
			cuDoubleComplex cv2v = Z_MUL(Z_MAKE(dcqg_b, (double)0.0e0), Z_ADD(cvv1a, cvv1b));
			cuDoubleComplex cv2s = Z_MUL(Z_MAKE(dcqg_b, (double)0.0e0), Z_ADD(csv1a, csv1b));
			// --------------- cv2 ------------------

			//res_cv = 0.5d0*dl*( 1.d0+DLOG(16.d0*an/dl) )/(pi*an)	! dla ifun=1 i ifun=2
			//res_cs = dl*( 1.d0+DLOG(16.d0*an/dl) )/(pi*an)	        ! dla ifun=3
			double res_cs = dl*((double)1.0e0 + log((double)16.0e0*an / dl)) / (pi*an);
			double res_cv = (double)0.5e0*res_cs;
			//cvval1 = cv1v + cv2v + DCMPLX(res_cv,0.d0)     
			//csval1 = cv1s + cv2s + DCMPLX(res_cs,0.d0)
			cvval1 = Z_ADD(cv2v, Z_MAKE(res_cv, (double)0.0e0));
			cvval1 = Z_ADD(cv1v, cvval1);
			csval1 = Z_ADD(cv2s, Z_MAKE(res_cs, (double)0.0e0));
			csval1 = Z_ADD(cv1s, csval1);
		}
		else {
			//call dcqg(0.d0,dl,ckmn,cv,2)
			//cvval1 = dcmplx(0.0d0,0.0d0)
			//csval1 = dcmplx(0.0d0,0.0d0)
			cvval1 = Z_MAKE((double)0.0e0, (double)0.0e0);
			csval1 = Z_MAKE((double)0.0e0, (double)0.0e0);
		}
		// ------------ Zamiast prcodeury cwire1 -------------------

		//tvecl1 = dc%x*lmc%x+dc%y*lmc%y+dc%z*lmc%z       ! obliczane !lm+ * 1n+ dla k=1, obliczane !lm- * 1n+ dla k=1
		double tvecl1 = dc.x*lmc.x + dc.y*lmc.y + dc.z*lmc.z;
		//csval1 = csval1/dl
		csval1 = Z_DIV(csval1, Z_MAKE(dl, (double)0.0e0));

		//cvval = sgnenv*dcmplx(tvecl1,0.d0)*cvval1
		cuDoubleComplex cvval = Z_MUL(cvval1, Z_MAKE((double)sgnenv*tvecl1, (double)0.0e0));
		//csval = csval1*sgnenv
		cuDoubleComplex csval = Z_MUL(Z_MAKE((double)sgnenv, (double)0.0e0), csval1);

		//cz_ij = cz_ij+ka*ka*cvval + znak * csval
		cuDoubleComplex cz_ij = Z_MUL(Z_MAKE((double)znak, (double)0.0e0), csval);
		cz_ij = Z_ADD(Z_MUL(Z_MAKE(ka*ka, (double)0.0e0), cvval), cz_ij);

		//cz_ij = cz_ij*cj*eta0/(4.0e0*pi*ka)
		cz_ij = Z_MUL(cz_ij, Z_MAKE(eta0 / ((double)4.0e0*pi*ka), (double)0.0e0));
		cz_ij = Z_MUL(cz_ij, Z_ONE_IMG);


		//cz_d[i_m+j_n*n] = Z_MAKE(10.0f,0.0f);
		//cz_d[i_m+j_n*n] = Z_MUL ( cz_d[i_m+j_n*n] , (Z_MAKE(10.0f,0.0f));
		//cz_d[i_m+j_n*n] = Z_MAKE((double)i_m,(double)j_n);
		//cz_d[i_m+j_n*n]=make_cuDoubleComplex( (i_m),(j_n) );
		//cz_d[ i+j*n ] =  Z_MAKE(float(i*n+j),0.0f) ;
		//cz[ i+j*n ] = Z_MAKE(x[mp-1],z[mp-1]);
		//cz[ i+j*n ] = Z_MAKE(dc.x,dc.z);
		//cz[ i+j*n ] = cknn_cbnd1;
		//cz[ i+j*n ] = Z_MAKE(cknn_res,cknn_bnd2);
		cz[i + (j + kk)*n] = Z_ADD(cz[i + (j + kk)*n], cz_ij);
		//cz[ i+j*n ] =  Z_ADD(c1.z,Z_MAKE(float(i*n+j),0.0f)) ;
	}
}

//----------------------------------------------------------------------

__global__ void wim_gpu_cknn_cwire1_im2(cuDoubleComplex *cz, float *x, float *y, float *z,
	float *rad, int *ipc, int *ips, int *ise, int nsmax,
	int n, double ka, double eta0,
	float sgnx, float sgny, float sgnz, float sgnenv,
	int realSubSliceSize, int k, int kk)
{

	double pi = acos(double(-1.e0));
	cuDoubleComplex cvval1, csval1;

	// Block index
	int bx = blockIdx.x;
	int by = blockIdx.y;

	// Thread index
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	// Accumulate row i of A and column j of B
	int i = tx + bx * blockDim.x;
	int j = ty + by * blockDim.y;

	//c16vec c1 = make_c16vec(Z_ONE,Z_MAKE(0.0f, 1.0f),Z_MUL(Z_MAKE(2.0f, 0.0f),Z_ONE));

	if (i < n  && j < realSubSliceSize) {
		int mp = abs(ipc[i]);
		//int ms = ips[i];
		//int mlr = ise[ms-1];
		int ms = ips[i + nsmax];
		int mlr = ise[ms - 1 + nsmax];
		//double am = rad[ms-1];
		double3 rmp = make_double3((double)x[mp - 1], (double)y[mp - 1], (double)z[mp - 1]);
		double3 rmlr = make_double3((double)x[mlr - 1], (double)y[mlr - 1], (double)z[mlr - 1]);
		double3 rmc = make_double3((double)0.5e0*(rmlr.x + rmp.x), (double)0.5e0*(rmlr.y + rmp.y), (double)0.5e0*(rmlr.z + rmp.z));
		float znak = 1.0f;
		double3 lmc = make_double3((double)znak*(rmc.x - rmp.x), (double)znak*(rmc.y - rmp.y), (double)znak*(rmc.z - rmp.z));

		// ------------ Zamiast prcodeury cwire1 -------------------
		//ifun=1 i ifun=2=> res_cv
		//fun=3 => res_cs		  
		//call cwire1_gpu(ms,n,cvval1)             ! obliczanie calki na Sn+  pot.wek.
		//call cwire1(ms,n,csval1)

		//int np = abs(ipc[k-1+kk-1+j]);
		//int ns = ips[k-1+kk-1+j];
		int np = abs(ipc[k + kk + j]);
		int ns = ips[k + kk + j];
		double3 rnp = make_double3((double)sgnx*(double)x[np - 1], (double)sgny*(double)y[np - 1], (double)sgnz*(double)z[np - 1]);
		int nlr = ise[ns - 1];
		double an = rad[ns - 1];
		double3 rn = make_double3((double)sgnx*(double)x[nlr - 1], (double)sgny*(double)y[nlr - 1], (double)sgnz*(double)z[nlr - 1]);
		double3 rd = make_double3(rnp.x - rn.x, rnp.y - rn.y, rnp.z - rn.z);
		double dl = sqrt(rd.x*rd.x + rd.y*rd.y + rd.z*rd.z);
		double dlh = (double)0.5e0*dl;
		double3 dc = make_double3(rd.x / dl, rd.y / dl, rd.z / dl);

		int ise_ms1 = ise[ms - 1];
		int ise_ms2 = ise[ms - 1 + nsmax];
		int ise_ns1 = ise[ns - 1];
		int ise_ns2 = ise[ns - 1 + nsmax];
		//float z_ise_ms1 = z[ise[ms-1]];
		//float z_ise_ms2 = z[ise[ms-1+nsmax]];
		float z_ise_ms1 = z[ise[ms - 1] - 1];
		float z_ise_ms2 = z[ise[ms - 1 + nsmax] - 1];

		if (war_gpu_pec_im(ise_ms1, ise_ms2, ise_ns1, ise_ns2, z_ise_ms1, z_ise_ms2)) {
			//call dcqg(0.d0,dlh,cknn,cv1,4)
			//nord=4
			//  res_cv = 0.5d0*dl*( 1.d0+DLOG(16.d0*an/dl) )/(pi*an)	! dla ifun=1
			//  res_cs = dl*( 1.d0+DLOG(16.d0*an/dl) )/(pi*an)	        ! dla ifun=3		    
			// --------------- cv1 ------------------
			//dcqg_a = 0.5d0*(0.0d0+dlh)
			//dcqg_b = dlh-0.0d0
			double dcqg_a = (double)0.5e0*dlh;
			double dcqg_b = dlh;
			//
			double dcqg_c = (double)0.43056815579702629e0 * dcqg_b;
			//y=.17392742256872693d0*(f(a+c)+f(a-c))
			//s=a+c
			double cknn_s = dcqg_a + dcqg_c;
			// ---- cknn ----
			double cknn_r = (double)0.5e0*dl - cknn_s;
			double cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			///cknn_cbnd1 = ( dcos(ka*cknn_rr) - cj*dsin(ka*cknn_rr) - cone )/cknn_rr
			///cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cuDoubleComplex cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			double cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			double cknn_res = delick_gpu(cknn_bet);
			double cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))        !1 segment n-tej f.bazowej
			cuDoubleComplex cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//cv = Z_MUL( Z_MAKE(cknn_s/dl,(double)0.0e0),cv );
			cuDoubleComplex cvv1a = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)          ! pot. scalarny 
			//cv = Z_ADD( cknn_cbnd1,Z_MAKE(cknn_bnd2,(double)0.0e0) );
			cuDoubleComplex csv1a = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			// ---- cknn ----
			//s=a-c
			cknn_s = dcqg_a - dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			//cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//!cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))        !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			cv = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			cvv1a = Z_MUL(Z_MAKE((double)0.17392742256872693e0, (double)0.0e0), Z_ADD(cv, cvv1a));
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)            ! pot. scalarny 
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//csv1a = .17392742256872693d0*( cv + csv1a )
			csv1a = Z_MUL(Z_MAKE((double)0.17392742256872693e0, (double)0.0e0), Z_ADD(cv, csv1a));
			// ---- cknn ---- 

			dcqg_c = (double)0.16999052179242813e0 * dcqg_b;
			//y=b*(y+.32607257743127307d0*(f(a+c)+f(a-c))) 
			//s=a+c
			cknn_s = dcqg_a + dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			//cknn_cbnd1 = ( dcos(ka*cknn_rr) - cj*dsin(ka*cknn_rr) - cone )/cknn_rr
			//cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))         !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//cv = Z_MUL( Z_MAKE(cknn_s/dl,(double)0.0e0),cv );
			cuDoubleComplex cvv1b = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)               ! pot. scalarny 
			//cv = Z_ADD( cknn_cbnd1,Z_MAKE(cknn_bnd2,(double)0.0e0) );
			cuDoubleComplex csv1b = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			// ---- cknn ----
			//s=a-c
			cknn_s = dcqg_a - dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			//cknn_cbnd1 = ( dcos(ka*cknn_rr) - cj*dsin(ka*cknn_rr) - cone )/cknn_rr
			//cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))            !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			cv = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			cvv1b = Z_MUL(Z_MAKE((double)0.32607257743127307e0, (double)0.0e0), Z_ADD(cv, cvv1b));
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)                      ! pot. scalarny
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//csv1b = .32607257743127307d0*( cv + csv1b ) 
			csv1b = Z_MUL(Z_MAKE((double)0.32607257743127307e0, (double)0.0e0), Z_ADD(cv, csv1b));
			// ---- cknn ---- 
			//cv1v = dcqg_b*( cvv1a+cvv1b )
			//cv1s = dcqg_b*( csv1a+csv1b ) 
			cuDoubleComplex cv1v = Z_MUL(Z_MAKE(dcqg_b, (double)0.0e0), Z_ADD(cvv1a, cvv1b));
			cuDoubleComplex cv1s = Z_MUL(Z_MAKE(dcqg_b, (double)0.0e0), Z_ADD(csv1a, csv1b));
			// --------------- cv1 ------------------


			//call dcqg(dlh,dl,cknn,cv2,4)	      
			//nord=4
			// --------------- cv2 ------------------
			dcqg_a = (double)0.5e0*(dlh + dl);
			dcqg_b = dl - dlh;
			//
			dcqg_c = (double)0.43056815579702629e0 * dcqg_b;
			//y=.17392742256872693d0*(f(a+c)+f(a-c))
			//s=a+c
			cknn_s = dcqg_a + dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			///cknn_cbnd1 = ( dcos(ka*cknn_rr) - cj*dsin(ka*cknn_rr) - cone )/cknn_rr
			///cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))        !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			cv = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			cvv1a = cv;
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)          ! pot. scalarny 
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			csv1a = cv;
			// ---- cknn ----
			//s=a-c
			cknn_s = dcqg_a - dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			//cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//!cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))        !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			cv = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			cvv1a = Z_MUL(Z_MAKE((double)0.17392742256872693e0, (double)0.0e0), Z_ADD(cv, cvv1a));
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)            ! pot. scalarny 
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//csv1a = .17392742256872693d0*( cv + csv1a )
			csv1a = Z_MUL(Z_MAKE((double)0.17392742256872693e0, (double)0.0e0), Z_ADD(cv, csv1a));
			// ---- cknn ---- 

			dcqg_c = (double)0.16999052179242813e0 * dcqg_b;
			//y=b*(y+.32607257743127307d0*(f(a+c)+f(a-c))) 
			//s=a+c
			cknn_s = dcqg_a + dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			//cknn_cbnd1 = ( dcos(ka*cknn_rr) - cj*dsin(ka*cknn_rr) - cone )/cknn_rr
			//cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))         !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//cv = Z_MUL( Z_MAKE(cknn_s/dl,(double)0.0e0),cv );
			cvv1b = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)               ! pot. scalarny 
			//cv = Z_ADD( cknn_cbnd1,Z_MAKE(cknn_bnd2,(double)0.0e0) );
			csv1b = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			// ---- cknn ----
			//s=a-c
			cknn_s = dcqg_a - dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			//cknn_cbnd1 = ( dcos(ka*cknn_rr) - cj*dsin(ka*cknn_rr) - cone )/cknn_rr
			//cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))            !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			cv = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			cvv1b = Z_MUL(Z_MAKE((double)0.32607257743127307e0, (double)0.0e0), Z_ADD(cv, cvv1b));
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)                      ! pot. scalarny
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//csv1b = .32607257743127307d0*( cv + csv1b ) 
			csv1b = Z_MUL(Z_MAKE((double)0.32607257743127307e0, (double)0.0e0), Z_ADD(cv, csv1b));
			// ---- cknn ---- 
			//cv2v = dcqg_b*( cvv1a+cvv1b )
			//cv2s = dcqg_b*( csv1a+csv1b  
			cuDoubleComplex cv2v = Z_MUL(Z_MAKE(dcqg_b, (double)0.0e0), Z_ADD(cvv1a, cvv1b));
			cuDoubleComplex cv2s = Z_MUL(Z_MAKE(dcqg_b, (double)0.0e0), Z_ADD(csv1a, csv1b));
			// --------------- cv2 ------------------

			//res_cv = 0.5d0*dl*( 1.d0+DLOG(16.d0*an/dl) )/(pi*an)	! dla ifun=1 i ifun=2
			//res_cs = dl*( 1.d0+DLOG(16.d0*an/dl) )/(pi*an)	        ! dla ifun=3
			double res_cs = dl*((double)1.0e0 + log((double)16.0e0*an / dl)) / (pi*an);
			double res_cv = (double)0.5e0*res_cs;
			//cvval1 = cv1v + cv2v + DCMPLX(res_cv,0.d0)     
			//csval1 = cv1s + cv2s + DCMPLX(res_cs,0.d0)
			cvval1 = Z_ADD(cv2v, Z_MAKE(res_cv, (double)0.0e0));
			cvval1 = Z_ADD(cv1v, cvval1);
			csval1 = Z_ADD(cv2s, Z_MAKE(res_cs, (double)0.0e0));
			csval1 = Z_ADD(cv1s, csval1);
		}
		else {
			//call dcqg(0.d0,dl,ckmn,cv,2)
			//cvval1 = dcmplx(0.0d0,0.0d0)
			//csval1 = dcmplx(0.0d0,0.0d0)
			cvval1 = Z_MAKE((double)0.0e0, (double)0.0e0);
			csval1 = Z_MAKE((double)0.0e0, (double)0.0e0);
		}
		// ------------ Zamiast prcodeury cwire1 -------------------

		//tvecl1 = dc%x*lmc%x+dc%y*lmc%y+dc%z*lmc%z       ! obliczane !lm+ * 1n+ dla k=1, obliczane !lm- * 1n+ dla k=1
		double tvecl1 = dc.x*lmc.x + dc.y*lmc.y + dc.z*lmc.z;
		//csval1 = csval1/dl
		csval1 = Z_DIV(csval1, Z_MAKE(dl, (double)0.0e0));

		//cvval = sgnenv*dcmplx(tvecl1,0.d0)*cvval1
		cuDoubleComplex cvval = Z_MUL(cvval1, Z_MAKE((double)sgnenv*tvecl1, (double)0.0e0));
		//csval = csval1*sgnenv
		cuDoubleComplex csval = Z_MUL(Z_MAKE((double)sgnenv, (double)0.0e0), csval1);

		//cz_ij = cz_ij+ka*ka*cvval + znak * csval
		cuDoubleComplex cz_ij = Z_MUL(Z_MAKE((double)znak, (double)0.0e0), csval);
		cz_ij = Z_ADD(Z_MUL(Z_MAKE(ka*ka, (double)0.0e0), cvval), cz_ij);

		//cz_ij = cz_ij*cj*eta0/(4.0e0*pi*ka)
		cz_ij = Z_MUL(cz_ij, Z_MAKE(eta0 / ((double)4.0e0*pi*ka), (double)0.0e0));
		cz_ij = Z_MUL(cz_ij, Z_ONE_IMG);


		//cz_d[i_m+j_n*n] = Z_MAKE(10.0f,0.0f);
		//cz_d[i_m+j_n*n] = Z_MUL ( cz_d[i_m+j_n*n] , (Z_MAKE(10.0f,0.0f));
		//cz_d[i_m+j_n*n] = Z_MAKE((double)i_m,(double)j_n);
		//cz_d[i_m+j_n*n]=make_cuDoubleComplex( (i_m),(j_n) );
		//cz_d[ i+j*n ] =  Z_MAKE(float(i*n+j),0.0f) ;
		//cz[ i+j*n ] = Z_MAKE(x[mp-1],z[mp-1]);
		//cz[ i+j*n ] = Z_MAKE(lmc.x,lmc.z);
		//cz[ i+j*n ] = cknn_cbnd1;
		//cz[ i+j*n ] = Z_MAKE(cknn_res,cknn_bnd2);
		cz[i + (j + kk)*n] = Z_ADD(cz[i + (j + kk)*n], cz_ij);
		//cz[ i+j*n ] = Z_ADD( cz[ i+j*n ],Z_MAKE((double)10.0f,(double)0.0f) );
		//cz[ i+j*n ] =  Z_ADD(c1.z,Z_MAKE(float(i*n+j),0.0f)) ;
	}
}

//----------------------------------------------------------------------

__global__ void wim_gpu_cknn_cwire2_re1(cuDoubleComplex *cz, float *x, float *y, float *z,
	float *rad, int *ipc, int *ips, int *ise, int nsmax,
	int n, double ka, double eta0,
	int realSubSliceSize, int k, int kk)
{

	double pi = acos(double(-1.e0));
	cuDoubleComplex cvval2, csval2;

	// Block index
	int bx = blockIdx.x;
	int by = blockIdx.y;

	// Thread index
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	// Accumulate row i of A and column j of B
	int i = tx + bx * blockDim.x;
	int j = ty + by * blockDim.y;

	//c16vec c1 = make_c16vec(Z_ONE,Z_MAKE(0.0f, 1.0f),Z_MUL(Z_MAKE(2.0f, 0.0f),Z_ONE));

	if (i < n  && j < realSubSliceSize) {
		int mp = abs(ipc[i]);
		int ms = ips[i];
		int mlr = ise[ms - 1];
		//double am = rad[ms-1];
		double3 rmp = make_double3((double)x[mp - 1], (double)y[mp - 1], (double)z[mp - 1]);
		double3 rmlr = make_double3((double)x[mlr - 1], (double)y[mlr - 1], (double)z[mlr - 1]);
		double3 rmc = make_double3((double)0.5e0*(rmlr.x + rmp.x), (double)0.5e0*(rmlr.y + rmp.y), (double)0.5e0*(rmlr.z + rmp.z));
		float znak = -1.0f;
		double3 lmc = make_double3((double)znak*(rmc.x - rmp.x), (double)znak*(rmc.y - rmp.y), (double)znak*(rmc.z - rmp.z));

		// ------------ Zamiast prcodeury cwire1 -------------------
		//ifun=1 i ifun=2=> res_cv
		//fun=3 => res_cs		  
		//call cwire1_gpu(ms,n,cvval1)             ! obliczanie calki na Sn+  pot.wek.
		//call cwire1(ms,n,csval1)

		//int np = abs(ipc[k-1+kk-1+j]);
		//int ns = ips[k-1+kk-1+j];
		int np = abs(ipc[k + kk + j]);
		int ns = ips[k + kk + j + nsmax];
		int nlr = ise[ns - 1 + nsmax];
		double an = rad[ns - 1];
		double3 rnr = make_double3((double)x[nlr - 1], (double)y[nlr - 1], (double)z[nlr - 1]);
		double3 rn = make_double3((double)x[np - 1], (double)y[np - 1], (double)z[np - 1]);
		double3 rd = make_double3(rnr.x - rn.x, rnr.y - rn.y, rnr.z - rn.z);
		double dl = sqrt(rd.x*rd.x + rd.y*rd.y + rd.z*rd.z);
		double dlh = (double)0.5e0*dl;
		double3 dc = make_double3(rd.x / dl, rd.y / dl, rd.z / dl);

		int ise_ms1 = ise[ms - 1];
		int ise_ms2 = ise[ms - 1 + nsmax];
		int ise_ns1 = ise[ns - 1];
		int ise_ns2 = ise[ns - 1 + nsmax];


		if (war_gpu_pec_re(ise_ms1, ise_ms2, ise_ns1, ise_ns2)) {
			//call dcqg(0.d0,dlh,cknn,cv1,4)
			//nord=4
			//  res_cv = 0.5d0*dl*( 1.d0+DLOG(16.d0*an/dl) )/(pi*an)	! dla ifun=1
			//  res_cs = dl*( 1.d0+DLOG(16.d0*an/dl) )/(pi*an)	        ! dla ifun=3		    
			// --------------- cv1 ------------------
			//dcqg_a = 0.5d0*(0.0d0+dlh)
			//dcqg_b = dlh-0.0d0
			double dcqg_a = (double)0.5e0*dlh;
			double dcqg_b = dlh;
			//
			double dcqg_c = (double)0.43056815579702629e0 * dcqg_b;
			//y=.17392742256872693d0*(f(a+c)+f(a-c))
			//s=a+c
			double cknn_s = dcqg_a + dcqg_c;
			// ---- cknn ----
			double cknn_r = (double)0.5e0*dl - cknn_s;
			double cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			///cknn_cbnd1 = ( dcos(ka*cknn_rr) - cj*dsin(ka*cknn_rr) - cone )/cknn_rr
			///cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cuDoubleComplex cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			double cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			double cknn_res = delick_gpu(cknn_bet);
			double cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))        !1 segment n-tej f.bazowej
			cuDoubleComplex cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//cv = Z_MUL( Z_MAKE(cknn_s/dl,(double)0.0e0),cv );
			cuDoubleComplex cvv1a = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)          ! pot. scalarny 
			//cv = Z_ADD( cknn_cbnd1,Z_MAKE(cknn_bnd2,(double)0.0e0) );
			cuDoubleComplex csv1a = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			// ---- cknn ----
			//s=a-c
			cknn_s = dcqg_a - dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			//cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//!cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))        !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			cv = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			cvv1a = Z_MUL(Z_MAKE((double)0.17392742256872693e0, (double)0.0e0), Z_ADD(cv, cvv1a));
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)            ! pot. scalarny 
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//csv1a = .17392742256872693d0*( cv + csv1a )
			csv1a = Z_MUL(Z_MAKE((double)0.17392742256872693e0, (double)0.0e0), Z_ADD(cv, csv1a));
			// ---- cknn ---- 

			dcqg_c = (double)0.16999052179242813e0 * dcqg_b;
			//y=b*(y+.32607257743127307d0*(f(a+c)+f(a-c))) 
			//s=a+c
			cknn_s = dcqg_a + dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			//cknn_cbnd1 = ( dcos(ka*cknn_rr) - cj*dsin(ka*cknn_rr) - cone )/cknn_rr
			//cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))         !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//cv = Z_MUL( Z_MAKE(cknn_s/dl,(double)0.0e0),cv );
			cuDoubleComplex cvv1b = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)               ! pot. scalarny 
			//cv = Z_ADD( cknn_cbnd1,Z_MAKE(cknn_bnd2,(double)0.0e0) );
			cuDoubleComplex csv1b = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			// ---- cknn ----
			//s=a-c
			cknn_s = dcqg_a - dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			//cknn_cbnd1 = ( dcos(ka*cknn_rr) - cj*dsin(ka*cknn_rr) - cone )/cknn_rr
			//cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))            !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			cv = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			cvv1b = Z_MUL(Z_MAKE((double)0.32607257743127307e0, (double)0.0e0), Z_ADD(cv, cvv1b));
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)                      ! pot. scalarny
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//csv1b = .32607257743127307d0*( cv + csv1b ) 
			csv1b = Z_MUL(Z_MAKE((double)0.32607257743127307e0, (double)0.0e0), Z_ADD(cv, csv1b));
			// ---- cknn ---- 
			//cv1v = dcqg_b*( cvv1a+cvv1b )
			//cv1s = dcqg_b*( csv1a+csv1b ) 
			cuDoubleComplex cv1v = Z_MUL(Z_MAKE(dcqg_b, (double)0.0e0), Z_ADD(cvv1a, cvv1b));
			cuDoubleComplex cv1s = Z_MUL(Z_MAKE(dcqg_b, (double)0.0e0), Z_ADD(csv1a, csv1b));
			//cvval2 = cvv1a;
			// --------------- cv1 ------------------


			//call dcqg(dlh,dl,cknn,cv2,4)	      
			//nord=4
			// --------------- cv2 ------------------
			dcqg_a = (double)0.5e0*(dlh + dl);
			dcqg_b = dl - dlh;
			//
			dcqg_c = (double)0.43056815579702629e0 * dcqg_b;
			//y=.17392742256872693d0*(f(a+c)+f(a-c))
			//s=a+c
			cknn_s = dcqg_a + dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			///cknn_cbnd1 = ( dcos(ka*cknn_rr) - cj*dsin(ka*cknn_rr) - cone )/cknn_rr
			///cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))        !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			cv = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			cvv1a = cv;
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)          ! pot. scalarny 
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			csv1a = cv;
			// ---- cknn ----
			//s=a-c
			cknn_s = dcqg_a - dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			//cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//!cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))        !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			cv = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			cvv1a = Z_MUL(Z_MAKE((double)0.17392742256872693e0, (double)0.0e0), Z_ADD(cv, cvv1a));
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)            ! pot. scalarny 
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//csv1a = .17392742256872693d0*( cv + csv1a )
			csv1a = Z_MUL(Z_MAKE((double)0.17392742256872693e0, (double)0.0e0), Z_ADD(cv, csv1a));
			// ---- cknn ---- 

			dcqg_c = (double)0.16999052179242813e0 * dcqg_b;
			//y=b*(y+.32607257743127307d0*(f(a+c)+f(a-c))) 
			//s=a+c
			cknn_s = dcqg_a + dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			//cknn_cbnd1 = ( dcos(ka*cknn_rr) - cj*dsin(ka*cknn_rr) - cone )/cknn_rr
			//cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))         !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//cv = Z_MUL( Z_MAKE(cknn_s/dl,(double)0.0e0),cv );
			cvv1b = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)               ! pot. scalarny 
			//cv = Z_ADD( cknn_cbnd1,Z_MAKE(cknn_bnd2,(double)0.0e0) );
			csv1b = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			// ---- cknn ----
			//s=a-c
			cknn_s = dcqg_a - dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			//cknn_cbnd1 = ( dcos(ka*cknn_rr) - cj*dsin(ka*cknn_rr) - cone )/cknn_rr
			//cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))            !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			cv = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			cvv1b = Z_MUL(Z_MAKE((double)0.32607257743127307e0, (double)0.0e0), Z_ADD(cv, cvv1b));
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)                      ! pot. scalarny
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//csv1b = .32607257743127307d0*( cv + csv1b ) 
			csv1b = Z_MUL(Z_MAKE((double)0.32607257743127307e0, (double)0.0e0), Z_ADD(cv, csv1b));
			// ---- cknn ---- 
			//cv2v = dcqg_b*( cvv1a+cvv1b )
			//cv2s = dcqg_b*( csv1a+csv1b  
			cuDoubleComplex cv2v = Z_MUL(Z_MAKE(dcqg_b, (double)0.0e0), Z_ADD(cvv1a, cvv1b));
			cuDoubleComplex cv2s = Z_MUL(Z_MAKE(dcqg_b, (double)0.0e0), Z_ADD(csv1a, csv1b));
			// --------------- cv2 ------------------

			//res_cv = 0.5d0*dl*( 1.d0+DLOG(16.d0*an/dl) )/(pi*an)	! dla ifun=1 i ifun=2
			//res_cs = dl*( 1.d0+DLOG(16.d0*an/dl) )/(pi*an)	        ! dla ifun=3
			double res_cs = dl*((double)1.0e0 + log((double)16.0e0*an / dl)) / (pi*an);
			double res_cv = (double)0.5e0*res_cs;
			//cvval2 = cv1v + cv2v + DCMPLX(res_cv,0.d0)     
			//csval2 = cv1s + cv2s + DCMPLX(res_cs,0.d0)
			cvval2 = Z_ADD(cv2v, Z_MAKE(res_cv, (double)0.0e0));
			cvval2 = Z_ADD(cv1v, cvval2);
			csval2 = Z_ADD(cv2s, Z_MAKE(res_cs, (double)0.0e0));
			csval2 = Z_ADD(cv1s, csval2);
		}
		else {
			//call dcqg(0.d0,dl,ckmn,cv,2)
			//cvval2 = dcmplx(0.0d0,0.0d0)
			//csval2 = dcmplx(0.0d0,0.0d0)
			cvval2 = Z_MAKE((double)0.0e0, (double)0.0e0);
			csval2 = Z_MAKE((double)0.0e0, (double)0.0e0);
		}
		// ------------ Zamiast prcodeury cwire2 -------------------

		//tvecl2 = dc%x*lmc%x+dc%y*lmc%y+dc%z*lmc%z      ! obliczane !lm+ * 1n- dla k=1,	obliczane  !lm- * 1n- dla k=2 
		//!cvval = dcmplx(tvecl1,0.d0)*cvval1+dcmplx(tvecl2,0.d0)*cvval2   !A(rm+)	dla k=1, A(rm-) dla k=2
		double tvecl2 = dc.x*lmc.x + dc.y*lmc.y + dc.z*lmc.z;
		//csval2 = csval2/dl
		csval2 = Z_DIV(csval2, Z_MAKE(dl, (double)0.0e0));

		//cvval = dcmplx(tvecl2,0.d0)*cvval2
		cuDoubleComplex cvval = Z_MUL(cvval2, Z_MAKE(tvecl2, (double)0.0e0));
		//csval = -csval2
		cuDoubleComplex csval = Z_MUL(Z_MAKE((double)-1.0e0, (double)0.0e0), csval2);

		//cz_ij = cz_ij+ka*ka*cvval + znak * csval
		cuDoubleComplex cz_ij = Z_MUL(Z_MAKE((double)znak, (double)0.0e0), csval);
		cz_ij = Z_ADD(Z_MUL(Z_MAKE(ka*ka, (double)0.0e0), cvval), cz_ij);

		//cz_ij = cz_ij*cj*eta0/(4.0e0*pi*ka)
		cz_ij = Z_MUL(cz_ij, Z_MAKE(eta0 / ((double)4.0e0*pi*ka), (double)0.0e0));
		cz_ij = Z_MUL(cz_ij, Z_ONE_IMG);


		//cz_d[i_m+j_n*n] = Z_MAKE(10.0f,0.0f);
		//cz_d[i_m+j_n*n] = Z_MUL ( cz_d[i_m+j_n*n] , (Z_MAKE(10.0f,0.0f));
		//cz_d[i_m+j_n*n] = Z_MAKE((double)i_m,(double)j_n);
		//cz_d[i_m+j_n*n]=make_cuDoubleComplex( (i_m),(j_n) );
		//cz_d[ i+j*n ] =  Z_MAKE(float(i*n+j),0.0f) ;
		//cz[ i+j*n ] = Z_MAKE(x[mp-1],z[mp-1]);
		//cz[ i+j*n ] = Z_MAKE(dc.x,dc.z);
		//cz[ i+j*n ] = cknn_cbnd1;
		//cz[ i+j*n ] = Z_MAKE(cknn_res,cknn_bnd2);
		cz[i + (j + kk)*n] = Z_ADD(cz[i + (j + kk)*n], cz_ij);
		//cz[i + (j + kk)*n] = cvval2;
		//cz[ i+j*n ] =  Z_ADD(c1.z,Z_MAKE(float(i*n+j),0.0f)) ;
	}
}

//----------------------------------------------------------------------

__global__ void wim_gpu_cknn_cwire2_re2(cuDoubleComplex *cz, float *x, float *y, float *z,
	float *rad, int *ipc, int *ips, int *ise, int nsmax,
	int n, double ka, double eta0,
	int realSubSliceSize, int k, int kk)
{

	double pi = acos(double(-1.e0));
	cuDoubleComplex cvval2, csval2;

	// Block index
	int bx = blockIdx.x;
	int by = blockIdx.y;

	// Thread index
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	// Accumulate row i of A and column j of B
	int i = tx + bx * blockDim.x;
	int j = ty + by * blockDim.y;

	//c16vec c1 = make_c16vec(Z_ONE,Z_MAKE(0.0f, 1.0f),Z_MUL(Z_MAKE(2.0f, 0.0f),Z_ONE));

	if (i < n  && j < realSubSliceSize) {
		int mp = abs(ipc[i]);
		//int ms = ips[i];
		//int mlr = ise[ms-1];
		int ms = ips[i + nsmax];
		int mlr = ise[ms - 1 + nsmax];
		//double am = rad[ms-1];
		double3 rmp = make_double3((double)x[mp - 1], (double)y[mp - 1], (double)z[mp - 1]);
		double3 rmlr = make_double3((double)x[mlr - 1], (double)y[mlr - 1], (double)z[mlr - 1]);
		double3 rmc = make_double3((double)0.5e0*(rmlr.x + rmp.x), (double)0.5e0*(rmlr.y + rmp.y), (double)0.5e0*(rmlr.z + rmp.z));
		float znak = 1.0f;
		double3 lmc = make_double3((double)znak*(rmc.x - rmp.x), (double)znak*(rmc.y - rmp.y), (double)znak*(rmc.z - rmp.z));

		// ------------ Zamiast prcodeury cwire1 -------------------
		//ifun=1 i ifun=2=> res_cv
		//fun=3 => res_cs		  
		//call cwire1_gpu(ms,n,cvval1)             ! obliczanie calki na Sn+  pot.wek.
		//call cwire1(ms,n,csval1)

		//int np = abs(ipc[k-1+kk-1+j]);
		//int ns = ips[k-1+kk-1+j];
		int np = abs(ipc[k + kk + j]);
		int ns = ips[k + kk + j + nsmax];
		int nlr = ise[ns - 1 + nsmax];
		double an = rad[ns - 1];
		double3 rnr = make_double3((double)x[nlr - 1], (double)y[nlr - 1], (double)z[nlr - 1]);
		double3 rn = make_double3((double)x[np - 1], (double)y[np - 1], (double)z[np - 1]);
		double3 rd = make_double3(rnr.x - rn.x, rnr.y - rn.y, rnr.z - rn.z);
		double dl = sqrt(rd.x*rd.x + rd.y*rd.y + rd.z*rd.z);
		double dlh = (double)0.5e0*dl;
		double3 dc = make_double3(rd.x / dl, rd.y / dl, rd.z / dl);

		int ise_ms1 = ise[ms - 1];
		int ise_ms2 = ise[ms - 1 + nsmax];
		int ise_ns1 = ise[ns - 1];
		int ise_ns2 = ise[ns - 1 + nsmax];


		if (war_gpu_pec_re(ise_ms1, ise_ms2, ise_ns1, ise_ns2)) {
			//call dcqg(0.d0,dlh,cknn,cv1,4)
			//nord=4
			//  res_cv = 0.5d0*dl*( 1.d0+DLOG(16.d0*an/dl) )/(pi*an)	! dla ifun=1
			//  res_cs = dl*( 1.d0+DLOG(16.d0*an/dl) )/(pi*an)	        ! dla ifun=3		    
			// --------------- cv1 ------------------
			//dcqg_a = 0.5d0*(0.0d0+dlh)
			//dcqg_b = dlh-0.0d0
			double dcqg_a = (double)0.5e0*dlh;
			double dcqg_b = dlh;
			//
			double dcqg_c = (double)0.43056815579702629e0 * dcqg_b;
			//y=.17392742256872693d0*(f(a+c)+f(a-c))
			//s=a+c
			double cknn_s = dcqg_a + dcqg_c;
			// ---- cknn ----
			double cknn_r = (double)0.5e0*dl - cknn_s;
			double cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			///cknn_cbnd1 = ( dcos(ka*cknn_rr) - cj*dsin(ka*cknn_rr) - cone )/cknn_rr
			///cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cuDoubleComplex cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			double cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			double cknn_res = delick_gpu(cknn_bet);
			double cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))        !1 segment n-tej f.bazowej
			cuDoubleComplex cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//cv = Z_MUL( Z_MAKE(cknn_s/dl,(double)0.0e0),cv );
			cuDoubleComplex cvv1a = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)          ! pot. scalarny 
			//cv = Z_ADD( cknn_cbnd1,Z_MAKE(cknn_bnd2,(double)0.0e0) );
			cuDoubleComplex csv1a = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			// ---- cknn ----
			//s=a-c
			cknn_s = dcqg_a - dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			//cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//!cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))        !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			cv = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			cvv1a = Z_MUL(Z_MAKE((double)0.17392742256872693e0, (double)0.0e0), Z_ADD(cv, cvv1a));
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)            ! pot. scalarny 
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//csv1a = .17392742256872693d0*( cv + csv1a )
			csv1a = Z_MUL(Z_MAKE((double)0.17392742256872693e0, (double)0.0e0), Z_ADD(cv, csv1a));
			// ---- cknn ---- 

			dcqg_c = (double)0.16999052179242813e0 * dcqg_b;
			//y=b*(y+.32607257743127307d0*(f(a+c)+f(a-c))) 
			//s=a+c
			cknn_s = dcqg_a + dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			//cknn_cbnd1 = ( dcos(ka*cknn_rr) - cj*dsin(ka*cknn_rr) - cone )/cknn_rr
			//cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))         !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//cv = Z_MUL( Z_MAKE(cknn_s/dl,(double)0.0e0),cv );
			cuDoubleComplex cvv1b = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)               ! pot. scalarny 
			//cv = Z_ADD( cknn_cbnd1,Z_MAKE(cknn_bnd2,(double)0.0e0) );
			cuDoubleComplex csv1b = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			// ---- cknn ----
			//s=a-c
			cknn_s = dcqg_a - dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			//cknn_cbnd1 = ( dcos(ka*cknn_rr) - cj*dsin(ka*cknn_rr) - cone )/cknn_rr
			//cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))            !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			cv = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			cvv1b = Z_MUL(Z_MAKE((double)0.32607257743127307e0, (double)0.0e0), Z_ADD(cv, cvv1b));
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)                      ! pot. scalarny
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//csv1b = .32607257743127307d0*( cv + csv1b ) 
			csv1b = Z_MUL(Z_MAKE((double)0.32607257743127307e0, (double)0.0e0), Z_ADD(cv, csv1b));
			// ---- cknn ---- 
			//cv1v = dcqg_b*( cvv1a+cvv1b )
			//cv1s = dcqg_b*( csv1a+csv1b ) 
			cuDoubleComplex cv1v = Z_MUL(Z_MAKE(dcqg_b, (double)0.0e0), Z_ADD(cvv1a, cvv1b));
			cuDoubleComplex cv1s = Z_MUL(Z_MAKE(dcqg_b, (double)0.0e0), Z_ADD(csv1a, csv1b));
			// --------------- cv1 ------------------


			//call dcqg(dlh,dl,cknn,cv2,4)	      
			//nord=4
			// --------------- cv2 ------------------
			dcqg_a = (double)0.5e0*(dlh + dl);
			dcqg_b = dl - dlh;
			//
			dcqg_c = (double)0.43056815579702629e0 * dcqg_b;
			//y=.17392742256872693d0*(f(a+c)+f(a-c))
			//s=a+c
			cknn_s = dcqg_a + dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			///cknn_cbnd1 = ( dcos(ka*cknn_rr) - cj*dsin(ka*cknn_rr) - cone )/cknn_rr
			///cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))        !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			cv = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			cvv1a = cv;
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)          ! pot. scalarny 
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			csv1a = cv;
			// ---- cknn ----
			//s=a-c
			cknn_s = dcqg_a - dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			//cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//!cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))        !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			cv = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			cvv1a = Z_MUL(Z_MAKE((double)0.17392742256872693e0, (double)0.0e0), Z_ADD(cv, cvv1a));
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)            ! pot. scalarny 
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//csv1a = .17392742256872693d0*( cv + csv1a )
			csv1a = Z_MUL(Z_MAKE((double)0.17392742256872693e0, (double)0.0e0), Z_ADD(cv, csv1a));
			// ---- cknn ---- 

			dcqg_c = (double)0.16999052179242813e0 * dcqg_b;
			//y=b*(y+.32607257743127307d0*(f(a+c)+f(a-c))) 
			//s=a+c
			cknn_s = dcqg_a + dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			//cknn_cbnd1 = ( dcos(ka*cknn_rr) - cj*dsin(ka*cknn_rr) - cone )/cknn_rr
			//cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))         !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//cv = Z_MUL( Z_MAKE(cknn_s/dl,(double)0.0e0),cv );
			cvv1b = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)               ! pot. scalarny 
			//cv = Z_ADD( cknn_cbnd1,Z_MAKE(cknn_bnd2,(double)0.0e0) );
			csv1b = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			// ---- cknn ----
			//s=a-c
			cknn_s = dcqg_a - dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			//cknn_cbnd1 = ( dcos(ka*cknn_rr) - cj*dsin(ka*cknn_rr) - cone )/cknn_rr
			//cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))            !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			cv = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			cvv1b = Z_MUL(Z_MAKE((double)0.32607257743127307e0, (double)0.0e0), Z_ADD(cv, cvv1b));
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)                      ! pot. scalarny
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//csv1b = .32607257743127307d0*( cv + csv1b ) 
			csv1b = Z_MUL(Z_MAKE((double)0.32607257743127307e0, (double)0.0e0), Z_ADD(cv, csv1b));
			// ---- cknn ---- 
			//cv2v = dcqg_b*( cvv1a+cvv1b )
			//cv2s = dcqg_b*( csv1a+csv1b  
			cuDoubleComplex cv2v = Z_MUL(Z_MAKE(dcqg_b, (double)0.0e0), Z_ADD(cvv1a, cvv1b));
			cuDoubleComplex cv2s = Z_MUL(Z_MAKE(dcqg_b, (double)0.0e0), Z_ADD(csv1a, csv1b));
			// --------------- cv2 ------------------

			//res_cv = 0.5d0*dl*( 1.d0+DLOG(16.d0*an/dl) )/(pi*an)	! dla ifun=1 i ifun=2
			//res_cs = dl*( 1.d0+DLOG(16.d0*an/dl) )/(pi*an)	        ! dla ifun=3
			double res_cs = dl*((double)1.0e0 + log((double)16.0e0*an / dl)) / (pi*an);
			double res_cv = (double)0.5e0*res_cs;
			//cvval2 = cv1v + cv2v + DCMPLX(res_cv,0.d0)     
			//csval2 = cv1s + cv2s + DCMPLX(res_cs,0.d0)
			cvval2 = Z_ADD(cv2v, Z_MAKE(res_cv, (double)0.0e0));
			cvval2 = Z_ADD(cv1v, cvval2);
			csval2 = Z_ADD(cv2s, Z_MAKE(res_cs, (double)0.0e0));
			csval2 = Z_ADD(cv1s, csval2);
		}
		else {
			//call dcqg(0.d0,dl,ckmn,cv,2)
			//cvval2 = dcmplx(0.0d0,0.0d0)
			//csval2 = dcmplx(0.0d0,0.0d0)
			cvval2 = Z_MAKE((double)0.0e0, (double)0.0e0);
			csval2 = Z_MAKE((double)0.0e0, (double)0.0e0);
		}
		// ------------ Zamiast prcodeury cwire1 -------------------

		//tvecl2 = dc%x*lmc%x+dc%y*lmc%y+dc%z*lmc%z      ! obliczane !lm+ * 1n- dla k=1,	obliczane  !lm- * 1n- dla k=2 
		//!cvval = dcmplx(tvecl1,0.d0)*cvval1+dcmplx(tvecl2,0.d0)*cvval2   !A(rm+)	dla k=1, A(rm-) dla k=2
		double tvecl2 = dc.x*lmc.x + dc.y*lmc.y + dc.z*lmc.z;
		//csval2 = csval2/dl
		csval2 = Z_DIV(csval2, Z_MAKE(dl, (double)0.0e0));

		//cvval = dcmplx(tvecl2,0.d0)*cvval2
		cuDoubleComplex cvval = Z_MUL(cvval2, Z_MAKE(tvecl2, (double)0.0e0));
		//csval = -csval2
		cuDoubleComplex csval = Z_MUL(Z_MAKE((double)-1.0e0, (double)0.0e0), csval2);

		//cz_ij = cz_ij+ka*ka*cvval + znak * csval
		cuDoubleComplex cz_ij = Z_MUL(Z_MAKE((double)znak, (double)0.0e0), csval);
		cz_ij = Z_ADD(Z_MUL(Z_MAKE(ka*ka, (double)0.0e0), cvval), cz_ij);

		//cz_ij = cz_ij*cj*eta0/(4.0e0*pi*ka)
		cz_ij = Z_MUL(cz_ij, Z_MAKE(eta0 / ((double)4.0e0*pi*ka), (double)0.0e0));
		cz_ij = Z_MUL(cz_ij, Z_ONE_IMG);
		//cuDoubleComplex cz_ij1 = cz[ i+j*n ];

		//cz_d[i_m+j_n*n] = Z_MAKE(10.0f,0.0f);
		//cz_d[i_m+j_n*n] = Z_MUL ( cz_d[i_m+j_n*n] , (Z_MAKE(10.0f,0.0f));
		//cz_d[i_m+j_n*n] = Z_MAKE((double)i_m,(double)j_n);
		//cz_d[i_m+j_n*n]=make_cuDoubleComplex( (i_m),(j_n) );
		//cz_d[ i+j*n ] =  Z_MAKE(float(i*n+j),0.0f) ;
		//cz[ i+j*n ] = Z_MAKE(x[mp-1],z[mp-1]);
		//cz[ i+j*n ] = Z_MAKE(dc.z,lmc.z);
		//cz[ i+j*n ] = cvval;
		//cz[ i+j*n ] = Z_MAKE(tvecl2,(double)0.0e0);
		cz[i + (j + kk)*n] = Z_ADD(cz[i + (j + kk)*n], cz_ij);
		//cz[ i+j*n ] = cz_ij;
		//cz[ i+j*n ] = Z_ADD( cz[ i+j*n ],Z_MAKE((double)10.0f,(double)0.0f) );
		//cz[ i+j*n ] = Z_ADD(cz_ij1, cz_ij );
		//cz[ i+j*n ] =  Z_ADD(c1.z,Z_MAKE(float(i*n+j),0.0f)) ;
	}
}

//----------------------------------------------------------------------

__global__ void wim_gpu_cknn_cwire2_im1(cuDoubleComplex *cz, float *x, float *y, float *z,
	float *rad, int *ipc, int *ips, int *ise, int nsmax,
	int n, double ka, double eta0,
	float sgnx, float sgny, float sgnz, float sgnenv,
	int realSubSliceSize, int k, int kk)
{

	double pi = acos(double(-1.e0));
	cuDoubleComplex cvval2, csval2;

	// Block index
	int bx = blockIdx.x;
	int by = blockIdx.y;

	// Thread index
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	// Accumulate row i of A and column j of B
	int i = tx + bx * blockDim.x;
	int j = ty + by * blockDim.y;

	//c16vec c1 = make_c16vec(Z_ONE,Z_MAKE(0.0f, 1.0f),Z_MUL(Z_MAKE(2.0f, 0.0f),Z_ONE));

	if (i < n  && j < realSubSliceSize) {
		int mp = abs(ipc[i]);
		int ms = ips[i];
		int mlr = ise[ms - 1];
		//double am = rad[ms-1];
		double3 rmp = make_double3((double)x[mp - 1], (double)y[mp - 1], (double)z[mp - 1]);
		double3 rmlr = make_double3((double)x[mlr - 1], (double)y[mlr - 1], (double)z[mlr - 1]);
		double3 rmc = make_double3((double)0.5e0*(rmlr.x + rmp.x), (double)0.5e0*(rmlr.y + rmp.y), (double)0.5e0*(rmlr.z + rmp.z));
		float znak = -1.0f;
		double3 lmc = make_double3((double)znak*(rmc.x - rmp.x), (double)znak*(rmc.y - rmp.y), (double)znak*(rmc.z - rmp.z));

		// ------------ Zamiast prcodeury cwire1 -------------------
		//ifun=1 i ifun=2=> res_cv
		//fun=3 => res_cs		  
		//call cwire1_gpu(ms,n,cvval1)             ! obliczanie calki na Sn+  pot.wek.
		//call cwire1(ms,n,csval1)

		//int np = abs(ipc[k-1+kk-1+j]);
		//int ns = ips[k-1+kk-1+j];
		int np = abs(ipc[k + kk + j]);
		int ns = ips[k + kk + j + nsmax];
		int nlr = ise[ns - 1 + nsmax];
		double an = rad[ns - 1];
		double3 rnr = make_double3((double)sgnx*(double)x[nlr - 1], (double)sgny*(double)y[nlr - 1], (double)sgnz*(double)z[nlr - 1]);
		double3 rn = make_double3((double)sgnx*(double)x[np - 1], (double)sgny*(double)y[np - 1], (double)sgnz*(double)z[np - 1]);
		double3 rd = make_double3(rnr.x - rn.x, rnr.y - rn.y, rnr.z - rn.z);
		double dl = sqrt(rd.x*rd.x + rd.y*rd.y + rd.z*rd.z);
		double dlh = (double)0.5e0*dl;
		double3 dc = make_double3(rd.x / dl, rd.y / dl, rd.z / dl);

		int ise_ms1 = ise[ms - 1];
		int ise_ms2 = ise[ms - 1 + nsmax];
		int ise_ns1 = ise[ns - 1];
		int ise_ns2 = ise[ns - 1 + nsmax];
		//float z_ise_ms1 = z[ise[ms-1]];
		//float z_ise_ms2 = z[ise[ms-1+nsmax]];
		float z_ise_ms1 = z[ise[ms - 1] - 1];
		float z_ise_ms2 = z[ise[ms - 1 + nsmax] - 1];

		if (war_gpu_pec_im(ise_ms1, ise_ms2, ise_ns1, ise_ns2, z_ise_ms1, z_ise_ms2)) {
			//call dcqg(0.d0,dlh,cknn,cv1,4)
			//nord=4
			//  res_cv = 0.5d0*dl*( 1.d0+DLOG(16.d0*an/dl) )/(pi*an)	! dla ifun=1
			//  res_cs = dl*( 1.d0+DLOG(16.d0*an/dl) )/(pi*an)	        ! dla ifun=3		    
			// --------------- cv1 ------------------
			//dcqg_a = 0.5d0*(0.0d0+dlh)
			//dcqg_b = dlh-0.0d0
			double dcqg_a = (double)0.5e0*dlh;
			double dcqg_b = dlh;
			//
			double dcqg_c = (double)0.43056815579702629e0 * dcqg_b;
			//y=.17392742256872693d0*(f(a+c)+f(a-c))
			//s=a+c
			double cknn_s = dcqg_a + dcqg_c;
			// ---- cknn ----
			double cknn_r = (double)0.5e0*dl - cknn_s;
			double cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			///cknn_cbnd1 = ( dcos(ka*cknn_rr) - cj*dsin(ka*cknn_rr) - cone )/cknn_rr
			///cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cuDoubleComplex cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			double cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			double cknn_res = delick_gpu(cknn_bet);
			double cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))        !1 segment n-tej f.bazowej
			cuDoubleComplex cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//cv = Z_MUL( Z_MAKE(cknn_s/dl,(double)0.0e0),cv );
			cuDoubleComplex cvv1a = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)          ! pot. scalarny 
			//cv = Z_ADD( cknn_cbnd1,Z_MAKE(cknn_bnd2,(double)0.0e0) );
			cuDoubleComplex csv1a = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			// ---- cknn ----
			//s=a-c
			cknn_s = dcqg_a - dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			//cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//!cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))        !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			cv = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			cvv1a = Z_MUL(Z_MAKE((double)0.17392742256872693e0, (double)0.0e0), Z_ADD(cv, cvv1a));
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)            ! pot. scalarny 
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//csv1a = .17392742256872693d0*( cv + csv1a )
			csv1a = Z_MUL(Z_MAKE((double)0.17392742256872693e0, (double)0.0e0), Z_ADD(cv, csv1a));
			// ---- cknn ---- 

			dcqg_c = (double)0.16999052179242813e0 * dcqg_b;
			//y=b*(y+.32607257743127307d0*(f(a+c)+f(a-c))) 
			//s=a+c
			cknn_s = dcqg_a + dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			//cknn_cbnd1 = ( dcos(ka*cknn_rr) - cj*dsin(ka*cknn_rr) - cone )/cknn_rr
			//cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))         !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//cv = Z_MUL( Z_MAKE(cknn_s/dl,(double)0.0e0),cv );
			cuDoubleComplex cvv1b = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)               ! pot. scalarny 
			//cv = Z_ADD( cknn_cbnd1,Z_MAKE(cknn_bnd2,(double)0.0e0) );
			cuDoubleComplex csv1b = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			// ---- cknn ----
			//s=a-c
			cknn_s = dcqg_a - dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			//cknn_cbnd1 = ( dcos(ka*cknn_rr) - cj*dsin(ka*cknn_rr) - cone )/cknn_rr
			//cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))            !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			cv = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			cvv1b = Z_MUL(Z_MAKE((double)0.32607257743127307e0, (double)0.0e0), Z_ADD(cv, cvv1b));
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)                      ! pot. scalarny
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//csv1b = .32607257743127307d0*( cv + csv1b ) 
			csv1b = Z_MUL(Z_MAKE((double)0.32607257743127307e0, (double)0.0e0), Z_ADD(cv, csv1b));
			// ---- cknn ---- 
			//cv1v = dcqg_b*( cvv1a+cvv1b )
			//cv1s = dcqg_b*( csv1a+csv1b ) 
			cuDoubleComplex cv1v = Z_MUL(Z_MAKE(dcqg_b, (double)0.0e0), Z_ADD(cvv1a, cvv1b));
			cuDoubleComplex cv1s = Z_MUL(Z_MAKE(dcqg_b, (double)0.0e0), Z_ADD(csv1a, csv1b));
			// --------------- cv1 ------------------


			//call dcqg(dlh,dl,cknn,cv2,4)	      
			//nord=4
			// --------------- cv2 ------------------
			dcqg_a = (double)0.5e0*(dlh + dl);
			dcqg_b = dl - dlh;
			//
			dcqg_c = (double)0.43056815579702629e0 * dcqg_b;
			//y=.17392742256872693d0*(f(a+c)+f(a-c))
			//s=a+c
			cknn_s = dcqg_a + dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			///cknn_cbnd1 = ( dcos(ka*cknn_rr) - cj*dsin(ka*cknn_rr) - cone )/cknn_rr
			///cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))        !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			cv = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			cvv1a = cv;
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)          ! pot. scalarny 
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			csv1a = cv;
			// ---- cknn ----
			//s=a-c
			cknn_s = dcqg_a - dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			//cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//!cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))        !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			cv = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			cvv1a = Z_MUL(Z_MAKE((double)0.17392742256872693e0, (double)0.0e0), Z_ADD(cv, cvv1a));
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)            ! pot. scalarny 
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//csv1a = .17392742256872693d0*( cv + csv1a )
			csv1a = Z_MUL(Z_MAKE((double)0.17392742256872693e0, (double)0.0e0), Z_ADD(cv, csv1a));
			// ---- cknn ---- 

			dcqg_c = (double)0.16999052179242813e0 * dcqg_b;
			//y=b*(y+.32607257743127307d0*(f(a+c)+f(a-c))) 
			//s=a+c
			cknn_s = dcqg_a + dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			//cknn_cbnd1 = ( dcos(ka*cknn_rr) - cj*dsin(ka*cknn_rr) - cone )/cknn_rr
			//cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))         !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//cv = Z_MUL( Z_MAKE(cknn_s/dl,(double)0.0e0),cv );
			cvv1b = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)               ! pot. scalarny 
			//cv = Z_ADD( cknn_cbnd1,Z_MAKE(cknn_bnd2,(double)0.0e0) );
			csv1b = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			// ---- cknn ----
			//s=a-c
			cknn_s = dcqg_a - dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			//cknn_cbnd1 = ( dcos(ka*cknn_rr) - cj*dsin(ka*cknn_rr) - cone )/cknn_rr
			//cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))            !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			cv = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			cvv1b = Z_MUL(Z_MAKE((double)0.32607257743127307e0, (double)0.0e0), Z_ADD(cv, cvv1b));
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)                      ! pot. scalarny
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//csv1b = .32607257743127307d0*( cv + csv1b ) 
			csv1b = Z_MUL(Z_MAKE((double)0.32607257743127307e0, (double)0.0e0), Z_ADD(cv, csv1b));
			// ---- cknn ---- 
			//cv2v = dcqg_b*( cvv1a+cvv1b )
			//cv2s = dcqg_b*( csv1a+csv1b  
			cuDoubleComplex cv2v = Z_MUL(Z_MAKE(dcqg_b, (double)0.0e0), Z_ADD(cvv1a, cvv1b));
			cuDoubleComplex cv2s = Z_MUL(Z_MAKE(dcqg_b, (double)0.0e0), Z_ADD(csv1a, csv1b));
			// --------------- cv2 ------------------

			//res_cv = 0.5d0*dl*( 1.d0+DLOG(16.d0*an/dl) )/(pi*an)	! dla ifun=1 i ifun=2
			//res_cs = dl*( 1.d0+DLOG(16.d0*an/dl) )/(pi*an)	        ! dla ifun=3
			double res_cs = dl*((double)1.0e0 + log((double)16.0e0*an / dl)) / (pi*an);
			double res_cv = (double)0.5e0*res_cs;
			//cvval2 = cv1v + cv2v + DCMPLX(res_cv,0.d0)     
			//csval2 = cv1s + cv2s + DCMPLX(res_cs,0.d0)
			cvval2 = Z_ADD(cv2v, Z_MAKE(res_cv, (double)0.0e0));
			cvval2 = Z_ADD(cv1v, cvval2);
			csval2 = Z_ADD(cv2s, Z_MAKE(res_cs, (double)0.0e0));
			csval2 = Z_ADD(cv1s, csval2);
		}
		else {
			//call dcqg(0.d0,dl,ckmn,cv,2)
			//cvval2 = dcmplx(0.0d0,0.0d0)
			//csval2 = dcmplx(0.0d0,0.0d0)
			cvval2 = Z_MAKE((double)0.0e0, (double)0.0e0);
			csval2 = Z_MAKE((double)0.0e0, (double)0.0e0);
		}
		// ------------ Zamiast prcodeury cwire2 -------------------

		//tvecl2 = dc%x*lmc%x+dc%y*lmc%y+dc%z*lmc%z      ! obliczane !lm+ * 1n- dla k=1,	obliczane  !lm- * 1n- dla k=2 
		//!cvval = dcmplx(tvecl1,0.d0)*cvval1+dcmplx(tvecl2,0.d0)*cvval2   !A(rm+)	dla k=1, A(rm-) dla k=2
		double tvecl2 = dc.x*lmc.x + dc.y*lmc.y + dc.z*lmc.z;
		//csval2 = csval2/dl
		csval2 = Z_DIV(csval2, Z_MAKE(dl, (double)0.0e0));

		//cvval = sgnenv*dcmplx(tvecl2,0.d0)*cvval2
		cuDoubleComplex cvval = Z_MUL(cvval2, Z_MAKE((double)sgnenv*tvecl2, (double)0.0e0));
		//csval = -sgnenv*csval2
		cuDoubleComplex csval = Z_MUL(Z_MAKE((double)-1.0e0*(double)sgnenv, (double)0.0e0), csval2);

		//cz_ij = cz_ij+ka*ka*cvval + znak * csval
		cuDoubleComplex cz_ij = Z_MUL(Z_MAKE((double)znak, (double)0.0e0), csval);
		cz_ij = Z_ADD(Z_MUL(Z_MAKE(ka*ka, (double)0.0e0), cvval), cz_ij);

		//cz_ij = cz_ij*cj*eta0/(4.0e0*pi*ka)
		cz_ij = Z_MUL(cz_ij, Z_MAKE(eta0 / ((double)4.0e0*pi*ka), (double)0.0e0));
		cz_ij = Z_MUL(cz_ij, Z_ONE_IMG);


		//cz_d[i_m+j_n*n] = Z_MAKE(10.0f,0.0f);
		//cz_d[i_m+j_n*n] = Z_MUL ( cz_d[i_m+j_n*n] , (Z_MAKE(10.0f,0.0f));
		//cz_d[i_m+j_n*n] = Z_MAKE((double)i_m,(double)j_n);
		//cz_d[i_m+j_n*n]=make_cuDoubleComplex( (i_m),(j_n) );
		//cz_d[ i+j*n ] =  Z_MAKE(float(i*n+j),0.0f) ;
		//cz[ i+j*n ] = Z_MAKE(x[mp-1],z[mp-1]);
		//cz[ i+j*n ] = Z_MAKE(dc.x,dc.z);
		//cz[ i+j*n ] = cknn_cbnd1;
		//cz[ i+j*n ] = Z_MAKE(cknn_res,cknn_bnd2);
		cz[i + (j + kk)*n] = Z_ADD(cz[i + (j + kk)*n], cz_ij);
		//cz[ i+j*n ] =  Z_ADD(c1.z,Z_MAKE(float(i*n+j),0.0f)) ;
	}
}

//----------------------------------------------------------------------

__global__ void wim_gpu_cknn_cwire2_im2(cuDoubleComplex *cz, float *x, float *y, float *z,
	float *rad, int *ipc, int *ips, int *ise, int nsmax,
	int n, double ka, double eta0,
	float sgnx, float sgny, float sgnz, float sgnenv,
	int realSubSliceSize, int k, int kk)
{

	double pi = acos(double(-1.e0));
	cuDoubleComplex cvval2, csval2;

	// Block index
	int bx = blockIdx.x;
	int by = blockIdx.y;

	// Thread index
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	// Accumulate row i of A and column j of B
	int i = tx + bx * blockDim.x;
	int j = ty + by * blockDim.y;

	//c16vec c1 = make_c16vec(Z_ONE,Z_MAKE(0.0f, 1.0f),Z_MUL(Z_MAKE(2.0f, 0.0f),Z_ONE));

	if (i < n  && j < realSubSliceSize) {
		int mp = abs(ipc[i]);
		//int ms = ips[i];
		//int mlr = ise[ms-1];
		int ms = ips[i + nsmax];
		int mlr = ise[ms - 1 + nsmax];
		//double am = rad[ms-1];
		double3 rmp = make_double3((double)x[mp - 1], (double)y[mp - 1], (double)z[mp - 1]);
		double3 rmlr = make_double3((double)x[mlr - 1], (double)y[mlr - 1], (double)z[mlr - 1]);
		double3 rmc = make_double3((double)0.5e0*(rmlr.x + rmp.x), (double)0.5e0*(rmlr.y + rmp.y), (double)0.5e0*(rmlr.z + rmp.z));
		float znak = 1.0f;
		double3 lmc = make_double3((double)znak*(rmc.x - rmp.x), (double)znak*(rmc.y - rmp.y), (double)znak*(rmc.z - rmp.z));

		// ------------ Zamiast prcodeury cwire1 -------------------
		//ifun=1 i ifun=2=> res_cv
		//fun=3 => res_cs		  
		//call cwire1_gpu(ms,n,cvval1)             ! obliczanie calki na Sn+  pot.wek.
		//call cwire1(ms,n,csval1)

		//int np = abs(ipc[k-1+kk-1+j]);
		//int ns = ips[k-1+kk-1+j];
		int np = abs(ipc[k + kk + j]);
		int ns = ips[k + kk + j + nsmax];
		int nlr = ise[ns - 1 + nsmax];
		double an = rad[ns - 1];
		double3 rnr = make_double3((double)sgnx*(double)x[nlr - 1], (double)sgny*(double)y[nlr - 1], (double)sgnz*(double)z[nlr - 1]);
		double3 rn = make_double3((double)sgnx*(double)x[np - 1], (double)sgny*(double)y[np - 1], (double)sgnz*(double)z[np - 1]);
		double3 rd = make_double3(rnr.x - rn.x, rnr.y - rn.y, rnr.z - rn.z);
		double dl = sqrt(rd.x*rd.x + rd.y*rd.y + rd.z*rd.z);
		double dlh = (double)0.5e0*dl;
		double3 dc = make_double3(rd.x / dl, rd.y / dl, rd.z / dl);

		int ise_ms1 = ise[ms - 1];
		int ise_ms2 = ise[ms - 1 + nsmax];
		int ise_ns1 = ise[ns - 1];
		int ise_ns2 = ise[ns - 1 + nsmax];
		//float z_ise_ms1 = z[ise[ms-1]];
		//float z_ise_ms2 = z[ise[ms-1+nsmax]];
		float z_ise_ms1 = z[ise[ms - 1] - 1];
		float z_ise_ms2 = z[ise[ms - 1 + nsmax] - 1];

		if (war_gpu_pec_im(ise_ms1, ise_ms2, ise_ns1, ise_ns2, z_ise_ms1, z_ise_ms2)) {
			//call dcqg(0.d0,dlh,cknn,cv1,4)
			//nord=4
			//  res_cv = 0.5d0*dl*( 1.d0+DLOG(16.d0*an/dl) )/(pi*an)	! dla ifun=1
			//  res_cs = dl*( 1.d0+DLOG(16.d0*an/dl) )/(pi*an)	        ! dla ifun=3		    
			// --------------- cv1 ------------------
			//dcqg_a = 0.5d0*(0.0d0+dlh)
			//dcqg_b = dlh-0.0d0
			double dcqg_a = (double)0.5e0*dlh;
			double dcqg_b = dlh;
			//
			double dcqg_c = (double)0.43056815579702629e0 * dcqg_b;
			//y=.17392742256872693d0*(f(a+c)+f(a-c))
			//s=a+c
			double cknn_s = dcqg_a + dcqg_c;
			// ---- cknn ----
			double cknn_r = (double)0.5e0*dl - cknn_s;
			double cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			///cknn_cbnd1 = ( dcos(ka*cknn_rr) - cj*dsin(ka*cknn_rr) - cone )/cknn_rr
			///cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cuDoubleComplex cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			double cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			double cknn_res = delick_gpu(cknn_bet);
			double cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))        !1 segment n-tej f.bazowej
			cuDoubleComplex cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//cv = Z_MUL( Z_MAKE(cknn_s/dl,(double)0.0e0),cv );
			cuDoubleComplex cvv1a = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)          ! pot. scalarny 
			//cv = Z_ADD( cknn_cbnd1,Z_MAKE(cknn_bnd2,(double)0.0e0) );
			cuDoubleComplex csv1a = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			// ---- cknn ----
			//s=a-c
			cknn_s = dcqg_a - dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			//cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//!cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))        !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			cv = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			cvv1a = Z_MUL(Z_MAKE((double)0.17392742256872693e0, (double)0.0e0), Z_ADD(cv, cvv1a));
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)            ! pot. scalarny 
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//csv1a = .17392742256872693d0*( cv + csv1a )
			csv1a = Z_MUL(Z_MAKE((double)0.17392742256872693e0, (double)0.0e0), Z_ADD(cv, csv1a));
			// ---- cknn ---- 

			dcqg_c = (double)0.16999052179242813e0 * dcqg_b;
			//y=b*(y+.32607257743127307d0*(f(a+c)+f(a-c))) 
			//s=a+c
			cknn_s = dcqg_a + dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			//cknn_cbnd1 = ( dcos(ka*cknn_rr) - cj*dsin(ka*cknn_rr) - cone )/cknn_rr
			//cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))         !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//cv = Z_MUL( Z_MAKE(cknn_s/dl,(double)0.0e0),cv );
			cuDoubleComplex cvv1b = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)               ! pot. scalarny 
			//cv = Z_ADD( cknn_cbnd1,Z_MAKE(cknn_bnd2,(double)0.0e0) );
			cuDoubleComplex csv1b = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			// ---- cknn ----
			//s=a-c
			cknn_s = dcqg_a - dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			//cknn_cbnd1 = ( dcos(ka*cknn_rr) - cj*dsin(ka*cknn_rr) - cone )/cknn_rr
			//cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))            !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			cv = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			cvv1b = Z_MUL(Z_MAKE((double)0.32607257743127307e0, (double)0.0e0), Z_ADD(cv, cvv1b));
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)                      ! pot. scalarny
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//csv1b = .32607257743127307d0*( cv + csv1b ) 
			csv1b = Z_MUL(Z_MAKE((double)0.32607257743127307e0, (double)0.0e0), Z_ADD(cv, csv1b));
			// ---- cknn ---- 
			//cv1v = dcqg_b*( cvv1a+cvv1b )
			//cv1s = dcqg_b*( csv1a+csv1b ) 
			cuDoubleComplex cv1v = Z_MUL(Z_MAKE(dcqg_b, (double)0.0e0), Z_ADD(cvv1a, cvv1b));
			cuDoubleComplex cv1s = Z_MUL(Z_MAKE(dcqg_b, (double)0.0e0), Z_ADD(csv1a, csv1b));
			// --------------- cv1 ------------------


			//call dcqg(dlh,dl,cknn,cv2,4)	      
			//nord=4
			// --------------- cv2 ------------------
			dcqg_a = (double)0.5e0*(dlh + dl);
			dcqg_b = dl - dlh;
			//
			dcqg_c = (double)0.43056815579702629e0 * dcqg_b;
			//y=.17392742256872693d0*(f(a+c)+f(a-c))
			//s=a+c
			cknn_s = dcqg_a + dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			///cknn_cbnd1 = ( dcos(ka*cknn_rr) - cj*dsin(ka*cknn_rr) - cone )/cknn_rr
			///cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))        !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			cv = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			cvv1a = cv;
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)          ! pot. scalarny 
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			csv1a = cv;
			// ---- cknn ----
			//s=a-c
			cknn_s = dcqg_a - dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			//cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//!cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))        !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			cv = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			cvv1a = Z_MUL(Z_MAKE((double)0.17392742256872693e0, (double)0.0e0), Z_ADD(cv, cvv1a));
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)            ! pot. scalarny 
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//csv1a = .17392742256872693d0*( cv + csv1a )
			csv1a = Z_MUL(Z_MAKE((double)0.17392742256872693e0, (double)0.0e0), Z_ADD(cv, csv1a));
			// ---- cknn ---- 

			dcqg_c = (double)0.16999052179242813e0 * dcqg_b;
			//y=b*(y+.32607257743127307d0*(f(a+c)+f(a-c))) 
			//s=a+c
			cknn_s = dcqg_a + dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			//cknn_cbnd1 = ( dcos(ka*cknn_rr) - cj*dsin(ka*cknn_rr) - cone )/cknn_rr
			//cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))         !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//cv = Z_MUL( Z_MAKE(cknn_s/dl,(double)0.0e0),cv );
			cvv1b = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)               ! pot. scalarny 
			//cv = Z_ADD( cknn_cbnd1,Z_MAKE(cknn_bnd2,(double)0.0e0) );
			csv1b = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			// ---- cknn ----
			//s=a-c
			cknn_s = dcqg_a - dcqg_c;
			// ---- cknn ----
			cknn_r = (double)0.5e0*dl - cknn_s;
			cknn_rr = sqrt(an*an + cknn_r*cknn_r);
			//cknn_cbnd1 = ( dcos(ka*cknn_rr) - cj*dsin(ka*cknn_rr) - cone )/cknn_rr
			//cknn_cbnd1 = dcmplx( cos(ka*cknn_rr)-1.0d0, - sin(ka*cknn_rr) )/cknn_rr
			cknn_cbnd1 = Z_MAKE((cos(ka*cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka*cknn_rr) / cknn_rr);
			cknn_bet = (an + an) / sqrt((double)4.0e0*an*an + cknn_r*cknn_r);
			//call delick_gpu(cknn_bet,cknn_res)
			cknn_res = delick_gpu(cknn_bet);
			cknn_bnd2 = (cknn_bet*cknn_res + log(abs(cknn_r / ((double)8.e0*an)))) / (pi*an);
			//ifun=1
			//cv = dcmplx(cknn_s/dl,0.d0)*(cknn_cbnd1+DCMPLX(cknn_bnd2,0.d0))            !1 segment n-tej f.bazowej
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			cv = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
			cvv1b = Z_MUL(Z_MAKE((double)0.32607257743127307e0, (double)0.0e0), Z_ADD(cv, cvv1b));
			//ifun=3
			//cv = cknn_cbnd1+dcmplx(cknn_bnd2,0.d0)                      ! pot. scalarny
			cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
			//csv1b = .32607257743127307d0*( cv + csv1b ) 
			csv1b = Z_MUL(Z_MAKE((double)0.32607257743127307e0, (double)0.0e0), Z_ADD(cv, csv1b));
			// ---- cknn ---- 
			//cv2v = dcqg_b*( cvv1a+cvv1b )
			//cv2s = dcqg_b*( csv1a+csv1b  
			cuDoubleComplex cv2v = Z_MUL(Z_MAKE(dcqg_b, (double)0.0e0), Z_ADD(cvv1a, cvv1b));
			cuDoubleComplex cv2s = Z_MUL(Z_MAKE(dcqg_b, (double)0.0e0), Z_ADD(csv1a, csv1b));
			// --------------- cv2 ------------------

			//res_cv = 0.5d0*dl*( 1.d0+DLOG(16.d0*an/dl) )/(pi*an)	! dla ifun=1 i ifun=2
			//res_cs = dl*( 1.d0+DLOG(16.d0*an/dl) )/(pi*an)	        ! dla ifun=3
			double res_cs = dl*((double)1.0e0 + log((double)16.0e0*an / dl)) / (pi*an);
			double res_cv = (double)0.5e0*res_cs;
			//cvval2 = cv1v + cv2v + DCMPLX(res_cv,0.d0)     
			//csval2 = cv1s + cv2s + DCMPLX(res_cs,0.d0)
			cvval2 = Z_ADD(cv2v, Z_MAKE(res_cv, (double)0.0e0));
			cvval2 = Z_ADD(cv1v, cvval2);
			csval2 = Z_ADD(cv2s, Z_MAKE(res_cs, (double)0.0e0));
			csval2 = Z_ADD(cv1s, csval2);
		}
		else {
			//call dcqg(0.d0,dl,ckmn,cv,2)
			//cvval2 = dcmplx(0.0d0,0.0d0)
			//csval2 = dcmplx(0.0d0,0.0d0)
			cvval2 = Z_MAKE((double)0.0e0, (double)0.0e0);
			csval2 = Z_MAKE((double)0.0e0, (double)0.0e0);
		}
		// ------------ Zamiast prcodeury cwire1 -------------------

		//tvecl2 = dc%x*lmc%x+dc%y*lmc%y+dc%z*lmc%z      ! obliczane !lm+ * 1n- dla k=1,	obliczane  !lm- * 1n- dla k=2 
		//!cvval = dcmplx(tvecl1,0.d0)*cvval1+dcmplx(tvecl2,0.d0)*cvval2   !A(rm+)	dla k=1, A(rm-) dla k=2
		double tvecl2 = dc.x*lmc.x + dc.y*lmc.y + dc.z*lmc.z;
		//csval2 = csval2/dl
		csval2 = Z_DIV(csval2, Z_MAKE(dl, (double)0.0e0));

		//cvval = cvval = sgnenv*dcmplx(tvecl2,0.d0)*cvval2
		cuDoubleComplex cvval = Z_MUL(cvval2, Z_MAKE((double)sgnenv*tvecl2, (double)0.0e0));
		//csval = -sgnenv*csval2
		cuDoubleComplex csval = Z_MUL(Z_MAKE((double)-1.0e0*(double)sgnenv, (double)0.0e0), csval2);

		//cz_ij = cz_ij+ka*ka*cvval + znak * csval
		cuDoubleComplex cz_ij = Z_MUL(Z_MAKE((double)znak, (double)0.0e0), csval);
		cz_ij = Z_ADD(Z_MUL(Z_MAKE(ka*ka, (double)0.0e0), cvval), cz_ij);

		//cz_ij = cz_ij*cj*eta0/(4.0e0*pi*ka)
		cz_ij = Z_MUL(cz_ij, Z_MAKE(eta0 / ((double)4.0e0*pi*ka), (double)0.0e0));
		cz_ij = Z_MUL(cz_ij, Z_ONE_IMG);
		//cuDoubleComplex cz_ij1 = cz[ i+j*n ];

		//cz_d[i_m+j_n*n] = Z_MAKE(10.0f,0.0f);
		//cz_d[i_m+j_n*n] = Z_MUL ( cz_d[i_m+j_n*n] , (Z_MAKE(10.0f,0.0f));
		//cz_d[i_m+j_n*n] = Z_MAKE((double)i_m,(double)j_n);
		//cz_d[i_m+j_n*n]=make_cuDoubleComplex( (i_m),(j_n) );
		//cz_d[ i+j*n ] =  Z_MAKE(float(i*n+j),0.0f) ;
		//cz[ i+j*n ] = Z_MAKE(x[mp-1],z[mp-1]);
		//cz[ i+j*n ] = Z_MAKE(dc.z,lmc.z);
		//cz[ i+j*n ] = cvval;
		//cz[ i+j*n ] = Z_MAKE(tvecl2,(double)0.0e0);
		cz[i + (j + kk)*n] = Z_ADD(cz[i + (j + kk)*n], cz_ij);
		//cz[ i+j*n ] = cz_ij;
		//cz[ i+j*n ] = Z_ADD( cz[ i+j*n ],Z_MAKE((double)10.0f,(double)0.0f) );
		//cz[ i+j*n ] = Z_ADD(cz_ij1, cz_ij );
		//cz[ i+j*n ] =  Z_ADD(c1.z,Z_MAKE(float(i*n+j),0.0f)) ;
	}
}

//----------------------------------------------------------------------

__global__ void wim_gpu_ckmn_re1(cuDoubleComplex *cz, float *x, float *y, float *z,
	float *rad, int *ipc, int *ips, int *ise, int nsmax,
	int n, double ka, double eta0,
	int realSubSliceSize, int k, int kk)
{

	double pi = acos(double(-1.e0));
	cuDoubleComplex cvval1, csval1, cvval2, csval2;
	double tvecl1, tvecl2;

	// Block index
	int bx = blockIdx.x;
	int by = blockIdx.y;

	// Thread index
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	// Accumulate row i of A and column j of B
	int i = tx + bx * blockDim.x;
	int j = ty + by * blockDim.y;

	//c16vec c1 = make_c16vec(Z_ONE,Z_MAKE(0.0f, 1.0f),Z_MUL(Z_MAKE(2.0f, 0.0f),Z_ONE));

	if (i < n  && j < realSubSliceSize) {
		int mp = abs(ipc[i]);
		int ms = ips[i];
		int mlr = ise[ms - 1];
		//double am = rad[ms-1];
		double3 rmp = make_double3((double)x[mp - 1], (double)y[mp - 1], (double)z[mp - 1]);
		double3 rmlr = make_double3((double)x[mlr - 1], (double)y[mlr - 1], (double)z[mlr - 1]);
		double3 rmc = make_double3((double)0.5e0*(rmlr.x + rmp.x), (double)0.5e0*(rmlr.y + rmp.y), (double)0.5e0*(rmlr.z + rmp.z));
		float znak = -1.0f;
		double3 lmc = make_double3((double)znak*(rmc.x - rmp.x), (double)znak*(rmc.y - rmp.y), (double)znak*(rmc.z - rmp.z));

		// ------------ Zamiast prcodeury cwire1 -------------------
		//ifun=1 i ifun=2=> res_cv
		//fun=3 => res_cs		  
		//call cwire1_gpu(ms,n,cvval1)             ! obliczanie calki na Sn+  pot.wek.
		//call cwire1(ms,n,csval1)

		//int np = abs(ipc[k-1+kk-1+j]);
		//int ns = ips[k-1+kk-1+j];
		int np = abs(ipc[k + kk + j]);
		int ns = ips[k + kk + j];
		double3 rnp = make_double3((double)x[np - 1], (double)y[np - 1], (double)z[np - 1]);
		int nlr = ise[ns - 1];
		double an = rad[ns - 1];
		double3 rn = make_double3((double)x[nlr - 1], (double)y[nlr - 1], (double)z[nlr - 1]);
		double3 rd = make_double3(rnp.x - rn.x, rnp.y - rn.y, rnp.z - rn.z);
		double dl = sqrt(rd.x*rd.x + rd.y*rd.y + rd.z*rd.z);
		//double dlh = (double)0.5e0*dl;
		double3 dc = make_double3(rd.x / dl, rd.y / dl, rd.z / dl);

		int ise_ms1 = ise[ms - 1];
		int ise_ms2 = ise[ms - 1 + nsmax];
		int ise_ns1 = ise[ns - 1];
		int ise_ns2 = ise[ns - 1 + nsmax];


		if (war_gpu_pec_re(ise_ms1, ise_ms2, ise_ns1, ise_ns2)) {
			//call dcqg(0.d0,dlh,cknn,cv1,4)
			//nord=4		    
			// --------------- cv1 ------------------
			cvval1 = Z_MAKE((double)0.0e0, (double)0.0e0);
			csval1 = Z_MAKE((double)0.0e0, (double)0.0e0);
		}
		else {
			//call dcqg(0.d0,dl,ckmn,cv,2)
			//nord=2
			double dcqg_a = (double)0.5e0*dl;
			double dcqg_b = dl;
			double dcqg_c = (double)0.288675134594812882e0 * dcqg_b;
			//ifun=1
			//y=b*0.5d0*(f(a+c)+f(a-c))
			//dla s=a+c
			double ckmn_s = dcqg_a + dcqg_c;
			double3 ckmn_dr = make_double3(rmc.x - rn.x - ckmn_s*dc.x, rmc.y - rn.y - ckmn_s*dc.y, rmc.z - rn.z - ckmn_s*dc.z);
			double ckmn_r = sqrt(an*an + ckmn_dr.x*ckmn_dr.x + ckmn_dr.y*ckmn_dr.y + ckmn_dr.z*ckmn_dr.z);
			//stosujemy wz�r Eulera exp(-ja)=cos(a)-jsin(a) => a=ka*r 
			//c_jkr = cdexp(-cj*ka*ckmn_r)/ckmn_r
			cuDoubleComplex c_jkr = Z_MAKE(cos(ka*ckmn_r) / ckmn_r, -sin(ka*ckmn_r) / ckmn_r);
			//cvval1 = dcmplx(ckmn_s/dl,0.d0)*c_jkr	!ifun=1
			cvval1 = Z_MUL(Z_MAKE(ckmn_s / dl, (double)0.0e0), c_jkr);
			//csval1 = c_jkr	                        !ifun=3
			csval1 = c_jkr;
			//dla s=a-c
			ckmn_s = dcqg_a - dcqg_c;
			ckmn_dr = make_double3(rmc.x - rn.x - ckmn_s*dc.x, rmc.y - rn.y - ckmn_s*dc.y, rmc.z - rn.z - ckmn_s*dc.z);
			ckmn_r = sqrt(an*an + ckmn_dr.x*ckmn_dr.x + ckmn_dr.y*ckmn_dr.y + ckmn_dr.z*ckmn_dr.z);
			//c_jkr = exp(-cj*ka*ckmn_r)/ckmn_r
			//c_jkr = dcmplx( cos(ka*ckmn_r), -sin(ka*ckmn_r) )
			//c_jkr = c_jkr/ckmn_r
			c_jkr = Z_MAKE(cos(ka*ckmn_r) / ckmn_r, -sin(ka*ckmn_r) / ckmn_r);
			//cvval1 =  cvval1 + dcmplx(ckmn_s/dl,0.d0)*c_jkr	!ifun=1	
			//csval1 = csval1 + c_jkr							!ifun=3	
			cvval1 = Z_ADD(Z_MUL(Z_MAKE(ckmn_s / dl, (double)0.0e0), c_jkr), cvval1);
			csval1 = Z_ADD(csval1, c_jkr);

			cvval1 = Z_MUL(Z_MAKE(dcqg_b*(double)0.5e0, (double)0.0e0), cvval1);
			csval1 = Z_MUL(Z_MAKE(dcqg_b*(double)0.5e0, (double)0.0e0), csval1);

			//tvecl1= dc%x*lmc%x+dc%y*lmc%y+dc%z*lmc%z       ! obliczane !lm+ * 1n+ dla k=1, obliczane !lm- * 1n+ dla k=1
			tvecl1 = dc.x*lmc.x + dc.y*lmc.y + dc.z*lmc.z;
			//csval1 = csval1/dl
			csval1 = Z_DIV(csval1, Z_MAKE(dl, (double)0.0e0));
		}
		// ------------ Zamiast prcodeury cwire1 -------------------


		// ------------ Zamiast prcodeury cwire2 -------------------
		ns = ips[k + kk + j + nsmax];
		nlr = ise[ns - 1 + nsmax];
		an = rad[ns - 1];
		double3 rnr = make_double3((double)x[nlr - 1], (double)y[nlr - 1], (double)z[nlr - 1]);
		rn = make_double3((double)x[np - 1], (double)y[np - 1], (double)z[np - 1]);
		rd = make_double3(rnr.x - rn.x, rnr.y - rn.y, rnr.z - rn.z);
		dl = sqrt(rd.x*rd.x + rd.y*rd.y + rd.z*rd.z);
		//dlh = (double)0.5e0*dl;
		dc = make_double3(rd.x / dl, rd.y / dl, rd.z / dl);

		ise_ms1 = ise[ms - 1];
		ise_ms2 = ise[ms - 1 + nsmax];
		ise_ns1 = ise[ns - 1];
		ise_ns2 = ise[ns - 1 + nsmax];

		if (war_gpu_pec_re(ise_ms1, ise_ms2, ise_ns1, ise_ns2)) {
			//call dcqg(0.d0,dlh,cknn,cv1,4)
			//nord=4		    
			// --------------- cv1 ------------------
			cvval2 = Z_MAKE((double)0.0e0, (double)0.0e0);
			csval2 = Z_MAKE((double)0.0e0, (double)0.0e0);
		}
		else {
			//call dcqg(0.d0,dl,ckmn,cv,2)
			//nord=2
			double dcqg_a = (double)0.5e0*dl;
			double dcqg_b = dl;
			double dcqg_c = (double)0.288675134594812882e0 * dcqg_b;
			//ifun=2
			//y=b*0.5d0*(f(a+c)+f(a-c))
			//dla s=a+c
			double ckmn_s = dcqg_a + dcqg_c;
			double3 ckmn_dr = make_double3(rmc.x - rn.x - ckmn_s*dc.x, rmc.y - rn.y - ckmn_s*dc.y, rmc.z - rn.z - ckmn_s*dc.z);
			double ckmn_r = sqrt(an*an + ckmn_dr.x*ckmn_dr.x + ckmn_dr.y*ckmn_dr.y + ckmn_dr.z*ckmn_dr.z);
			//stosujemy wz�r Eulera exp(-ja)=cos(a)-jsin(a) => a=ka*r 
			//c_jkr = cdexp(-cj*ka*ckmn_r)/ckmn_r
			cuDoubleComplex c_jkr = Z_MAKE(cos(ka*ckmn_r) / ckmn_r, -sin(ka*ckmn_r) / ckmn_r);
			//cvval2 = dcmplx((dl-ckmn_s)/dl,0.d0)*c_jkr	!ifun=2
			cvval2 = Z_MUL(Z_MAKE((dl - ckmn_s) / dl, (double)0.0e0), c_jkr);
			//csval1 = c_jkr	                            !ifun=3
			csval2 = c_jkr;
			//dla s=a-c
			ckmn_s = dcqg_a - dcqg_c;
			ckmn_dr = make_double3(rmc.x - rn.x - ckmn_s*dc.x, rmc.y - rn.y - ckmn_s*dc.y, rmc.z - rn.z - ckmn_s*dc.z);
			ckmn_r = sqrt(an*an + ckmn_dr.x*ckmn_dr.x + ckmn_dr.y*ckmn_dr.y + ckmn_dr.z*ckmn_dr.z);
			//c_jkr = exp(-cj*ka*ckmn_r)/ckmn_r
			//c_jkr = dcmplx( cos(ka*ckmn_r), -sin(ka*ckmn_r) )
			//c_jkr = c_jkr/ckmn_r
			c_jkr = Z_MAKE(cos(ka*ckmn_r) / ckmn_r, -sin(ka*ckmn_r) / ckmn_r);
			//cvval2 =  cvval2 + dcmplx((dl-ckmn_s)/dl,0.d0)*c_jkr	!ifun=2	
			//csval2 = csval2 + c_jkr								!ifun=3		
			cvval2 = Z_ADD(Z_MUL(Z_MAKE((dl - ckmn_s) / dl, (double)0.0e0), c_jkr), cvval2);
			csval2 = Z_ADD(csval2, c_jkr);

			cvval2 = Z_MUL(Z_MAKE(dcqg_b*(double)0.5e0, (double)0.0e0), cvval2);
			csval2 = Z_MUL(Z_MAKE(dcqg_b*(double)0.5e0, (double)0.0e0), csval2);

			//tvecl2 = dc%x*lmc%x+dc%y*lmc%y+dc%z*lmc%z      ! obliczane !lm+ * 1n- dla k=1,	obliczane  !lm- * 1n- dla k=2
			tvecl2 = dc.x*lmc.x + dc.y*lmc.y + dc.z*lmc.z;
			//csval2 = csval2/dl
			csval2 = Z_DIV(csval2, Z_MAKE(dl, (double)0.0e0));
		}
		// ------------ Zamiast prcodeury cwire2 -------------------

		//cvval = dcmplx(tvecl1,0.d0)*cvval1+dcmplx(tvecl2,0.d0)*cvval2   !A(rm+)	dla k=1, A(rm-) dla k=2		  
		//csval = csval1 - csval2					!O(rm+) dla k=1, !O(rm-) dla k=2
		cuDoubleComplex cvval = Z_MUL(Z_MAKE(tvecl1, (double)0.0e0), cvval1);
		cvval = Z_ADD(Z_MUL(Z_MAKE(tvecl2, (double)0.0e0), cvval2), cvval);
		cuDoubleComplex csval = Z_SUB(csval1, csval2);


		//cz_ij = cz_ij+ka*ka*cvval + znak * csval
		cuDoubleComplex cz_ij = Z_MUL(Z_MAKE((double)znak, (double)0.0e0), csval);
		cz_ij = Z_ADD(Z_MUL(Z_MAKE(ka*ka, (double)0.0e0), cvval), cz_ij);

		//cz_ij = cz_ij*cj*eta0/(4.0e0*pi*ka)
		cz_ij = Z_MUL(cz_ij, Z_MAKE(eta0 / ((double)4.0e0*pi*ka), (double)0.0e0));
		cz_ij = Z_MUL(cz_ij, Z_ONE_IMG);



		//cz_d[i_m+j_n*n] = Z_MAKE(10.0f,0.0f);
		//cz_d[i_m+j_n*n] = Z_MUL ( cz_d[i_m+j_n*n] , (Z_MAKE(10.0f,0.0f));
		//cz_d[i_m+j_n*n] = Z_MAKE((double)i_m,(double)j_n);
		//cz_d[i_m+j_n*n]=make_cuDoubleComplex( (i_m),(j_n) );
		//cz_d[ i+j*n ] =  Z_MAKE(float(i*n+j),0.0f) ;
		//cz[ i+j*n ] = Z_MAKE(x[mp-1],z[mp-1]);
		//cz[ i+j*n ] = Z_MAKE(dc.x,dc.z);
		//cz[ i+j*n ] = cknn_cbnd1;
		//cz[ i+j*n ] = Z_MAKE(cknn_res,cknn_bnd2);
		//cz[ i+j*n ] = cz_ij;
		cz[i + (j + kk)*n] = Z_ADD(cz[i + (j + kk)*n], cz_ij);
		//cz[ i+j*n ] =  Z_ADD(c1.z,Z_MAKE(float(i*n+j),0.0f)) ;
	}
}

//----------------------------------------------------------------------

__global__ void wim_gpu_ckmn_re2(cuDoubleComplex *cz, float *x, float *y, float *z,
	float *rad, int *ipc, int *ips, int *ise, int nsmax,
	int n, double ka, double eta0,
	int realSubSliceSize, int k, int kk)
{

	double pi = acos(double(-1.e0));
	cuDoubleComplex cvval1, csval1, cvval2, csval2;
	double tvecl1, tvecl2;

	// Block index
	int bx = blockIdx.x;
	int by = blockIdx.y;

	// Thread index
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	// Accumulate row i of A and column j of B
	int i = tx + bx * blockDim.x;
	int j = ty + by * blockDim.y;

	//c16vec c1 = make_c16vec(Z_ONE,Z_MAKE(0.0f, 1.0f),Z_MUL(Z_MAKE(2.0f, 0.0f),Z_ONE));

	if (i < n  && j < realSubSliceSize) {
		int mp = abs(ipc[i]);
		int ms = ips[i + nsmax];
		int mlr = ise[ms - 1 + nsmax];
		//double am = rad[ms-1];
		double3 rmp = make_double3((double)x[mp - 1], (double)y[mp - 1], (double)z[mp - 1]);
		double3 rmlr = make_double3((double)x[mlr - 1], (double)y[mlr - 1], (double)z[mlr - 1]);
		double3 rmc = make_double3((double)0.5e0*(rmlr.x + rmp.x), (double)0.5e0*(rmlr.y + rmp.y), (double)0.5e0*(rmlr.z + rmp.z));
		float znak = 1.0f;
		double3 lmc = make_double3((double)znak*(rmc.x - rmp.x), (double)znak*(rmc.y - rmp.y), (double)znak*(rmc.z - rmp.z));

		// ------------ Zamiast prcodeury cwire1 -------------------
		//ifun=1 i ifun=2=> res_cv
		//fun=3 => res_cs		  
		//call cwire1_gpu(ms,n,cvval1)             ! obliczanie calki na Sn+  pot.wek.
		//call cwire1(ms,n,csval1)

		//int np = abs(ipc[k-1+kk-1+j]);
		//int ns = ips[k-1+kk-1+j];
		int np = abs(ipc[k + kk + j]);
		int ns = ips[k + kk + j];
		double3 rnp = make_double3((double)x[np - 1], (double)y[np - 1], (double)z[np - 1]);
		int nlr = ise[ns - 1];
		double an = rad[ns - 1];
		double3 rn = make_double3((double)x[nlr - 1], (double)y[nlr - 1], (double)z[nlr - 1]);
		double3 rd = make_double3(rnp.x - rn.x, rnp.y - rn.y, rnp.z - rn.z);
		double dl = sqrt(rd.x*rd.x + rd.y*rd.y + rd.z*rd.z);
		//double dlh = (double)0.5e0*dl;
		double3 dc = make_double3(rd.x / dl, rd.y / dl, rd.z / dl);

		int ise_ms1 = ise[ms - 1];
		int ise_ms2 = ise[ms - 1 + nsmax];
		int ise_ns1 = ise[ns - 1];
		int ise_ns2 = ise[ns - 1 + nsmax];


		if (war_gpu_pec_re(ise_ms1, ise_ms2, ise_ns1, ise_ns2)) {
			//call dcqg(0.d0,dlh,cknn,cv1,4)
			//nord=4		    
			// --------------- cv1 ------------------
			cvval1 = Z_MAKE((double)0.0e0, (double)0.0e0);
			csval1 = Z_MAKE((double)0.0e0, (double)0.0e0);
		}
		else {
			//call dcqg(0.d0,dl,ckmn,cv,2)
			//nord=2
			double dcqg_a = (double)0.5e0*dl;
			double dcqg_b = dl;
			double dcqg_c = (double)0.288675134594812882e0 * dcqg_b;
			//ifun=1
			//y=b*0.5d0*(f(a+c)+f(a-c))
			//dla s=a+c
			double ckmn_s = dcqg_a + dcqg_c;
			double3 ckmn_dr = make_double3(rmc.x - rn.x - ckmn_s*dc.x, rmc.y - rn.y - ckmn_s*dc.y, rmc.z - rn.z - ckmn_s*dc.z);
			double ckmn_r = sqrt(an*an + ckmn_dr.x*ckmn_dr.x + ckmn_dr.y*ckmn_dr.y + ckmn_dr.z*ckmn_dr.z);
			//stosujemy wz�r Eulera exp(-ja)=cos(a)-jsin(a) => a=ka*r 
			//c_jkr = cdexp(-cj*ka*ckmn_r)/ckmn_r
			cuDoubleComplex c_jkr = Z_MAKE(cos(ka*ckmn_r) / ckmn_r, -sin(ka*ckmn_r) / ckmn_r);
			//cvval1 = dcmplx(ckmn_s/dl,0.d0)*c_jkr	!ifun=1
			cvval1 = Z_MUL(Z_MAKE(ckmn_s / dl, (double)0.0e0), c_jkr);
			//csval1 = c_jkr	                        !ifun=3
			csval1 = c_jkr;
			//dla s=a-c
			ckmn_s = dcqg_a - dcqg_c;
			ckmn_dr = make_double3(rmc.x - rn.x - ckmn_s*dc.x, rmc.y - rn.y - ckmn_s*dc.y, rmc.z - rn.z - ckmn_s*dc.z);
			ckmn_r = sqrt(an*an + ckmn_dr.x*ckmn_dr.x + ckmn_dr.y*ckmn_dr.y + ckmn_dr.z*ckmn_dr.z);
			//c_jkr = exp(-cj*ka*ckmn_r)/ckmn_r
			//c_jkr = dcmplx( cos(ka*ckmn_r), -sin(ka*ckmn_r) )
			//c_jkr = c_jkr/ckmn_r
			c_jkr = Z_MAKE(cos(ka*ckmn_r) / ckmn_r, -sin(ka*ckmn_r) / ckmn_r);
			//cvval1 =  cvval1 + dcmplx(ckmn_s/dl,0.d0)*c_jkr	!ifun=1	
			//csval1 = csval1 + c_jkr							!ifun=3	
			cvval1 = Z_ADD(Z_MUL(Z_MAKE(ckmn_s / dl, (double)0.0e0), c_jkr), cvval1);
			csval1 = Z_ADD(csval1, c_jkr);

			cvval1 = Z_MUL(Z_MAKE(dcqg_b*(double)0.5e0, (double)0.0e0), cvval1);
			csval1 = Z_MUL(Z_MAKE(dcqg_b*(double)0.5e0, (double)0.0e0), csval1);

			//tvecl1= dc%x*lmc%x+dc%y*lmc%y+dc%z*lmc%z       ! obliczane !lm+ * 1n+ dla k=1, obliczane !lm- * 1n+ dla k=1
			tvecl1 = dc.x*lmc.x + dc.y*lmc.y + dc.z*lmc.z;
			//csval1 = csval1/dl
			csval1 = Z_DIV(csval1, Z_MAKE(dl, (double)0.0e0));
		}
		// ------------ Zamiast prcodeury cwire1 -------------------


		// ------------ Zamiast prcodeury cwire2 -------------------
		ns = ips[k + kk + j + nsmax];
		nlr = ise[ns - 1 + nsmax];
		an = rad[ns - 1];
		double3 rnr = make_double3((double)x[nlr - 1], (double)y[nlr - 1], (double)z[nlr - 1]);
		rn = make_double3((double)x[np - 1], (double)y[np - 1], (double)z[np - 1]);
		rd = make_double3(rnr.x - rn.x, rnr.y - rn.y, rnr.z - rn.z);
		dl = sqrt(rd.x*rd.x + rd.y*rd.y + rd.z*rd.z);
		//dlh = (double)0.5e0*dl;
		dc = make_double3(rd.x / dl, rd.y / dl, rd.z / dl);

		ise_ms1 = ise[ms - 1];
		ise_ms2 = ise[ms - 1 + nsmax];
		ise_ns1 = ise[ns - 1];
		ise_ns2 = ise[ns - 1 + nsmax];

		if (war_gpu_pec_re(ise_ms1, ise_ms2, ise_ns1, ise_ns2)) {
			//call dcqg(0.d0,dlh,cknn,cv1,4)
			//nord=4		    
			// --------------- cv1 ------------------
			cvval2 = Z_MAKE((double)0.0e0, (double)0.0e0);
			csval2 = Z_MAKE((double)0.0e0, (double)0.0e0);
		}
		else {
			//call dcqg(0.d0,dl,ckmn,cv,2)
			//nord=2
			double dcqg_a = (double)0.5e0*dl;
			double dcqg_b = dl;
			double dcqg_c = (double)0.288675134594812882e0 * dcqg_b;
			//ifun=2
			//y=b*0.5d0*(f(a+c)+f(a-c))
			//dla s=a+c
			double ckmn_s = dcqg_a + dcqg_c;
			double3 ckmn_dr = make_double3(rmc.x - rn.x - ckmn_s*dc.x, rmc.y - rn.y - ckmn_s*dc.y, rmc.z - rn.z - ckmn_s*dc.z);
			double ckmn_r = sqrt(an*an + ckmn_dr.x*ckmn_dr.x + ckmn_dr.y*ckmn_dr.y + ckmn_dr.z*ckmn_dr.z);
			//stosujemy wz�r Eulera exp(-ja)=cos(a)-jsin(a) => a=ka*r 
			//c_jkr = cdexp(-cj*ka*ckmn_r)/ckmn_r
			cuDoubleComplex c_jkr = Z_MAKE(cos(ka*ckmn_r) / ckmn_r, -sin(ka*ckmn_r) / ckmn_r);
			//cvval2 = dcmplx((dl-ckmn_s)/dl,0.d0)*c_jkr	!ifun=2
			cvval2 = Z_MUL(Z_MAKE((dl - ckmn_s) / dl, (double)0.0e0), c_jkr);
			//csval1 = c_jkr	                            !ifun=3
			csval2 = c_jkr;
			//dla s=a-c
			ckmn_s = dcqg_a - dcqg_c;
			ckmn_dr = make_double3(rmc.x - rn.x - ckmn_s*dc.x, rmc.y - rn.y - ckmn_s*dc.y, rmc.z - rn.z - ckmn_s*dc.z);
			ckmn_r = sqrt(an*an + ckmn_dr.x*ckmn_dr.x + ckmn_dr.y*ckmn_dr.y + ckmn_dr.z*ckmn_dr.z);
			//c_jkr = exp(-cj*ka*ckmn_r)/ckmn_r
			//c_jkr = dcmplx( cos(ka*ckmn_r), -sin(ka*ckmn_r) )
			//c_jkr = c_jkr/ckmn_r
			c_jkr = Z_MAKE(cos(ka*ckmn_r) / ckmn_r, -sin(ka*ckmn_r) / ckmn_r);
			//cvval2 =  cvval2 + dcmplx((dl-ckmn_s)/dl,0.d0)*c_jkr	!ifun=2	
			//csval2 = csval2 + c_jkr								!ifun=3		
			cvval2 = Z_ADD(Z_MUL(Z_MAKE((dl - ckmn_s) / dl, (double)0.0e0), c_jkr), cvval2);
			csval2 = Z_ADD(csval2, c_jkr);

			cvval2 = Z_MUL(Z_MAKE(dcqg_b*(double)0.5e0, (double)0.0e0), cvval2);
			csval2 = Z_MUL(Z_MAKE(dcqg_b*(double)0.5e0, (double)0.0e0), csval2);

			//tvecl2 = dc%x*lmc%x+dc%y*lmc%y+dc%z*lmc%z      ! obliczane !lm+ * 1n- dla k=1,	obliczane  !lm- * 1n- dla k=2
			tvecl2 = dc.x*lmc.x + dc.y*lmc.y + dc.z*lmc.z;
			//csval2 = csval2/dl
			csval2 = Z_DIV(csval2, Z_MAKE(dl, (double)0.0e0));
		}
		// ------------ Zamiast prcodeury cwire2 -------------------

		//cvval = dcmplx(tvecl1,0.d0)*cvval1+dcmplx(tvecl2,0.d0)*cvval2   !A(rm+)	dla k=1, A(rm-) dla k=2		  
		//csval = csval1 - csval2					!O(rm+) dla k=1, !O(rm-) dla k=2
		cuDoubleComplex cvval = Z_MUL(Z_MAKE(tvecl1, (double)0.0e0), cvval1);
		cvval = Z_ADD(Z_MUL(Z_MAKE(tvecl2, (double)0.0e0), cvval2), cvval);
		cuDoubleComplex csval = Z_SUB(csval1, csval2);


		//cz_ij = cz_ij+ka*ka*cvval + znak * csval
		cuDoubleComplex cz_ij = Z_MUL(Z_MAKE((double)znak, (double)0.0e0), csval);
		cz_ij = Z_ADD(Z_MUL(Z_MAKE(ka*ka, (double)0.0e0), cvval), cz_ij);

		//cz_ij = cz_ij*cj*eta0/(4.0e0*pi*ka)
		cz_ij = Z_MUL(cz_ij, Z_MAKE(eta0 / ((double)4.0e0*pi*ka), (double)0.0e0));
		cz_ij = Z_MUL(cz_ij, Z_ONE_IMG);



		//cz_d[i_m+j_n*n] = Z_MAKE(10.0f,0.0f);
		//cz_d[i_m+j_n*n] = Z_MUL ( cz_d[i_m+j_n*n] , (Z_MAKE(10.0f,0.0f));
		//cz_d[i_m+j_n*n] = Z_MAKE((double)i_m,(double)j_n);
		//cz_d[i_m+j_n*n]=make_cuDoubleComplex( (i_m),(j_n) );
		//cz_d[ i+j*n ] =  Z_MAKE(float(i*n+j),0.0f) ;
		//cz[ i+j*n ] = Z_MAKE(x[mp-1],z[mp-1]);
		//cz[ i+j*n ] = Z_MAKE(dc.x,dc.z);
		//cz[ i+j*n ] = cknn_cbnd1;
		//cz[ i+j*n ] = Z_MAKE(cknn_res,cknn_bnd2);
		//cz[ i+j*n ] = cz_ij;
		cz[i + (j + kk)*n] = Z_ADD(cz[i + (j + kk)*n], cz_ij);
		//cz[ i+j*n ] =  Z_ADD(c1.z,Z_MAKE(float(i*n+j),0.0f)) ;
	}
}

//----------------------------------------------------------------------

__global__ void wim_gpu_ckmn_im1(cuDoubleComplex *cz, float *x, float *y, float *z,
	float *rad, int *ipc, int *ips, int *ise, int nsmax,
	int n, double ka, double eta0,
	float sgnx, float sgny, float sgnz, float sgnenv,
	int realSubSliceSize, int k, int kk)
{

	double pi = acos(double(-1.e0));
	cuDoubleComplex cvval1, csval1, cvval2, csval2;
	double tvecl1, tvecl2;

	// Block index
	int bx = blockIdx.x;
	int by = blockIdx.y;

	// Thread index
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	// Accumulate row i of A and column j of B
	int i = tx + bx * blockDim.x;
	int j = ty + by * blockDim.y;

	//c16vec c1 = make_c16vec(Z_ONE,Z_MAKE(0.0f, 1.0f),Z_MUL(Z_MAKE(2.0f, 0.0f),Z_ONE));

	if (i < n  && j < realSubSliceSize) {
		int mp = abs(ipc[i]);
		int ms = ips[i];
		int mlr = ise[ms - 1];
		//double am = rad[ms-1];
		double3 rmp = make_double3((double)x[mp - 1], (double)y[mp - 1], (double)z[mp - 1]);
		double3 rmlr = make_double3((double)x[mlr - 1], (double)y[mlr - 1], (double)z[mlr - 1]);
		double3 rmc = make_double3((double)0.5e0*(rmlr.x + rmp.x), (double)0.5e0*(rmlr.y + rmp.y), (double)0.5e0*(rmlr.z + rmp.z));
		float znak = -1.0f;
		double3 lmc = make_double3((double)znak*(rmc.x - rmp.x), (double)znak*(rmc.y - rmp.y), (double)znak*(rmc.z - rmp.z));

		// ------------ Zamiast prcodeury cwire1 -------------------
		//ifun=1 i ifun=2=> res_cv
		//fun=3 => res_cs		  
		//call cwire1_gpu(ms,n,cvval1)             ! obliczanie calki na Sn+  pot.wek.
		//call cwire1(ms,n,csval1)

		//int np = abs(ipc[k-1+kk-1+j]);
		//int ns = ips[k-1+kk-1+j];
		int np = abs(ipc[k + kk + j]);
		int ns = ips[k + kk + j];
		double3 rnp = make_double3((double)sgnx*(double)x[np - 1], (double)sgny*(double)y[np - 1], (double)sgnz*(double)z[np - 1]);
		int nlr = ise[ns - 1];
		double an = rad[ns - 1];
		double3 rn = make_double3((double)sgnx*(double)x[nlr - 1], (double)sgny*(double)y[nlr - 1], (double)sgnz*(double)z[nlr - 1]);
		double3 rd = make_double3(rnp.x - rn.x, rnp.y - rn.y, rnp.z - rn.z);
		double dl = sqrt(rd.x*rd.x + rd.y*rd.y + rd.z*rd.z);
		//double dlh = (double)0.5e0*dl;
		double3 dc = make_double3(rd.x / dl, rd.y / dl, rd.z / dl);

		int ise_ms1 = ise[ms - 1];
		int ise_ms2 = ise[ms - 1 + nsmax];
		int ise_ns1 = ise[ns - 1];
		int ise_ns2 = ise[ns - 1 + nsmax];
		//float z_ise_ms1 = z[ise[ms-1]];
		//float z_ise_ms2 = z[ise[ms-1+nsmax]];
		float z_ise_ms1 = z[ise[ms - 1] - 1];
		float z_ise_ms2 = z[ise[ms - 1 + nsmax] - 1];

		if (war_gpu_pec_im(ise_ms1, ise_ms2, ise_ns1, ise_ns2, z_ise_ms1, z_ise_ms2)) {
			//call dcqg(0.d0,dlh,cknn,cv1,4)
			//nord=4		    
			// --------------- cv1 ------------------
			cvval1 = Z_MAKE((double)0.0e0, (double)0.0e0);
			csval1 = Z_MAKE((double)0.0e0, (double)0.0e0);
		}
		else {
			//call dcqg(0.d0,dl,ckmn,cv,2)
			//nord=2
			double dcqg_a = (double)0.5e0*dl;
			double dcqg_b = dl;
			double dcqg_c = (double)0.288675134594812882e0 * dcqg_b;
			//ifun=1
			//y=b*0.5d0*(f(a+c)+f(a-c))
			//dla s=a+c
			double ckmn_s = dcqg_a + dcqg_c;
			double3 ckmn_dr = make_double3(rmc.x - rn.x - ckmn_s*dc.x, rmc.y - rn.y - ckmn_s*dc.y, rmc.z - rn.z - ckmn_s*dc.z);
			double ckmn_r = sqrt(an*an + ckmn_dr.x*ckmn_dr.x + ckmn_dr.y*ckmn_dr.y + ckmn_dr.z*ckmn_dr.z);
			//stosujemy wz�r Eulera exp(-ja)=cos(a)-jsin(a) => a=ka*r 
			//c_jkr = cdexp(-cj*ka*ckmn_r)/ckmn_r
			cuDoubleComplex c_jkr = Z_MAKE(cos(ka*ckmn_r) / ckmn_r, -sin(ka*ckmn_r) / ckmn_r);
			//cvval1 = dcmplx(ckmn_s/dl,0.d0)*c_jkr	!ifun=1
			cvval1 = Z_MUL(Z_MAKE(ckmn_s / dl, (double)0.0e0), c_jkr);
			//csval1 = c_jkr	                        !ifun=3
			csval1 = c_jkr;
			//dla s=a-c
			ckmn_s = dcqg_a - dcqg_c;
			ckmn_dr = make_double3(rmc.x - rn.x - ckmn_s*dc.x, rmc.y - rn.y - ckmn_s*dc.y, rmc.z - rn.z - ckmn_s*dc.z);
			ckmn_r = sqrt(an*an + ckmn_dr.x*ckmn_dr.x + ckmn_dr.y*ckmn_dr.y + ckmn_dr.z*ckmn_dr.z);
			//c_jkr = exp(-cj*ka*ckmn_r)/ckmn_r
			//c_jkr = dcmplx( cos(ka*ckmn_r), -sin(ka*ckmn_r) )
			//c_jkr = c_jkr/ckmn_r
			c_jkr = Z_MAKE(cos(ka*ckmn_r) / ckmn_r, -sin(ka*ckmn_r) / ckmn_r);
			//cvval1 =  cvval1 + dcmplx(ckmn_s/dl,0.d0)*c_jkr	!ifun=1	
			//csval1 = csval1 + c_jkr							!ifun=3	
			cvval1 = Z_ADD(Z_MUL(Z_MAKE(ckmn_s / dl, (double)0.0e0), c_jkr), cvval1);
			csval1 = Z_ADD(csval1, c_jkr);

			cvval1 = Z_MUL(Z_MAKE(dcqg_b*(double)0.5e0, (double)0.0e0), cvval1);
			csval1 = Z_MUL(Z_MAKE(dcqg_b*(double)0.5e0, (double)0.0e0), csval1);

			//tvecl1= dc%x*lmc%x+dc%y*lmc%y+dc%z*lmc%z       ! obliczane !lm+ * 1n+ dla k=1, obliczane !lm- * 1n+ dla k=1
			tvecl1 = dc.x*lmc.x + dc.y*lmc.y + dc.z*lmc.z;
			//csval1 = csval1/dl
			csval1 = Z_DIV(csval1, Z_MAKE(dl, (double)0.0e0));
		}
		// ------------ Zamiast prcodeury cwire1 -------------------


		// ------------ Zamiast prcodeury cwire2 -------------------
		ns = ips[k + kk + j + nsmax];
		nlr = ise[ns - 1 + nsmax];
		an = rad[ns - 1];
		double3 rnr = make_double3((double)sgnx*(double)x[nlr - 1], (double)sgny*(double)y[nlr - 1], (double)sgnz*(double)z[nlr - 1]);
		rn = make_double3((double)sgnx*(double)x[np - 1], (double)sgny*(double)y[np - 1], (double)sgnz*(double)z[np - 1]);
		rd = make_double3(rnr.x - rn.x, rnr.y - rn.y, rnr.z - rn.z);
		dl = sqrt(rd.x*rd.x + rd.y*rd.y + rd.z*rd.z);
		//dlh = (double)0.5e0*dl;
		dc = make_double3(rd.x / dl, rd.y / dl, rd.z / dl);

		ise_ms1 = ise[ms - 1];
		ise_ms2 = ise[ms - 1 + nsmax];
		ise_ns1 = ise[ns - 1];
		ise_ns2 = ise[ns - 1 + nsmax];
		//z_ise_ms1 = z[ise[ms-1]];
		//z_ise_ms2 = z[ise[ms-1+nsmax]];
		z_ise_ms1 = z[ise[ms - 1] - 1];
		z_ise_ms2 = z[ise[ms - 1 + nsmax] - 1];

		if (war_gpu_pec_im(ise_ms1, ise_ms2, ise_ns1, ise_ns2, z_ise_ms1, z_ise_ms2)) {
			//call dcqg(0.d0,dlh,cknn,cv1,4)
			//nord=4		    
			// --------------- cv1 ------------------
			cvval2 = Z_MAKE((double)0.0e0, (double)0.0e0);
			csval2 = Z_MAKE((double)0.0e0, (double)0.0e0);
		}
		else {
			//call dcqg(0.d0,dl,ckmn,cv,2)
			//nord=2
			double dcqg_a = (double)0.5e0*dl;
			double dcqg_b = dl;
			double dcqg_c = (double)0.288675134594812882e0 * dcqg_b;
			//ifun=2
			//y=b*0.5d0*(f(a+c)+f(a-c))
			//dla s=a+c
			double ckmn_s = dcqg_a + dcqg_c;
			double3 ckmn_dr = make_double3(rmc.x - rn.x - ckmn_s*dc.x, rmc.y - rn.y - ckmn_s*dc.y, rmc.z - rn.z - ckmn_s*dc.z);
			double ckmn_r = sqrt(an*an + ckmn_dr.x*ckmn_dr.x + ckmn_dr.y*ckmn_dr.y + ckmn_dr.z*ckmn_dr.z);
			//stosujemy wz�r Eulera exp(-ja)=cos(a)-jsin(a) => a=ka*r 
			//c_jkr = cdexp(-cj*ka*ckmn_r)/ckmn_r
			cuDoubleComplex c_jkr = Z_MAKE(cos(ka*ckmn_r) / ckmn_r, -sin(ka*ckmn_r) / ckmn_r);
			//cvval2 = dcmplx((dl-ckmn_s)/dl,0.d0)*c_jkr	!ifun=2
			cvval2 = Z_MUL(Z_MAKE((dl - ckmn_s) / dl, (double)0.0e0), c_jkr);
			//csval1 = c_jkr	                            !ifun=3
			csval2 = c_jkr;
			//dla s=a-c
			ckmn_s = dcqg_a - dcqg_c;
			ckmn_dr = make_double3(rmc.x - rn.x - ckmn_s*dc.x, rmc.y - rn.y - ckmn_s*dc.y, rmc.z - rn.z - ckmn_s*dc.z);
			ckmn_r = sqrt(an*an + ckmn_dr.x*ckmn_dr.x + ckmn_dr.y*ckmn_dr.y + ckmn_dr.z*ckmn_dr.z);
			//c_jkr = exp(-cj*ka*ckmn_r)/ckmn_r
			//c_jkr = dcmplx( cos(ka*ckmn_r), -sin(ka*ckmn_r) )
			//c_jkr = c_jkr/ckmn_r
			c_jkr = Z_MAKE(cos(ka*ckmn_r) / ckmn_r, -sin(ka*ckmn_r) / ckmn_r);
			//cvval2 =  cvval2 + dcmplx((dl-ckmn_s)/dl,0.d0)*c_jkr	!ifun=2	
			//csval2 = csval2 + c_jkr								!ifun=3		
			cvval2 = Z_ADD(Z_MUL(Z_MAKE((dl - ckmn_s) / dl, (double)0.0e0), c_jkr), cvval2);
			csval2 = Z_ADD(csval2, c_jkr);

			cvval2 = Z_MUL(Z_MAKE(dcqg_b*(double)0.5e0, (double)0.0e0), cvval2);
			csval2 = Z_MUL(Z_MAKE(dcqg_b*(double)0.5e0, (double)0.0e0), csval2);

			//tvecl2 = dc%x*lmc%x+dc%y*lmc%y+dc%z*lmc%z      ! obliczane !lm+ * 1n- dla k=1,	obliczane  !lm- * 1n- dla k=2
			tvecl2 = dc.x*lmc.x + dc.y*lmc.y + dc.z*lmc.z;
			//csval2 = csval2/dl
			csval2 = Z_DIV(csval2, Z_MAKE(dl, (double)0.0e0));
		}
		// ------------ Zamiast prcodeury cwire2 -------------------

		//cvval = sgnenv*dcmplx(tvecl1,0.d0)*cvval1 + sgnenv*dcmplx(tvecl2,0.d0)*cvval2   !A(rm+)	dla k=1, A(rm-) dla k=2		  
		//csval = csval1*sgnenv - sgnenv*csval2                      					!O(rm+) dla k=1, !O(rm-) dla k=2
		//
		//cuDoubleComplex cvval = Z_MUL( Z_MAKE((double)sgnenv*tvecl1,(double)0.0e0),cvval1 );
		//cvval = Z_ADD( Z_MUL(Z_MAKE((double)sgnenv*tvecl2,(double)0.0e0),cvval2),cvval );
		cuDoubleComplex cvval = Z_MUL(Z_MAKE(tvecl1, (double)0.0e0), cvval1);
		cvval = Z_ADD(Z_MUL(Z_MAKE(tvecl2, (double)0.0e0), cvval2), cvval);
		cvval = Z_MUL(Z_MAKE((double)sgnenv, (double)0.0e0), cvval);
		//cuDoubleComplex csval = Z_MUL( Z_MAKE((double)sgnenv,(double)0.0e0),csval2 );
		//csval = Z_SUB( Z_MUL( Z_MAKE((double)sgnenv,(double)0.0e0),csval1 ),csval );
		cuDoubleComplex csval = Z_MUL(Z_MAKE((double)sgnenv, (double)0.0e0), Z_SUB(csval1, csval2));

		//cz_ij = cz_ij+ka*ka*cvval + znak * csval
		cuDoubleComplex cz_ij = Z_MUL(Z_MAKE((double)znak, (double)0.0e0), csval);
		cz_ij = Z_ADD(Z_MUL(Z_MAKE(ka*ka, (double)0.0e0), cvval), cz_ij);

		//cz_ij = cz_ij*cj*eta0/(4.0e0*pi*ka)
		cz_ij = Z_MUL(cz_ij, Z_MAKE(eta0 / ((double)4.0e0*pi*ka), (double)0.0e0));
		cz_ij = Z_MUL(cz_ij, Z_ONE_IMG);



		//cz_d[i_m+j_n*n] = Z_MAKE(10.0f,0.0f);
		//cz_d[i_m+j_n*n] = Z_MUL ( cz_d[i_m+j_n*n] , (Z_MAKE(10.0f,0.0f));
		//cz_d[i_m+j_n*n] = Z_MAKE((double)i_m,(double)j_n);
		//cz_d[i_m+j_n*n]=make_cuDoubleComplex( (i_m),(j_n) );
		//cz_d[ i+j*n ] =  Z_MAKE(float(i*n+j),0.0f) ;
		//cz[ i+j*n ] = Z_MAKE(x[mp-1],z[mp-1]);
		//cz[ i+j*n ] = Z_MAKE(dc.x,dc.z);
		//cz[ i+j*n ] = cknn_cbnd1;
		//cz[ i+j*n ] = Z_MAKE((double)sgnenv,tvecl2);
		//cz[ i+j*n ] = cvval;
		//cz[ i+j*n ] = cz_ij;
		cz[i + (j + kk)*n] = Z_ADD(cz[i + (j + kk)*n], cz_ij);
		//cz[ i+j*n ] =  Z_ADD(c1.z,Z_MAKE(float(i*n+j),0.0f)) ;
	}
}

//----------------------------------------------------------------------

__global__ void wim_gpu_ckmn_im2(cuDoubleComplex *cz, float *x, float *y, float *z,
	float *rad, int *ipc, int *ips, int *ise, int nsmax,
	int n, double ka, double eta0,
	float sgnx, float sgny, float sgnz, float sgnenv,
	int realSubSliceSize, int k, int kk)
{

	double pi = acos(double(-1.e0));
	cuDoubleComplex cvval1, csval1, cvval2, csval2;
	double tvecl1, tvecl2;

	// Block index
	int bx = blockIdx.x;
	int by = blockIdx.y;

	// Thread index
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	// Accumulate row i of A and column j of B
	int i = tx + bx * blockDim.x;
	int j = ty + by * blockDim.y;

	//c16vec c1 = make_c16vec(Z_ONE,Z_MAKE(0.0f, 1.0f),Z_MUL(Z_MAKE(2.0f, 0.0f),Z_ONE));

	if (i < n  && j < realSubSliceSize) {
		int mp = abs(ipc[i]);
		int ms = ips[i + nsmax];
		int mlr = ise[ms - 1 + nsmax];
		//double am = rad[ms-1];
		double3 rmp = make_double3((double)x[mp - 1], (double)y[mp - 1], (double)z[mp - 1]);
		double3 rmlr = make_double3((double)x[mlr - 1], (double)y[mlr - 1], (double)z[mlr - 1]);
		double3 rmc = make_double3((double)0.5e0*(rmlr.x + rmp.x), (double)0.5e0*(rmlr.y + rmp.y), (double)0.5e0*(rmlr.z + rmp.z));
		float znak = 1.0f;
		double3 lmc = make_double3((double)znak*(rmc.x - rmp.x), (double)znak*(rmc.y - rmp.y), (double)znak*(rmc.z - rmp.z));

		// ------------ Zamiast prcodeury cwire1 -------------------
		//ifun=1 i ifun=2=> res_cv
		//fun=3 => res_cs		  
		//call cwire1_gpu(ms,n,cvval1)             ! obliczanie calki na Sn+  pot.wek.
		//call cwire1(ms,n,csval1)

		//int np = abs(ipc[k-1+kk-1+j]);
		//int ns = ips[k-1+kk-1+j];
		int np = abs(ipc[k + kk + j]);
		int ns = ips[k + kk + j];
		double3 rnp = make_double3((double)sgnx*(double)x[np - 1], (double)sgny*(double)y[np - 1], (double)sgnz*(double)z[np - 1]);
		int nlr = ise[ns - 1];
		double an = rad[ns - 1];
		double3 rn = make_double3((double)sgnx*(double)x[nlr - 1], (double)sgny*(double)y[nlr - 1], (double)sgnz*(double)z[nlr - 1]);
		double3 rd = make_double3(rnp.x - rn.x, rnp.y - rn.y, rnp.z - rn.z);
		double dl = sqrt(rd.x*rd.x + rd.y*rd.y + rd.z*rd.z);
		//double dlh = (double)0.5e0*dl;
		double3 dc = make_double3(rd.x / dl, rd.y / dl, rd.z / dl);

		int ise_ms1 = ise[ms - 1];
		int ise_ms2 = ise[ms - 1 + nsmax];
		int ise_ns1 = ise[ns - 1];
		int ise_ns2 = ise[ns - 1 + nsmax];
		//float z_ise_ms1 = z[ise[ms-1]];
		//float z_ise_ms2 = z[ise[ms-1+nsmax]];
		float z_ise_ms1 = z[ise[ms - 1] - 1];
		float z_ise_ms2 = z[ise[ms - 1 + nsmax] - 1];

		if (war_gpu_pec_im(ise_ms1, ise_ms2, ise_ns1, ise_ns2, z_ise_ms1, z_ise_ms2)) {
			//call dcqg(0.d0,dlh,cknn,cv1,4)
			//nord=4		    
			// --------------- cv1 ------------------
			cvval1 = Z_MAKE((double)0.0e0, (double)0.0e0);
			csval1 = Z_MAKE((double)0.0e0, (double)0.0e0);
		}
		else {
			//call dcqg(0.d0,dl,ckmn,cv,2)
			//nord=2
			double dcqg_a = (double)0.5e0*dl;
			double dcqg_b = dl;
			double dcqg_c = (double)0.288675134594812882e0 * dcqg_b;
			//ifun=1
			//y=b*0.5d0*(f(a+c)+f(a-c))
			//dla s=a+c
			double ckmn_s = dcqg_a + dcqg_c;
			double3 ckmn_dr = make_double3(rmc.x - rn.x - ckmn_s*dc.x, rmc.y - rn.y - ckmn_s*dc.y, rmc.z - rn.z - ckmn_s*dc.z);
			double ckmn_r = sqrt(an*an + ckmn_dr.x*ckmn_dr.x + ckmn_dr.y*ckmn_dr.y + ckmn_dr.z*ckmn_dr.z);
			//stosujemy wz�r Eulera exp(-ja)=cos(a)-jsin(a) => a=ka*r 
			//c_jkr = cdexp(-cj*ka*ckmn_r)/ckmn_r
			cuDoubleComplex c_jkr = Z_MAKE(cos(ka*ckmn_r) / ckmn_r, -sin(ka*ckmn_r) / ckmn_r);
			//cvval1 = dcmplx(ckmn_s/dl,0.d0)*c_jkr	!ifun=1
			cvval1 = Z_MUL(Z_MAKE(ckmn_s / dl, (double)0.0e0), c_jkr);
			//csval1 = c_jkr	                        !ifun=3
			csval1 = c_jkr;
			//dla s=a-c
			ckmn_s = dcqg_a - dcqg_c;
			ckmn_dr = make_double3(rmc.x - rn.x - ckmn_s*dc.x, rmc.y - rn.y - ckmn_s*dc.y, rmc.z - rn.z - ckmn_s*dc.z);
			ckmn_r = sqrt(an*an + ckmn_dr.x*ckmn_dr.x + ckmn_dr.y*ckmn_dr.y + ckmn_dr.z*ckmn_dr.z);
			//c_jkr = exp(-cj*ka*ckmn_r)/ckmn_r
			//c_jkr = dcmplx( cos(ka*ckmn_r), -sin(ka*ckmn_r) )
			//c_jkr = c_jkr/ckmn_r
			c_jkr = Z_MAKE(cos(ka*ckmn_r) / ckmn_r, -sin(ka*ckmn_r) / ckmn_r);
			//cvval1 =  cvval1 + dcmplx(ckmn_s/dl,0.d0)*c_jkr	!ifun=1	
			//csval1 = csval1 + c_jkr							!ifun=3	
			cvval1 = Z_ADD(Z_MUL(Z_MAKE(ckmn_s / dl, (double)0.0e0), c_jkr), cvval1);
			csval1 = Z_ADD(csval1, c_jkr);

			cvval1 = Z_MUL(Z_MAKE(dcqg_b*(double)0.5e0, (double)0.0e0), cvval1);
			csval1 = Z_MUL(Z_MAKE(dcqg_b*(double)0.5e0, (double)0.0e0), csval1);

			//tvecl1= dc%x*lmc%x+dc%y*lmc%y+dc%z*lmc%z       ! obliczane !lm+ * 1n+ dla k=1, obliczane !lm- * 1n+ dla k=1
			tvecl1 = dc.x*lmc.x + dc.y*lmc.y + dc.z*lmc.z;
			//csval1 = csval1/dl
			csval1 = Z_DIV(csval1, Z_MAKE(dl, (double)0.0e0));
		}
		// ------------ Zamiast prcodeury cwire1 -------------------


		// ------------ Zamiast prcodeury cwire2 -------------------
		ns = ips[k + kk + j + nsmax];
		nlr = ise[ns - 1 + nsmax];
		an = rad[ns - 1];
		double3 rnr = make_double3((double)sgnx*(double)x[nlr - 1], (double)sgny*(double)y[nlr - 1], (double)sgnz*(double)z[nlr - 1]);
		rn = make_double3((double)sgnx*(double)x[np - 1], (double)sgny*(double)y[np - 1], (double)sgnz*(double)z[np - 1]);
		rd = make_double3(rnr.x - rn.x, rnr.y - rn.y, rnr.z - rn.z);
		dl = sqrt(rd.x*rd.x + rd.y*rd.y + rd.z*rd.z);
		//dlh = (double)0.5e0*dl;
		dc = make_double3(rd.x / dl, rd.y / dl, rd.z / dl);

		ise_ms1 = ise[ms - 1];
		ise_ms2 = ise[ms - 1 + nsmax];
		ise_ns1 = ise[ns - 1];
		ise_ns2 = ise[ns - 1 + nsmax];
		//z_ise_ms1 = z[ise[ms-1]];
		//z_ise_ms2 = z[ise[ms-1+nsmax]];
		z_ise_ms1 = z[ise[ms - 1] - 1];
		z_ise_ms2 = z[ise[ms - 1 + nsmax] - 1];

		if (war_gpu_pec_im(ise_ms1, ise_ms2, ise_ns1, ise_ns2, z_ise_ms1, z_ise_ms2)) {
			//call dcqg(0.d0,dlh,cknn,cv1,4)
			//nord=4		    
			// --------------- cv1 ------------------
			cvval2 = Z_MAKE((double)0.0e0, (double)0.0e0);
			csval2 = Z_MAKE((double)0.0e0, (double)0.0e0);
		}
		else {
			//call dcqg(0.d0,dl,ckmn,cv,2)
			//nord=2
			double dcqg_a = (double)0.5e0*dl;
			double dcqg_b = dl;
			double dcqg_c = (double)0.288675134594812882e0 * dcqg_b;
			//ifun=2
			//y=b*0.5d0*(f(a+c)+f(a-c))
			//dla s=a+c
			double ckmn_s = dcqg_a + dcqg_c;
			double3 ckmn_dr = make_double3(rmc.x - rn.x - ckmn_s*dc.x, rmc.y - rn.y - ckmn_s*dc.y, rmc.z - rn.z - ckmn_s*dc.z);
			double ckmn_r = sqrt(an*an + ckmn_dr.x*ckmn_dr.x + ckmn_dr.y*ckmn_dr.y + ckmn_dr.z*ckmn_dr.z);
			//stosujemy wz�r Eulera exp(-ja)=cos(a)-jsin(a) => a=ka*r 
			//c_jkr = cdexp(-cj*ka*ckmn_r)/ckmn_r
			cuDoubleComplex c_jkr = Z_MAKE(cos(ka*ckmn_r) / ckmn_r, -sin(ka*ckmn_r) / ckmn_r);
			//cvval2 = dcmplx((dl-ckmn_s)/dl,0.d0)*c_jkr	!ifun=2
			cvval2 = Z_MUL(Z_MAKE((dl - ckmn_s) / dl, (double)0.0e0), c_jkr);
			//csval1 = c_jkr	                            !ifun=3
			csval2 = c_jkr;
			//dla s=a-c
			ckmn_s = dcqg_a - dcqg_c;
			ckmn_dr = make_double3(rmc.x - rn.x - ckmn_s*dc.x, rmc.y - rn.y - ckmn_s*dc.y, rmc.z - rn.z - ckmn_s*dc.z);
			ckmn_r = sqrt(an*an + ckmn_dr.x*ckmn_dr.x + ckmn_dr.y*ckmn_dr.y + ckmn_dr.z*ckmn_dr.z);
			//c_jkr = exp(-cj*ka*ckmn_r)/ckmn_r
			//c_jkr = dcmplx( cos(ka*ckmn_r), -sin(ka*ckmn_r) )
			//c_jkr = c_jkr/ckmn_r
			c_jkr = Z_MAKE(cos(ka*ckmn_r) / ckmn_r, -sin(ka*ckmn_r) / ckmn_r);
			//cvval2 =  cvval2 + dcmplx((dl-ckmn_s)/dl,0.d0)*c_jkr	!ifun=2	
			//csval2 = csval2 + c_jkr								!ifun=3		
			cvval2 = Z_ADD(Z_MUL(Z_MAKE((dl - ckmn_s) / dl, (double)0.0e0), c_jkr), cvval2);
			csval2 = Z_ADD(csval2, c_jkr);

			cvval2 = Z_MUL(Z_MAKE(dcqg_b*(double)0.5e0, (double)0.0e0), cvval2);
			csval2 = Z_MUL(Z_MAKE(dcqg_b*(double)0.5e0, (double)0.0e0), csval2);

			//tvecl2 = dc%x*lmc%x+dc%y*lmc%y+dc%z*lmc%z      ! obliczane !lm+ * 1n- dla k=1,	obliczane  !lm- * 1n- dla k=2
			tvecl2 = dc.x*lmc.x + dc.y*lmc.y + dc.z*lmc.z;
			//csval2 = csval2/dl
			csval2 = Z_DIV(csval2, Z_MAKE(dl, (double)0.0e0));
		}
		// ------------ Zamiast prcodeury cwire2 -------------------

		//cvval = sgnenv*dcmplx(tvecl1,0.d0)*cvval1 + sgnenv*dcmplx(tvecl2,0.d0)*cvval2    !A(rm+)	dla k=1, A(rm-) dla k=2		  
		//csval = csval1*sgnenv - sgnenv*csval2											!O(rm+) dla k=1, !O(rm-) dla k=2
		cuDoubleComplex cvval = Z_MUL(Z_MAKE(tvecl1, (double)0.0e0), cvval1);
		cvval = Z_ADD(Z_MUL(Z_MAKE(tvecl2, (double)0.0e0), cvval2), cvval);
		cvval = Z_MUL(Z_MAKE((double)sgnenv, (double)0.0e0), cvval);
		cuDoubleComplex csval = Z_MUL(Z_MAKE((double)sgnenv, (double)0.0e0), Z_SUB(csval1, csval2));


		//cz_ij = cz_ij+ka*ka*cvval + znak * csval
		cuDoubleComplex cz_ij = Z_MUL(Z_MAKE((double)znak, (double)0.0e0), csval);
		cz_ij = Z_ADD(Z_MUL(Z_MAKE(ka*ka, (double)0.0e0), cvval), cz_ij);

		//cz_ij = cz_ij*cj*eta0/(4.0e0*pi*ka)
		cz_ij = Z_MUL(cz_ij, Z_MAKE(eta0 / ((double)4.0e0*pi*ka), (double)0.0e0));
		cz_ij = Z_MUL(cz_ij, Z_ONE_IMG);



		//cz_d[i_m+j_n*n] = Z_MAKE(10.0f,0.0f);
		//cz_d[i_m+j_n*n] = Z_MUL ( cz_d[i_m+j_n*n] , (Z_MAKE(10.0f,0.0f));
		//cz_d[i_m+j_n*n] = Z_MAKE((double)i_m,(double)j_n);
		//cz_d[i_m+j_n*n]=make_cuDoubleComplex( (i_m),(j_n) );
		//cz_d[ i+j*n ] =  Z_MAKE(float(i*n+j),0.0f) ;
		//cz[ i+j*n ] = Z_MAKE(x[mp-1],z[mp-1]);
		//cz[ i+j*n ] = Z_MAKE(dc.x,dc.z);
		//cz[ i+j*n ] = cknn_cbnd1;
		//cz[ i+j*n ] = Z_MAKE(cknn_res,cknn_bnd2);
		//cz[ i+j*n ] = cz_ij;
		cz[i + (j + kk)*n] = Z_ADD(cz[i + (j + kk)*n], cz_ij);
		//cz[ i+j*n ] =  Z_ADD(c1.z,Z_MAKE(float(i*n+j),0.0f)) ;
	}
}

//-------------------------------------------------------------------------------

static inline int iabs(int x)
{
	if (x >= 0) {
		return x;
	}
	else {
		return -x;
	}
}

//-------------------------------------------------------------------------------

static inline double delick_cpu(double bet) {

	double a0 = (double)1.38629436112e0;
	double a1 = (double)0.09666344259e0;
	double a2 = (double)0.03590092383e0;
	double a3 = (double)0.03742563713e0;
	double a4 = (double)0.01451196212e0;
	double b0 = (double)0.5e0;
	double b1 = (double)0.12498593597e0;
	double b2 = (double)0.06880248576e0;
	double b3 = (double)0.03328355346e0;
	double b4 = (double)0.00441787012e0;
	double res, a, b;
	double am1 = (double)1.0e0 - bet * bet;
	double am12, am13, am14;

	a = a0 + a1 * am1;
	b = b0 + b1 * am1;

	if (am1 >= (double)1.0e-18) {
		am12 = am1 * am1;
		a = a + a2 * am12;
		b = b + b2 * am12;
		if (am1 >= (double)1.0e-12) {
			am13 = am12 * am1;
			a = a + a3 * am13;
			b = b + b3 * am13;
			if (am1 >= (double)1.0e-9) {
				am14 = am13 * am1;
				a = a + a4 * am14;
				b = b + b4 * am14;
			}
		}
	}

	res = a - b * log(am1);
	return res;
}

//-------------------------------------------------------------------------------

static inline bool war_cpu_pec_re(int ise_ms1, int ise_ms2, int ise_ns1, int ise_ns2)
{
	if ((ise_ns1 == ise_ms1 && ise_ns2 == ise_ms2) || (ise_ns1 == ise_ms2 && ise_ns2 == ise_ms1))
	{
		return true;
	}
	else {
		return false;
	}
}

//-------------------------------------------------------------------------------

static inline bool war_cpu_pec_im(int ise_ms1, int ise_ms2, int ise_ns1, int ise_ns2, float z_ise_ms1, float z_ise_ms2)
{
	if (((ise_ns1 == ise_ms1 && ise_ns2 == ise_ms2) && (z_ise_ms1 == 0.0f && z_ise_ms2 == 0.0f))
		|| ((ise_ns1 == ise_ms2 && ise_ns2 == ise_ms1) && (z_ise_ms1 == 0.0f && z_ise_ms2 == 0.0f)))
	{
		return true;
	}
	else {
		return false;
	}
}

//-------------------------------------------------------------------------------

void cpu_cknn_cwire1_re(clDoubleComplex* cz, float* x, float* y, float* z, float* rad,
	int* ipc, int* ips, int* ise, int nsmax, int n, double ka, double eta0,
	int realSubSliceSize, int k, int kk)
{

	double pi = 3.14159265358979323846e0;
	clDoubleComplex cvval1, csval1;
	//double tvecl1;


	for (int k1 = 0; k1 < 2; k1++)
	{
		for (int i = 0; i < n; i++)
		{
			int mp = iabs(ipc[i]);
			//int ms = ips[i];
			int ms = ips[i + k1 * nsmax];
			int mlr = ise[ms - 1 + k1 * nsmax];
			//int mlr = ise[ms-1];
			double rmp_x = (double)x[mp - 1];
			double rmp_y = (double)y[mp - 1];
			double rmp_z = (double)z[mp - 1];
			double rmlr_x = (double)x[mlr - 1];
			double rmlr_y = (double)y[mlr - 1];
			double rmlr_z = (double)z[mlr - 1];
			double rmc_x = (double)0.5e0 * (rmlr_x + rmp_x);
			double rmc_y = (double)0.5e0 * (rmlr_y + rmp_y);
			double rmc_z = (double)0.5e0 * (rmlr_z + rmp_z);
			//double_3 rmp = make_double_3( (double)x[mp-1], (double)y[mp-1], (double)z[mp-1] );
			//double_3 rmlr = make_double_3( (double)x[mlr-1], (double)y[mlr-1], (double)z[mlr-1] );
			//double_3 rmc = make_double_3( (double)0.5e0*(rmlr.x+rmp.x), (double)0.5e0*(rmlr.y+rmp.y), (double)0.5e0*(rmlr.z+rmp.z) );
			//float znak = -1.0f;
			float znak = 2.0f * (float)k1 - 1.0f;
			double lmc_x = (double)znak * (rmc_x - rmp_x);
			double lmc_y = (double)znak * (rmc_y - rmp_y);
			double lmc_z = (double)znak * (rmc_z - rmp_z);
			//double_3 lmc = make_double_3( (double)znak*(rmc.x-rmp.x), (double)znak*(rmc.y-rmp.y), (double)znak*(rmc.z-rmp.z) );

			for (int j = 0; j < realSubSliceSize; j++)
			{
				// ------------ Zamiast prcodeury cwire1 -------------------
				// Obliczanie calki na Sn+  pot.wek.
				int np = iabs(ipc[k + kk + j]);
				int ns = ips[k + kk + j];
				//int np = iabs(ipc[j]);
				//int ns = ips[j];
				double rnp_x = (double)x[np - 1];
				double rnp_y = (double)y[np - 1];
				double rnp_z = (double)z[np - 1];
				//double_3 rnp = make_double_3( (double)x[np-1],(double)y[np-1],(double)z[np-1] );
				int nlr = ise[ns - 1];
				double an = rad[ns - 1];
				double rn_x = (double)x[nlr - 1];
				double rn_y = (double)y[nlr - 1];
				double rn_z = (double)z[nlr - 1];
				double rd_x = rnp_x - rn_x;
				double rd_y = rnp_y - rn_y;
				double rd_z = rnp_z - rn_z;
				//double_3 rn = make_double_3( (double)x[nlr-1],(double)y[nlr-1],(double)z[nlr-1] );
				//double_3 rd = make_double_3( rnp.x-rn.x, rnp.y-rn.y, rnp.z-rn.z );
				double dl = sqrt(rd_x * rd_x + rd_y * rd_y + rd_z * rd_z);
				double dlh = (double)0.5e0 * dl;
				double dc_x = rd_x / dl;
				double dc_y = rd_y / dl;
				double dc_z = rd_z / dl;
				//double_3 dc = make_double_3( rd.x/dl, rd.y/dl, rd.z/dl );

				int ise_ms1 = ise[ms - 1];
				int ise_ms2 = ise[ms - 1 + nsmax];
				int ise_ns1 = ise[ns - 1];
				int ise_ns2 = ise[ns - 1 + nsmax];

				if (war_cpu_pec_re(ise_ms1, ise_ms2, ise_ns1, ise_ns2)) {
					// --------------- cv1 ------------------
					double dcqg_a = (double)0.5e0 * dlh;
					double dcqg_b = dlh;
					double dcqg_c = (double)0.43056815579702629e0 * dcqg_b;
					//s=a+c
					double cknn_s = dcqg_a + dcqg_c;
					// ---- cknn ----
					double cknn_r = (double)0.5e0 * dl - cknn_s;
					double cknn_rr = sqrt(an * an + cknn_r * cknn_r);
					clDoubleComplex cknn_cbnd1 = Z_MAKE((cos(ka * cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka * cknn_rr) / cknn_rr);
					double cknn_bet = (an + an) / sqrt((double)4.0e0 * an * an + cknn_r * cknn_r);
					double cknn_res = delick_cpu(cknn_bet);
					double cknn_bnd2 = (cknn_bet * cknn_res + log(fabs(cknn_r / ((double)8.e0 * an)))) / (pi * an);
					//ifun=1
					clDoubleComplex cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));   // 1 segment n-tej f.bazowej
					clDoubleComplex cvv1a = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
					//ifun=3
					clDoubleComplex csv1a = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));  // pot. scalarny
					// ---- cknn ----
					//s=a-c
					cknn_s = dcqg_a - dcqg_c;
					// ---- cknn ----
					cknn_r = (double)0.5e0 * dl - cknn_s;
					cknn_rr = sqrt(an * an + cknn_r * cknn_r);
					cknn_cbnd1 = Z_MAKE((cos(ka * cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka * cknn_rr) / cknn_rr);
					cknn_bet = (an + an) / sqrt((double)4.0e0 * an * an + cknn_r * cknn_r);
					cknn_res = delick_cpu(cknn_bet);
					cknn_bnd2 = (cknn_bet * cknn_res + log(fabs(cknn_r / ((double)8.e0 * an)))) / (pi * an);
					//ifun=1
					cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));   // 1 segment n-tej f.bazowej
					cv = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
					cvv1a = Z_MUL(Z_MAKE((double)0.17392742256872693e0, (double)0.0e0), Z_ADD(cv, cvv1a));
					//ifun=3
					cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));   // pot. scalarny
					csv1a = Z_MUL(Z_MAKE((double)0.17392742256872693e0, (double)0.0e0), Z_ADD(cv, csv1a));
					// ---- cknn ----
					dcqg_c = (double)0.16999052179242813e0 * dcqg_b;
					//s=a+c
					cknn_s = dcqg_a + dcqg_c;
					// ---- cknn ----
					cknn_r = (double)0.5e0 * dl - cknn_s;
					cknn_rr = sqrt(an * an + cknn_r * cknn_r);
					cknn_cbnd1 = Z_MAKE((cos(ka * cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka * cknn_rr) / cknn_rr);
					cknn_bet = (an + an) / sqrt((double)4.0e0 * an * an + cknn_r * cknn_r);
					cknn_res = delick_cpu(cknn_bet);
					cknn_bnd2 = (cknn_bet * cknn_res + log(fabs(cknn_r / ((double)8.e0 * an)))) / (pi * an);
					//ifun=1
					cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));   // 1 segment n-tej f.bazowej
					clDoubleComplex cvv1b = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
					//ifun=3
					clDoubleComplex csv1b = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
					// ---- cknn ----
					//s=a-c
					cknn_s = dcqg_a - dcqg_c;
					// ---- cknn ----
					cknn_r = (double)0.5e0 * dl - cknn_s;
					cknn_rr = sqrt(an * an + cknn_r * cknn_r);
					cknn_cbnd1 = Z_MAKE((cos(ka * cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka * cknn_rr) / cknn_rr);
					cknn_bet = (an + an) / sqrt((double)4.0e0 * an * an + cknn_r * cknn_r);
					cknn_res = delick_cpu(cknn_bet);
					cknn_bnd2 = (cknn_bet * cknn_res + log(fabs(cknn_r / ((double)8.e0 * an)))) / (pi * an);
					//ifun=1
					cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));   // !1 segment n-tej f.bazowej
					cv = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
					cvv1b = Z_MUL(Z_MAKE((double)0.32607257743127307e0, (double)0.0e0), Z_ADD(cv, cvv1b));
					//ifun=3
					cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));   // pot. scalarny
					csv1b = Z_MUL(Z_MAKE((double)0.32607257743127307e0, (double)0.0e0), Z_ADD(cv, csv1b));
					// ---- cknn ----
					clDoubleComplex cv1v = Z_MUL(Z_MAKE(dcqg_b, (double)0.0e0), Z_ADD(cvv1a, cvv1b));
					clDoubleComplex cv1s = Z_MUL(Z_MAKE(dcqg_b, (double)0.0e0), Z_ADD(csv1a, csv1b));
					// --------------- cv1 ------------------

					// --------------- cv2 ------------------
					dcqg_a = (double)0.5e0 * (dlh + dl);
					dcqg_b = dl - dlh;
					dcqg_c = (double)0.43056815579702629e0 * dcqg_b;
					//s=a+c
					cknn_s = dcqg_a + dcqg_c;
					// ---- cknn ----
					cknn_r = (double)0.5e0 * dl - cknn_s;
					cknn_rr = sqrt(an * an + cknn_r * cknn_r);
					cknn_cbnd1 = Z_MAKE((cos(ka * cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka * cknn_rr) / cknn_rr);
					cknn_bet = (an + an) / sqrt((double)4.0e0 * an * an + cknn_r * cknn_r);
					cknn_res = delick_cpu(cknn_bet);
					cknn_bnd2 = (cknn_bet * cknn_res + log(fabs(cknn_r / ((double)8.e0 * an)))) / (pi * an);
					//ifun=1
					cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));   // 1 segment n-tej f.bazowej
					cvv1a = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
					//ifun=3
					csv1a = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));   // pot. scalarny
					// ---- cknn ----
					//s=a-c
					cknn_s = dcqg_a - dcqg_c;
					// ---- cknn ----
					cknn_r = (double)0.5e0 * dl - cknn_s;
					cknn_rr = sqrt(an * an + cknn_r * cknn_r);
					cknn_cbnd1 = Z_MAKE((cos(ka * cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka * cknn_rr) / cknn_rr);
					cknn_bet = (an + an) / sqrt((double)4.0e0 * an * an + cknn_r * cknn_r);
					cknn_res = delick_cpu(cknn_bet);
					cknn_bnd2 = (cknn_bet * cknn_res + log(fabs(cknn_r / ((double)8.e0 * an)))) / (pi * an);
					//ifun=1
					cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));   // 1 segment n-tej f.bazowej
					cv = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
					cvv1a = Z_MUL(Z_MAKE((double)0.17392742256872693e0, (double)0.0e0), Z_ADD(cv, cvv1a));
					//ifun=3
					cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));   // pot. scalarny
					csv1a = Z_MUL(Z_MAKE((double)0.17392742256872693e0, (double)0.0e0), Z_ADD(cv, csv1a));
					// ---- cknn ----
					dcqg_c = (double)0.16999052179242813e0 * dcqg_b;
					//s=a+c
					cknn_s = dcqg_a + dcqg_c;
					// ---- cknn ----
					cknn_r = (double)0.5e0 * dl - cknn_s;
					cknn_rr = sqrt(an * an + cknn_r * cknn_r);
					cknn_cbnd1 = Z_MAKE((cos(ka * cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka * cknn_rr) / cknn_rr);
					cknn_bet = (an + an) / sqrt((double)4.0e0 * an * an + cknn_r * cknn_r);
					cknn_res = delick_cpu(cknn_bet);
					cknn_bnd2 = (cknn_bet * cknn_res + log(fabs(cknn_r / ((double)8.e0 * an)))) / (pi * an);
					//ifun=1
					cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));   // 1 segment n-tej f.bazowej
					cvv1b = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
					//ifun=3
					csv1b = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));   // pot. scalarny
					// ---- cknn ----
					//s=a-c
					cknn_s = dcqg_a - dcqg_c;
					// ---- cknn ----
					cknn_r = (double)0.5e0 * dl - cknn_s;
					cknn_rr = sqrt(an * an + cknn_r * cknn_r);
					cknn_cbnd1 = Z_MAKE((cos(ka * cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka * cknn_rr) / cknn_rr);
					cknn_bet = (an + an) / sqrt((double)4.0e0 * an * an + cknn_r * cknn_r);
					cknn_res = delick_cpu(cknn_bet);
					cknn_bnd2 = (cknn_bet * cknn_res + log(fabs(cknn_r / ((double)8.e0 * an)))) / (pi * an);
					//ifun=1
					cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));   // 1 segment n-tej f.bazowej
					cv = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
					cvv1b = Z_MUL(Z_MAKE((double)0.32607257743127307e0, (double)0.0e0), Z_ADD(cv, cvv1b));
					//ifun=3
					cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));   // pot. scalarny
					csv1b = Z_MUL(Z_MAKE((double)0.32607257743127307e0, (double)0.0e0), Z_ADD(cv, csv1b));
					// ---- cknn ----
					clDoubleComplex cv2v = Z_MUL(Z_MAKE(dcqg_b, (double)0.0e0), Z_ADD(cvv1a, cvv1b));
					clDoubleComplex cv2s = Z_MUL(Z_MAKE(dcqg_b, (double)0.0e0), Z_ADD(csv1a, csv1b));
					// --------------- cv2 ------------------

					double res_cs = dl * ((double)1.0e0 + log((double)16.0e0 * an / dl)) / (pi * an);
					double res_cv = (double)0.5e0 * res_cs;
					cvval1 = Z_ADD(cv2v, Z_MAKE(res_cv, (double)0.0e0));
					cvval1 = Z_ADD(cv1v, cvval1);
					csval1 = Z_ADD(cv2s, Z_MAKE(res_cs, (double)0.0e0));
					csval1 = Z_ADD(cv1s, csval1);
				}
				else {
					cvval1 = Z_ZERO;
					csval1 = Z_ZERO;
				}
				// ------------ Zamiast prcodeury cwire1 -------------------

				double tvecl1 = dc_x * lmc_x + dc_y * lmc_y + dc_z * lmc_z;   // obliczane !lm+ * 1n+ dla k=1, obliczane !lm- * 1n+ dla k=1
				csval1 = Z_DIV(csval1, Z_MAKE(dl, (double)0.0e0));
				clDoubleComplex cvval = Z_MUL(cvval1, Z_MAKE(tvecl1, (double)0.0e0));
				clDoubleComplex csval = csval1;
				//cz_ij = cz_ij+ka*ka*cvval + znak * csval
				clDoubleComplex cz_ij = Z_MUL(Z_MAKE((double)znak, (double)0.0e0), csval);
				cz_ij = Z_ADD(Z_MUL(Z_MAKE(ka * ka, (double)0.0e0), cvval), cz_ij);
				//cz_ij = cz_ij*cj*eta0/(4.0e0*pi*ka)
				cz_ij = Z_MUL(cz_ij, Z_MAKE(eta0 / ((double)4.0e0 * pi * ka), (double)0.0e0));
				cz_ij = Z_MUL(cz_ij, Z_ONE_IMG);

				//cz[ i+(j+kk)*n ] = Z_ADD(cz[i + (j + kk)*n], cz_ij);


				//clDoubleComplex cz_ij = Z_MAKE((double)0.0e0,(double)0.0e0);  //do testow
				//clDoubleComplex cz_ij	= cvval;
				//clDoubleComplex cz_ij = Z_MAKE((double)i,(double)j);
				//clDoubleComplex cz_ij = Z_MAKE(tvecl1,tvecl2);
				//clDoubleComplex cz_ij = Z_MAKE(rmp_z,rmlr_z);
				//cz[i + (j + kk)*n] = Z_MAKE(rmc_z,lmc_z);
				//cz[i + (j + kk)*n] = Z_MAKE(tvecl1,tvecl2);
				//cz[i + (j + kk)*n] = Z_MAKE((double)i,(double)j);
				//clDoubleComplex cz_ij = Z_MAKE((double)i,(double)j);
				cz[i + (j + kk) * n] = Z_ADD(cz[i + (j + kk) * n], cz_ij);
			}
		}
	}


	//double sin_v = sin(pi/4);
	//printf("Got sin value: %.14lf \n", sin_v);

}

//-------------------------------------------------------------------------------

void cpu_cknn_cwire2_re(clDoubleComplex* cz, float* x, float* y, float* z, float* rad,
	int* ipc, int* ips, int* ise, int nsmax, int n, double ka, double eta0,
	int realSubSliceSize, int k, int kk)
{

	double pi = 3.14159265358979323846e0;
	clDoubleComplex cvval2, csval2;
	//double tvecl1;


	for (int k1 = 0; k1 < 2; k1++)
	{
		for (int i = 0; i < n; i++)
		{
			int mp = iabs(ipc[i]);
			//int ms = ips[i];
			int ms = ips[i + k1 * nsmax];
			int mlr = ise[ms - 1 + k1 * nsmax];
			//int mlr = ise[ms-1];
			double rmp_x = (double)x[mp - 1];
			double rmp_y = (double)y[mp - 1];
			double rmp_z = (double)z[mp - 1];
			double rmlr_x = (double)x[mlr - 1];
			double rmlr_y = (double)y[mlr - 1];
			double rmlr_z = (double)z[mlr - 1];
			double rmc_x = (double)0.5e0 * (rmlr_x + rmp_x);
			double rmc_y = (double)0.5e0 * (rmlr_y + rmp_y);
			double rmc_z = (double)0.5e0 * (rmlr_z + rmp_z);
			//double_3 rmp = make_double_3( (double)x[mp-1], (double)y[mp-1], (double)z[mp-1] );
			//double_3 rmlr = make_double_3( (double)x[mlr-1], (double)y[mlr-1], (double)z[mlr-1] );
			//double_3 rmc = make_double_3( (double)0.5e0*(rmlr.x+rmp.x), (double)0.5e0*(rmlr.y+rmp.y), (double)0.5e0*(rmlr.z+rmp.z) );
			//float znak = -1.0f;
			float znak = 2.0f * (float)k1 - 1.0f;
			double lmc_x = (double)znak * (rmc_x - rmp_x);
			double lmc_y = (double)znak * (rmc_y - rmp_y);
			double lmc_z = (double)znak * (rmc_z - rmp_z);
			//double_3 lmc = make_double_3( (double)znak*(rmc.x-rmp.x), (double)znak*(rmc.y-rmp.y), (double)znak*(rmc.z-rmp.z) );

			for (int j = 0; j < realSubSliceSize; j++)
			{
				// ------------ Zamiast prcodeury cwire2 -------------------
				// Obliczanie calki na Sn+  pot.wek.
				int np = iabs(ipc[k + kk + j]);
				int ns = ips[k + kk + j + nsmax];
				//int np = iabs(ipc[j]);
				//int ns = ips[j];
				int nlr = ise[ns - 1 + nsmax];
				double an = rad[ns - 1];
				//double_3 rn = make_double_3( (double)x[np-1],(double)y[np-1],(double)z[np-1] );
				double rn_x = (double)x[np - 1];
				double rn_y = (double)y[np - 1];
				double rn_z = (double)z[np - 1];
				double rnr_x = (double)x[nlr - 1];
				double rnr_y = (double)y[nlr - 1];
				double rnr_z = (double)z[nlr - 1];
				double rd_x = rnr_x - rn_x;
				double rd_y = rnr_y - rn_y;
				double rd_z = rnr_z - rn_z;
				//double_3 rn = make_double_3( (double)x[nlr-1],(double)y[nlr-1],(double)z[nlr-1] );
				//double_3 rd = make_double_3( rnp.x-rn.x, rnp.y-rn.y, rnp.z-rn.z );
				double dl = sqrt(rd_x * rd_x + rd_y * rd_y + rd_z * rd_z);
				double dlh = (double)0.5e0 * dl;
				double dc_x = rd_x / dl;
				double dc_y = rd_y / dl;
				double dc_z = rd_z / dl;
				//double_3 dc = make_double_3( rd.x/dl, rd.y/dl, rd.z/dl );

				int ise_ms1 = ise[ms - 1];
				int ise_ms2 = ise[ms - 1 + nsmax];
				int ise_ns1 = ise[ns - 1];
				int ise_ns2 = ise[ns - 1 + nsmax];

				if (war_cpu_pec_re(ise_ms1, ise_ms2, ise_ns1, ise_ns2)) {
					// --------------- cv1 ------------------
					double dcqg_a = (double)0.5e0 * dlh;
					double dcqg_b = dlh;
					double dcqg_c = (double)0.43056815579702629e0 * dcqg_b;
					//s=a+c
					double cknn_s = dcqg_a + dcqg_c;
					// ---- cknn ----
					double cknn_r = (double)0.5e0 * dl - cknn_s;
					double cknn_rr = sqrt(an * an + cknn_r * cknn_r);
					clDoubleComplex cknn_cbnd1 = Z_MAKE((cos(ka * cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka * cknn_rr) / cknn_rr);
					double cknn_bet = (an + an) / sqrt((double)4.0e0 * an * an + cknn_r * cknn_r);
					double cknn_res = delick_cpu(cknn_bet);
					double cknn_bnd2 = (cknn_bet * cknn_res + log(fabs(cknn_r / ((double)8.e0 * an)))) / (pi * an);
					//ifun=1
					clDoubleComplex cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));   // 1 segment n-tej f.bazowej
					clDoubleComplex cvv1a = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
					//ifun=3
					clDoubleComplex csv1a = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));  // pot. scalarny
					// ---- cknn ----
					//s=a-c
					cknn_s = dcqg_a - dcqg_c;
					// ---- cknn ----
					cknn_r = (double)0.5e0 * dl - cknn_s;
					cknn_rr = sqrt(an * an + cknn_r * cknn_r);
					cknn_cbnd1 = Z_MAKE((cos(ka * cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka * cknn_rr) / cknn_rr);
					cknn_bet = (an + an) / sqrt((double)4.0e0 * an * an + cknn_r * cknn_r);
					cknn_res = delick_cpu(cknn_bet);
					cknn_bnd2 = (cknn_bet * cknn_res + log(fabs(cknn_r / ((double)8.e0 * an)))) / (pi * an);
					//ifun=1
					cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));   // 1 segment n-tej f.bazowej
					cv = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
					cvv1a = Z_MUL(Z_MAKE((double)0.17392742256872693e0, (double)0.0e0), Z_ADD(cv, cvv1a));
					//ifun=3
					cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));   // pot. scalarny
					csv1a = Z_MUL(Z_MAKE((double)0.17392742256872693e0, (double)0.0e0), Z_ADD(cv, csv1a));
					// ---- cknn ----
					dcqg_c = (double)0.16999052179242813e0 * dcqg_b;
					//s=a+c
					cknn_s = dcqg_a + dcqg_c;
					// ---- cknn ----
					cknn_r = (double)0.5e0 * dl - cknn_s;
					cknn_rr = sqrt(an * an + cknn_r * cknn_r);
					cknn_cbnd1 = Z_MAKE((cos(ka * cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka * cknn_rr) / cknn_rr);
					cknn_bet = (an + an) / sqrt((double)4.0e0 * an * an + cknn_r * cknn_r);
					cknn_res = delick_cpu(cknn_bet);
					cknn_bnd2 = (cknn_bet * cknn_res + log(fabs(cknn_r / ((double)8.e0 * an)))) / (pi * an);
					//ifun=1
					cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));   // 1 segment n-tej f.bazowej
					clDoubleComplex cvv1b = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
					//ifun=3
					clDoubleComplex csv1b = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));
					// ---- cknn ----
					//s=a-c
					cknn_s = dcqg_a - dcqg_c;
					// ---- cknn ----
					cknn_r = (double)0.5e0 * dl - cknn_s;
					cknn_rr = sqrt(an * an + cknn_r * cknn_r);
					cknn_cbnd1 = Z_MAKE((cos(ka * cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka * cknn_rr) / cknn_rr);
					cknn_bet = (an + an) / sqrt((double)4.0e0 * an * an + cknn_r * cknn_r);
					cknn_res = delick_cpu(cknn_bet);
					cknn_bnd2 = (cknn_bet * cknn_res + log(fabs(cknn_r / ((double)8.e0 * an)))) / (pi * an);
					//ifun=1
					cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));   // !1 segment n-tej f.bazowej
					cv = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
					cvv1b = Z_MUL(Z_MAKE((double)0.32607257743127307e0, (double)0.0e0), Z_ADD(cv, cvv1b));
					//ifun=3
					cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));   // pot. scalarny
					csv1b = Z_MUL(Z_MAKE((double)0.32607257743127307e0, (double)0.0e0), Z_ADD(cv, csv1b));
					// ---- cknn ----
					clDoubleComplex cv1v = Z_MUL(Z_MAKE(dcqg_b, (double)0.0e0), Z_ADD(cvv1a, cvv1b));
					clDoubleComplex cv1s = Z_MUL(Z_MAKE(dcqg_b, (double)0.0e0), Z_ADD(csv1a, csv1b));
					// --------------- cv1 ------------------

					// --------------- cv2 ------------------
					dcqg_a = (double)0.5e0 * (dlh + dl);
					dcqg_b = dl - dlh;
					dcqg_c = (double)0.43056815579702629e0 * dcqg_b;
					//s=a+c
					cknn_s = dcqg_a + dcqg_c;
					// ---- cknn ----
					cknn_r = (double)0.5e0 * dl - cknn_s;
					cknn_rr = sqrt(an * an + cknn_r * cknn_r);
					cknn_cbnd1 = Z_MAKE((cos(ka * cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka * cknn_rr) / cknn_rr);
					cknn_bet = (an + an) / sqrt((double)4.0e0 * an * an + cknn_r * cknn_r);
					cknn_res = delick_cpu(cknn_bet);
					cknn_bnd2 = (cknn_bet * cknn_res + log(fabs(cknn_r / ((double)8.e0 * an)))) / (pi * an);
					//ifun=1
					cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));   // 1 segment n-tej f.bazowej
					cvv1a = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
					//ifun=3
					csv1a = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));   // pot. scalarny
					// ---- cknn ----
					//s=a-c
					cknn_s = dcqg_a - dcqg_c;
					// ---- cknn ----
					cknn_r = (double)0.5e0 * dl - cknn_s;
					cknn_rr = sqrt(an * an + cknn_r * cknn_r);
					cknn_cbnd1 = Z_MAKE((cos(ka * cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka * cknn_rr) / cknn_rr);
					cknn_bet = (an + an) / sqrt((double)4.0e0 * an * an + cknn_r * cknn_r);
					cknn_res = delick_cpu(cknn_bet);
					cknn_bnd2 = (cknn_bet * cknn_res + log(fabs(cknn_r / ((double)8.e0 * an)))) / (pi * an);
					//ifun=1
					cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));   // 1 segment n-tej f.bazowej
					cv = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
					cvv1a = Z_MUL(Z_MAKE((double)0.17392742256872693e0, (double)0.0e0), Z_ADD(cv, cvv1a));
					//ifun=3
					cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));   // pot. scalarny
					csv1a = Z_MUL(Z_MAKE((double)0.17392742256872693e0, (double)0.0e0), Z_ADD(cv, csv1a));
					// ---- cknn ----
					dcqg_c = (double)0.16999052179242813e0 * dcqg_b;
					//s=a+c
					cknn_s = dcqg_a + dcqg_c;
					// ---- cknn ----
					cknn_r = (double)0.5e0 * dl - cknn_s;
					cknn_rr = sqrt(an * an + cknn_r * cknn_r);
					cknn_cbnd1 = Z_MAKE((cos(ka * cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka * cknn_rr) / cknn_rr);
					cknn_bet = (an + an) / sqrt((double)4.0e0 * an * an + cknn_r * cknn_r);
					cknn_res = delick_cpu(cknn_bet);
					cknn_bnd2 = (cknn_bet * cknn_res + log(fabs(cknn_r / ((double)8.e0 * an)))) / (pi * an);
					//ifun=1
					cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));   // 1 segment n-tej f.bazowej
					cvv1b = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
					//ifun=3
					csv1b = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));   // pot. scalarny
					// ---- cknn ----
					//s=a-c
					cknn_s = dcqg_a - dcqg_c;
					// ---- cknn ----
					cknn_r = (double)0.5e0 * dl - cknn_s;
					cknn_rr = sqrt(an * an + cknn_r * cknn_r);
					cknn_cbnd1 = Z_MAKE((cos(ka * cknn_rr) - (double)1.0e0) / cknn_rr, -sin(ka * cknn_rr) / cknn_rr);
					cknn_bet = (an + an) / sqrt((double)4.0e0 * an * an + cknn_r * cknn_r);
					cknn_res = delick_cpu(cknn_bet);
					cknn_bnd2 = (cknn_bet * cknn_res + log(fabs(cknn_r / ((double)8.e0 * an)))) / (pi * an);
					//ifun=1
					cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));   // 1 segment n-tej f.bazowej
					cv = Z_MUL(Z_MAKE(cknn_s / dl, (double)0.0e0), cv);
					cvv1b = Z_MUL(Z_MAKE((double)0.32607257743127307e0, (double)0.0e0), Z_ADD(cv, cvv1b));
					//ifun=3
					cv = Z_ADD(cknn_cbnd1, Z_MAKE(cknn_bnd2, (double)0.0e0));   // pot. scalarny
					csv1b = Z_MUL(Z_MAKE((double)0.32607257743127307e0, (double)0.0e0), Z_ADD(cv, csv1b));
					// ---- cknn ----
					clDoubleComplex cv2v = Z_MUL(Z_MAKE(dcqg_b, (double)0.0e0), Z_ADD(cvv1a, cvv1b));
					clDoubleComplex cv2s = Z_MUL(Z_MAKE(dcqg_b, (double)0.0e0), Z_ADD(csv1a, csv1b));
					// --------------- cv2 ------------------

					double res_cs = dl * ((double)1.0e0 + log((double)16.0e0 * an / dl)) / (pi * an);
					double res_cv = (double)0.5e0 * res_cs;
					cvval2 = Z_ADD(cv2v, Z_MAKE(res_cv, (double)0.0e0));
					cvval2 = Z_ADD(cv1v, cvval2);
					csval2 = Z_ADD(cv2s, Z_MAKE(res_cs, (double)0.0e0));
					csval2 = Z_ADD(cv1s, csval2);
				}
				else {
					cvval2 = Z_ZERO;
					csval2 = Z_ZERO;
				}
				// ------------ Zamiast prcodeury cwire2 -------------------

				double tvecl2 = dc_x * lmc_x + dc_y * lmc_y + dc_z * lmc_z;   // obliczane !lm+ * 1n+ dla k=1, obliczane !lm- * 1n+ dla k=2
				csval2 = Z_DIV(csval2, Z_MAKE(dl, (double)0.0e0));
				clDoubleComplex cvval = Z_MUL(cvval2, Z_MAKE(tvecl2, (double)0.0e0));
				//csval = -csval2
				clDoubleComplex csval = Z_MUL(Z_MAKE((double)-1.0e0, (double)0.0e0), csval2);
				//cz_ij = cz_ij+ka*ka*cvval + znak * csval
				clDoubleComplex cz_ij = Z_MUL(Z_MAKE((double)znak, (double)0.0e0), csval);
				cz_ij = Z_ADD(Z_MUL(Z_MAKE(ka * ka, (double)0.0e0), cvval), cz_ij);
				//cz_ij = cz_ij*cj*eta0/(4.0e0*pi*ka)
				cz_ij = Z_MUL(cz_ij, Z_MAKE(eta0 / ((double)4.0e0 * pi * ka), (double)0.0e0));
				cz_ij = Z_MUL(cz_ij, Z_ONE_IMG);

				//cz[ i+(j+kk)*n ] = Z_ADD(cz[i + (j + kk)*n], cz_ij);
				//clDoubleComplex cz_ij = Z_MAKE((double)0.0e0,(double)0.0e0);  //do testow
				//clDoubleComplex cz_ij	= cvval;
				//clDoubleComplex cz_ij = Z_MAKE((double)i,(double)j);
				//clDoubleComplex cz_ij = Z_MAKE(tvecl1,tvecl2);
				//clDoubleComplex cz_ij = Z_MAKE(rmp_z,rmlr_z);
				//cz[i + (j + kk)*n] = Z_MAKE(rmc_z,lmc_z);
				//cz[i + (j + kk)*n] = Z_MAKE(tvecl1,tvecl2);
				//cz[i + (j + kk)*n] = Z_MAKE((double)i,(double)j);
				//clDoubleComplex cz_ij = Z_MAKE((double)i,(double)j);
				cz[i + (j + kk) * n] = Z_ADD(cz[i + (j + kk) * n], cz_ij);
			}
		}
	}


	//double sin_v = sin(pi/4);
	//printf("Got sin value: %.14lf \n", sin_v);

}

//----------------------------------------------------------------------

void zmatrix_gpu(cuDoubleComplex *cz_h, float *x_h, float *y_h, float *z_h, float *rad_h,
	int *ipc_h, int *ips_h, int *ise_h, int nsmax, double ka, int npls,
	int ienv, float sgnx, float sgny, float sgnz, float sgnenv, int rank)
{

	double pi = acos(-1.0);
	double eta0 = (double)120.*pi;
	cuDoubleComplex cjq = Z_MAKE(0.0, eta0 / (4.0*pi*ka));

	//std::cout << std::endl;
	//std::cout << "  --------------------------------------------------------------------------" << std::endl;
	//std::cout << "  Initializing hardware-accelerated zmatrix assembly phase of wire-grid MoM " << std::endl;
	//std::cout << "  --------------------------------------------------------------------------" << std::endl;

	float *x_d, *y_d, *z_d, *rad_d;
	int *ipc_d, *ips_d, *ise_d;
	cuDoubleComplex *cz_d;
	dim3 threads, grid;
	float elapsedTime_kernel;
	float elapsedTime;
	cudaEvent_t start, stop;
	cudaEvent_t start_kernel, stop_kernel;

	cudaMalloc((void**)&x_d, nsmax * sizeof(float));
	cudaMalloc((void**)&y_d, nsmax * sizeof(float));
	cudaMalloc((void**)&z_d, nsmax * sizeof(float));
	cudaMalloc((void**)&rad_d, nsmax * sizeof(float));
	cudaMalloc((void**)&ipc_d, nsmax * sizeof(int));
	cudaMalloc((void**)&ips_d, 2 * nsmax * sizeof(int));
	cudaMalloc((void**)&ise_d, 2 * nsmax * sizeof(int));
	cudaMalloc((void**)&cz_d, npls * npls * sizeof(cuDoubleComplex));

	//Open file where results are saved
	FILE *f = fopen("GPU time.dat", "w+");
	if (f == NULL)
	{
		printf("Error opening file!\n");
		exit(1);
	}

	//Prepare events
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventCreate(&start_kernel);
	cudaEventCreate(&stop_kernel);

	// Start record (data transfer and kernel execution)
	cudaEventRecord(start, 0);

	cudaMemcpy(x_d, x_h, nsmax * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(y_d, y_h, nsmax * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(z_d, z_h, nsmax * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(rad_d, rad_h, nsmax * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(ipc_d, ipc_h, nsmax * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(ips_d, ips_h, 2 * nsmax * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(ise_d, ise_h, 2 * nsmax * sizeof(int), cudaMemcpyHostToDevice);

	// Setup execution parameters
	threads = dim3(BLOCK_SIZE, BLOCK_SIZE);
	//threads = dim3(2*BLOCK_SIZE,BLOCK_SIZE);
	//grid = dim3( 1+(int)ceil(float(N-1)/float(threads.x)), 1+(int)ceil(float(N-1)/float(threads.y)) );
	grid = dim3(((npls - 1) / threads.x) + 1, ((npls - 1) / threads.y) + 1);

	//printf(" |\n");
	//printf("  --> set grid size %d  %d \n",grid.x,grid.y); 
	//printf("\n");

	// Start record (data transfer and kernel execution)
	//cudaEventRecord(start, 0);

	// Start record (kernel execution only)
	cudaEventRecord(start_kernel, 0);


	// Execute the kernels //  
	wim_gpu_cknn_cwire1_re1 <<< grid, threads >>>(cz_d, x_d, y_d, z_d, rad_d, ipc_d, ips_d, ise_d, nsmax, npls, ka, eta0, npls, 0, 0);   // device 1 when multi-GPU or device 1 and device 2 with another SubSlice level
	wim_gpu_cknn_cwire1_re2 <<< grid, threads >>>(cz_d, x_d, y_d, z_d, rad_d, ipc_d, ips_d, ise_d, nsmax, npls, ka, eta0, npls, 0, 0);   // device 2 
	// ----- PEC -----
	if (ienv == 2){
		wim_gpu_cknn_cwire1_im1 <<< grid, threads >>>(cz_d, x_d, y_d, z_d, rad_d, ipc_d, ips_d, ise_d, nsmax, npls, ka, eta0, sgnx, sgny, sgnz, sgnenv, npls, 0, 0);
		wim_gpu_cknn_cwire1_im2 <<< grid, threads >>>(cz_d, x_d, y_d, z_d, rad_d, ipc_d, ips_d, ise_d, nsmax, npls, ka, eta0, sgnx, sgny, sgnz, sgnenv, npls, 0, 0);
	}
	// ----- PEC -----
	wim_gpu_cknn_cwire2_re1 <<< grid, threads >>>(cz_d, x_d, y_d, z_d, rad_d, ipc_d, ips_d, ise_d, nsmax, npls, ka, eta0, npls, 0, 0);
	wim_gpu_cknn_cwire2_re2 <<< grid, threads >>>(cz_d, x_d, y_d, z_d, rad_d, ipc_d, ips_d, ise_d, nsmax, npls, ka, eta0, npls, 0, 0);
	if (ienv == 2){
		wim_gpu_cknn_cwire2_im1 <<< grid, threads >>>(cz_d, x_d, y_d, z_d, rad_d, ipc_d, ips_d, ise_d, nsmax, npls, ka, eta0, sgnx, sgny, sgnz, sgnenv, npls, 0, 0);
		wim_gpu_cknn_cwire2_im2 <<< grid, threads >>>(cz_d, x_d, y_d, z_d, rad_d, ipc_d, ips_d, ise_d, nsmax, npls, ka, eta0, sgnx, sgny, sgnz, sgnenv, npls, 0, 0);
	}
	// ----- PEC -----
	wim_gpu_ckmn_re1 <<<grid, threads >>>(cz_d, x_d, y_d, z_d, rad_d, ipc_d, ips_d, ise_d, nsmax, npls, ka, eta0, npls, 0, 0);
	wim_gpu_ckmn_re2 <<<grid, threads >>>(cz_d, x_d, y_d, z_d, rad_d, ipc_d, ips_d, ise_d, nsmax, npls, ka, eta0, npls, 0, 0);
	// ----- PEC -----
	if (ienv == 2){
		wim_gpu_ckmn_im1 <<<grid, threads >>>(cz_d, x_d, y_d, z_d, rad_d, ipc_d, ips_d, ise_d, nsmax, npls, ka, eta0, sgnx, sgny, sgnz, sgnenv, npls, 0, 0);
		wim_gpu_ckmn_im2 <<<grid, threads >>>(cz_d, x_d, y_d, z_d, rad_d, ipc_d, ips_d, ise_d, nsmax, npls, ka, eta0, sgnx, sgny, sgnz, sgnenv, npls, 0, 0);
	}
	// ----- PEC -----  


	// Stop event (kernel execution only)
	cudaEventRecord(stop_kernel, 0);
	cudaEventSynchronize(stop_kernel);
	cudaEventElapsedTime(&elapsedTime_kernel, start_kernel, stop_kernel);
	std::cout << std::endl << " Process - " << rank << " zmatrix mix GPU execution time = " << std::setprecision(4) << elapsedTime_kernel / 1.0e3 << " s" << std::endl;

	// Copy impedance matrix from GPU to CPU
	cudaMemcpy(cz_h, cz_d, sizeof(cuDoubleComplex) * npls*npls, cudaMemcpyDeviceToHost);

	// Stop event (data transfer and kernel execution)
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&elapsedTime, start, stop);

	std::cout << std::endl << " Process - " << rank << " zmatrix mix GPU execution time including data transfer = " << std::setprecision(4) << elapsedTime / 1.0e3 << " s" << std::endl;
	std::cout << std::endl << " Process - " << rank << " zmatrix mix GPU Data transfer time = " << std::setprecision(4) << (elapsedTime - elapsedTime_kernel) / 1.0e3 << " s" << std::endl;

	fprintf(f, "GPU execution time = %.3f ms\n", elapsedTime_kernel / 1.0e0);
	fprintf(f, "GPU execution time including data transfer = %.3f ms\n", elapsedTime / 1.0e0);
	fprintf(f, "Data transfer time = %.3f ms\n", (elapsedTime - elapsedTime_kernel) / 1.0e0);
	fclose(f);


	//cudaDeviceSynchronize();

	//Destroy events
	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	cudaEventDestroy(start_kernel);
	cudaEventDestroy(stop_kernel);


	cudaFree(cz_d);
	cudaFree(x_d);
	cudaFree(y_d);
	cudaFree(z_d);
	cudaFree(rad_d);
	cudaFree(ipc_d);
	cudaFree(ips_d);
	cudaFree(ise_d);


} 
