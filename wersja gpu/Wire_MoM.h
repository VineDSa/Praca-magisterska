#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include "Configuration.h"
#include "Zmatrix.h"
#include "Complex.h"

#include <cublasXt.h>
#include <mpi.h>
#include <chrono>


// Thread block size
#define BLOCK_SIZE 8
#define BROADCASTER_RANK 0

#ifndef min
#define min(a,b)  (((a)<(b))?(a):(b))
#endif

#ifndef max
#define max(a,b)  (((a)<(b))?(b):(a))
#endif

#ifndef roundup
#define roundup(a, b) (b <= 0) ? (a) : (((a) + (b)-1) & ~((b)-1))
#endif

#ifndef CUDA_CALL
#define CUDA_CALL(res, str) { if (res != cudaSuccess) {std::cout << "CUDA Error : " << str " : " << __FILE__ <<" " << __LINE__ <<" : ERR " << cudaGetErrorName(res) << std::endl;} }
#endif

#ifndef CUBLAS_CALL
#define CUBLAS_CALL(res, str) { if (res != CUBLAS_STATUS_SUCCESS) {std::cout << "cuBLAS Error : " << str " : " << __FILE__ <<" " << __LINE__ <<" : ERR " << int(res) << std::endl;} }
#endif

void LU_cuBLAS(int n, int nrhs, cuDoubleComplex* ci_cu, cuDoubleComplex* cu_cu, cuDoubleComplex* cz, int batchSize, int rank);

void zmatrix_gpu(cuDoubleComplex* cz_h, float* x_h, float* y_h, float* z_h, float* rad_h,
	int* ipc_h, int* ips_h, int* ise_h, int nsmax, double ka, int npls,
	int ienv, float sgnx, float sgny, float sgnz, float sgnenv, int rank);

void CalculateByFreq(int i, int npls, int nfreq, float freq,
	float* x, float* y, float* z, float* rad, int* ipc, int* ips, int* ise, const int nsmax,
	int ienv, float sgnx, float sgny, float sgnz, float sgnenv, int iload,
	float phas, int ngap, float ampl, cuDoubleComplex* cz, cuDoubleComplex* cu_cu, cuDoubleComplex* ci_cu, int rank);
