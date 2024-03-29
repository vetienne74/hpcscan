
//-------------------------------------------------------------------------------------------------------
// This grid is activated with command line option -testMode CUDA
// Derived class from Grid
// CUDA implementation (target GPU)
//-------------------------------------------------------------------------------------------------------

#include "grid_Cuda.h"

#include <algorithm> // for min and max
#include <cassert>
#include <cfloat>  // for FLT_MAX ;
#include <cmath>   // for fabs
#include <cstddef> // for NULL
#include <fstream>
#include <stdio.h>

#include "mpi.h"

#include "config.h"
#include "constant.h"
#include "fdm.h"
#include "global.h"
#include "output_report.h"

#include <thrust/device_vector.h>
#include <thrust/reduce.h>

#include <cooperative_groups.h>
// #include <cooperative_groups/reduce.h>

namespace cg = cooperative_groups;

using namespace std;

namespace hpcscan {

//-------------------------------------------------------------------------------------------------------

//Macro for checking cuda errors following a cuda launch or api call
#define cudaCheckError() {                                          \
	cudaError_t e=cudaGetLastError();                                 \
	if(e!=cudaSuccess) {                                              \
	  printf("Cuda failure %s:%d: '%s'\n",__FILE__,__LINE__,cudaGetErrorString(e));           \
	  exit(0); \
	}                                                                 \
   }

//-------------------------------------------------------------------------------------------------------

// there is probably an easier way to implement this (3d blocks?)
__global__ void kernel_fill_const(Myfloat *data, Myfloat val, const int n1, const int n2, const int n3, Myint64 i1Start, Myint64 i1End, Myint64 i2Start, Myint64 i2End, Myint64 i3Start, Myint64 i3End)
{
	int size = n1*n2*n3;
	int tid = threadIdx.x + blockIdx.x*blockDim.x;

	// 1D to 3D index (thanks StackOverflow)
	// public int[] to3D( int idx ) {
	// 	final int z = idx / (xMax * yMax);
	// 	idx -= (z * xMax * yMax);
	// 	final int y = idx / xMax;
	// 	final int x = idx % xMax;
	// 	return new int[]{ x, y, z };
	// }

	while (tid < size)
	{
		int i3 = tid / (n1*n2);
		int idx = tid-i3*n1*n2;
		int i2 = idx/n1;
		int i1 = idx%n1;

		if (i1 >= i1Start && i1 <= i1End &&
			i2 >= i2Start && i2 <= i2End &&
			i3 >= i3Start && i3 <= i3End   )
		{
			data[tid] = val;
		}

		tid += blockDim.x * gridDim.x;
	}
}

//-------------------------------------------------------------------------------------------------------

__global__ void kernel_fill_sine(Myfloat *data, Myfloat64 param1, Myfloat64 param2, Myfloat64 param3, Myfloat64 amp, int n1, int n2, int n3, Myint64 i1Start, Myint64 i1End, Myint64 i2Start, Myint64 i2End, Myint64 i3Start, Myint64 i3End, Myfloat Orig1, Myfloat Orig2, Myfloat Orig3, Myfloat64 d1, Myfloat64 d2, Myfloat64 d3 )
{
	// printf("sine %f %f %f %f %f\n",param1,param2,param3,amp,Orig1);

	int size = n1*n2*n3;
	int tid = threadIdx.x + blockIdx.x*blockDim.x;

	while (tid < size)
	{
		int i3 = tid / (n1*n2);
		int idx = tid-i3*n1*n2;
		int i2 = idx/n1;
		int i1 = idx%n1;

		if (i1 >= i1Start && i1 <= i1End &&
			i2 >= i2Start && i2 <= i2End &&
			i3 >= i3Start && i3 <= i3End   )
		{
			Myfloat64 coord1 = Myfloat64(Orig1 + i1 * d1);
			Myfloat64 coord2 = Myfloat64(Orig2 + i2 * d2);
			Myfloat64 coord3 = Myfloat64(Orig3 + i3 * d3);

			Myfloat val = amp * sin(coord1 * param1) * sin(coord2 * param2) * sin(coord3 * param3);

			data[tid] = val;
			//printf("data[%d]=%f\n",tid,val);
		}

		tid += blockDim.x * gridDim.x;
	}
}

//-------------------------------------------------------------------------------------------------------

__global__ void kernel_fill_linear(Myfloat *data, Myfloat64 param1, Myfloat64 param2, Myfloat64 param3, Myfloat64 amp, int n1, int n2, int n3, Myint64 i1Start, Myint64 i1End, Myint64 i2Start, Myint64 i2End, Myint64 i3Start, Myint64 i3End, Myfloat Orig1, Myfloat Orig2, Myfloat Orig3, Myfloat64 d1, Myfloat64 d2, Myfloat64 d3)
{
	// printf("linear %f %f %f %f\n",param1,param2,param3,amp);

	int size = n1*n2*n3;
	int tid = threadIdx.x + blockIdx.x*blockDim.x;

	while (tid < size)
	{
		int i3 = tid / (n1*n2);
		int idx = tid-i3*n1*n2;
		int i2 = idx/n1;
		int i1 = idx%n1;

		if (i1 >= i1Start && i1 <= i1End &&
			i2 >= i2Start && i2 <= i2End &&
			i3 >= i3Start && i3 <= i3End   )
		{
			Myfloat64 coord1 = Myfloat64(Orig1 + i1 * d1);
			Myfloat64 coord2 = Myfloat64(Orig2 + i2 * d2);
			Myfloat64 coord3 = Myfloat64(Orig3 + i3 * d3);

			Myfloat val = amp * coord1 * coord2 * coord3;

			data[tid] = val;
			//printf("data[%d]=%f\n",tid,val);
		}

		tid += blockDim.x * gridDim.x;
	}
}

//-------------------------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------------------------

__global__ void kernel_diff(Myfloat *data1, Myfloat *data2, Myfloat *dataOut, int n1, int n2, int n3, Myint64 i1Start, Myint64 i1End, Myint64 i2Start, Myint64 i2End, Myint64 i3Start, Myint64 i3End)
{
	int size = n1*n2*n3;
	int tid = threadIdx.x + blockIdx.x*blockDim.x;

	cg::thread_block cta = cg::this_thread_block();
	/*extern*/ __shared__ float sdata[256];

	sdata[threadIdx.x]=0;

	while (tid < size)
	{
		int i3 = tid / (n1*n2);
		int idx = tid-i3*n1*n2;
		int i2 = idx/n1;
		int i1 = idx%n1;


		if (i1 >= i1Start && i1 <= i1End &&
			i2 >= i2Start && i2 <= i2End &&
			i3 >= i3Start && i3 <= i3End   )
		{
			sdata[threadIdx.x] += fabsf(data1[tid]-data2[tid]);
		}

		tid += blockDim.x * gridDim.x;
	}

	cg::sync(cta);
	for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1) 
	{
		if (threadIdx.x < s) 
		{
		  sdata[threadIdx.x] += sdata[threadIdx.x + s];
		}	
		cg::sync(cta);
	  }
	
	  // write result for this block to global mem
	  if (threadIdx.x == 0) dataOut[blockIdx.x] = sdata[0];
}

//-------------------------------------------------------------------------------------------------------

__global__ void kernel_fabsf(Myfloat *data, Myfloat *dataOut, int n1, int n2, int n3, Myint64 i1Start, Myint64 i1End, Myint64 i2Start, Myint64 i2End, Myint64 i3Start, Myint64 i3End)
{
	int size = n1*n2*n3;
	int tid = threadIdx.x + blockIdx.x*blockDim.x;

	while (tid < size)
	{
		int i3 = tid / (n1*n2);
		int idx = tid-i3*n1*n2;
		int i2 = idx/n1;
		int i1 = idx%n1;

		dataOut[tid]=0;

		if (i1 >= i1Start && i1 <= i1End &&
			i2 >= i2Start && i2 <= i2End &&
			i3 >= i3Start && i3 <= i3End   )
		{
			dataOut[tid] = fabsf(data[tid]);
		}

		tid += blockDim.x * gridDim.x;
	}
}

//-------------------------------------------------------------------------------------------------------

__global__ void kernel_min(Myfloat *data, Myfloat *dataOut,
		int n1, int n2, int n3, Myint64 i1Start, Myint64 i1End, Myint64 i2Start, Myint64 i2End, Myint64 i3Start, Myint64 i3End)
{
	int size = n1*n2*n3;
	int tid = threadIdx.x + blockIdx.x*blockDim.x;

	cg::thread_block cta = cg::this_thread_block();
	/*extern*/ __shared__ float sdata[256];
	sdata[threadIdx.x]=99999;

	while (tid < size)
	{
		int i3 = tid / (n1*n2);
		int idx = tid-i3*n1*n2;
		int i2 = idx/n1;
		int i1 = idx%n1;


		if (i1 >= i1Start && i1 <= i1End &&
			i2 >= i2Start && i2 <= i2End &&
			i3 >= i3Start && i3 <= i3End   )
		{
			Myfloat val = data[tid];
			if (val < sdata[threadIdx.x]) sdata[threadIdx.x] = val;
		}

		tid += blockDim.x * gridDim.x;
	}

	cg::sync(cta);
	for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1) 
	{
		if (threadIdx.x < s) 
		{
			Myfloat val = sdata[threadIdx.x + s];
			if (val < sdata[threadIdx.x]) sdata[threadIdx.x] = val;
		}	
		cg::sync(cta);
	  }
	
	  // write result for this block to global mem
	  if (threadIdx.x == 0) dataOut[blockIdx.x] = sdata[0];
}

//-------------------------------------------------------------------------------------------------------
__global__ void kernel_mask(Myfloat *data, Myfloat *dataOut, Myfloat val, int n1, int n2, int n3, Myint64 i1Start, Myint64 i1End, Myint64 i2Start, Myint64 i2End, Myint64 i3Start, Myint64 i3End)
{
	int size = n1*n2*n3;
	int tid = threadIdx.x + blockIdx.x*blockDim.x;

	while (tid < size)
	{
		int i3 = tid / (n1*n2);
		int idx = tid-i3*n1*n2;
		int i2 = idx/n1;
		int i1 = idx%n1;

		dataOut[tid]=val;

		if (i1 >= i1Start && i1 <= i1End &&
			i2 >= i2Start && i2 <= i2End &&
			i3 >= i3Start && i3 <= i3End   )
		{
			dataOut[tid] = data[tid];
		}

		tid += blockDim.x * gridDim.x;
	}
}

//-------------------------------------------------------------------------------------------------------

// #pragma omp parallel for collapse(2)
// for (Myint64 i3 = i3Start; i3<= i3End; i3++)
// 	for (Myint64 i2 = i2Start; i2<= i2End; i2++)
// 		for (Myint64 i1 = i1Start; i1<= i1End; i1++)
// 			prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
// 					coef[i1+i2*n1+i3*n1*n2] * lapla[i1+i2*n1+i3*n1*n2] ;

__global__ void kernel_updatePressure(Myfloat *prn, Myfloat *prc, Myfloat *coef, Myfloat *lapla, int n1, int n2, int n3, Myint64 i1Start, Myint64 i1End, Myint64 i2Start, Myint64 i2End, Myint64 i3Start, Myint64 i3End)
{
	Myint64 size = n1*n2*n3;
	Myint64 tid = threadIdx.x + blockIdx.x*blockDim.x;

	while (tid < size)
	{
		int i3 = tid / (n1*n2);
		int idx = tid-i3*n1*n2;
		int i2 = idx/n1;
		int i1 = idx%n1;

		if (i1 >= i1Start && i1 <= i1End &&
			i2 >= i2Start && i2 <= i2End &&
			i3 >= i3Start && i3 <= i3End)
		{
			prn[tid] = Myfloat(2.0) * prc[tid] - prn[tid] + coef[tid] * lapla[tid];
		}

		tid += blockDim.x * gridDim.x;
	}
}

//-------------------------------------------------------------------------------------------------------

__global__ void kernel_applyBoundaryCondition(Myfloat *data, int n1, int n2, int n3, 
			Myint64 i1halo1_i1Start, Myint64 i1halo1_i1End, Myint64 i1halo1_i2Start, Myint64 i1halo1_i2End, Myint64 i1halo1_i3Start, Myint64 i1halo1_i3End,
			Myint64 i1halo2_i1Start, Myint64 i1halo2_i1End, Myint64 i1halo2_i2Start, Myint64 i1halo2_i2End, Myint64 i1halo2_i3Start, Myint64 i1halo2_i3End,
			Myint64 i2halo1_i1Start, Myint64 i2halo1_i1End, Myint64 i2halo1_i2Start, Myint64 i2halo1_i2End, Myint64 i2halo1_i3Start, Myint64 i2halo1_i3End,
			Myint64 i2halo2_i1Start, Myint64 i2halo2_i1End, Myint64 i2halo2_i2Start, Myint64 i2halo2_i2End, Myint64 i2halo2_i3Start, Myint64 i2halo2_i3End,
			Myint64 i3halo1_i1Start, Myint64 i3halo1_i1End, Myint64 i3halo1_i2Start, Myint64 i3halo1_i2End, Myint64 i3halo1_i3Start, Myint64 i3halo1_i3End,
			Myint64 i3halo2_i1Start, Myint64 i3halo2_i1End, Myint64 i3halo2_i2Start, Myint64 i3halo2_i2End, Myint64 i3halo2_i3Start, Myint64 i3halo2_i3End)
			
{
	int size = n1*n2*n3;
	int tid = threadIdx.x + blockIdx.x*blockDim.x;

	while (tid < size)
	{
		int i3 = tid / (n1*n2);
		int idx = tid-i3*n1*n2;
		int i2 = idx/n1;
		int i1 = idx%n1;

		// I1HALO1
		Myint64 iInner = i1halo1_i1End+1;
		if (tid == iInner+i2*n1+i3*n1*n2) data[tid] = 0.0 ;

		if (i1 >= i1halo1_i1Start && i1 <= i1halo1_i1End &&
			i2 >= i1halo1_i2Start && i2 <= i1halo1_i2End &&
			i3 >= i1halo1_i3Start && i3 <= i1halo1_i3End   )
		{
			data[tid] = -data[(iInner+iInner-i1)+i2*n1+i3*n1*n2];
		}

		// I1HALO2
		iInner = i1halo2_i1Start-1;
		if (tid == iInner+i2*n1+i3*n1*n2) data[tid] = 0.0 ;

		if (i1 >= i1halo2_i1Start && i1 <= i1halo2_i1End &&
			i2 >= i1halo2_i2Start && i2 <= i1halo2_i2End &&
			i3 >= i1halo2_i3Start && i3 <= i1halo2_i3End   )
		{
			data[tid] = -data[(iInner-(i1-iInner))+i2*n1+i3*n1*n2];
		}

		// I2HALO1
		iInner = i2halo1_i2End+1;
		if (tid == i1+iInner*n1+i3*n1*n2) data[tid] = 0.0 ;

		if (i1 >= i2halo1_i1Start && i1 <= i2halo1_i1End &&
			i2 >= i2halo1_i2Start && i2 <= i2halo1_i2End &&
			i3 >= i2halo1_i3Start && i3 <= i2halo1_i3End   )
		{
			data[tid] = -data[i1+(iInner+iInner-i2)*n1+i3*n1*n2];
		}

		// I2HALO2
		iInner = i2halo2_i2Start-1;
		if (tid == i1+iInner*n1+i3*n1*n2) data[tid] = 0.0 ;

		if (i1 >= i2halo2_i1Start && i1 <= i2halo2_i1End &&
			i2 >= i2halo2_i2Start && i2 <= i2halo2_i2End &&
			i3 >= i2halo2_i3Start && i3 <= i2halo2_i3End   )
		{
			data[tid] = -data[i1+(iInner-(i2-iInner))*n1+i3*n1*n2];
		}

		// I3HALO1
		iInner = i3halo1_i3End+1;
		if (tid == i1+i2*n1+iInner*n1*n2) data[tid] = 0.0 ;

		if (i1 >= i3halo1_i1Start && i1 <= i3halo1_i1End &&
			i2 >= i3halo1_i2Start && i2 <= i3halo1_i2End &&
			i3 >= i3halo1_i3Start && i3 <= i3halo1_i3End   )
		{
			data[tid] = -data[i1+i2*n1+(iInner+iInner-i3)*n1*n2];			                  
		}

		// I3HALO2
		iInner = i3halo2_i3Start-1;
		if (tid == i1+i2*n1+iInner*n1*n2) data[tid] = 0.0 ;

		if (i1 >= i3halo2_i1Start && i1 <= i3halo2_i1End &&
			i2 >= i3halo2_i2Start && i2 <= i3halo2_i2End &&
			i3 >= i3halo2_i3Start && i3 <= i3halo2_i3End   )
		{
			data[tid] = -data[i1+i2*n1+(iInner-(i3-iInner))*n1*n2];
		}

		tid += blockDim.x * gridDim.x;
	}
}

//-------------------------------------------------------------------------------------------------------


// 	for (Myint64 i3 = i3Start; i3<= i3End; i3++)
// 	  for (Myint64 i2 = i2Start; i2<= i2End; i2++)
// 	    for (Myint64 i1 = i1Start; i1<= i1End; i1++)
// 		  prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] + coef[i1+i2*n1+i3*n1*n2] *
// 			(FD_D2_O4_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
// 		   + FD_D2_O4_N2(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
// 		   + FD_D2_O4_N3(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)) ;

// #define FD_D2_O4_N1(U, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) \
// 		(        (FD_D2_O4_A0 *  U[i1   + i2*n1 + i3*n2*n1] \
// 				+ FD_D2_O4_A1 * (U[i1+1 + i2*n1 + i3*n2*n1] + U[i1-1 + i2*n1 + i3*n2*n1])  \
// 				+ FD_D2_O4_A2 * (U[i1+2 + i2*n1 + i3*n2*n1] + U[i1-2 + i2*n1 + i3*n2*n1])) \
// 				* inv2_d1)

// #define FD_D2_O4_N2(U, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) \
// 		(        (FD_D2_O4_A0 *  U[i1 + i2    *n1 + i3*n2*n1] \
// 				+ FD_D2_O4_A1 * (U[i1 + (i2+1)*n1 + i3*n2*n1] + U[i1 + (i2-1)*n1 + i3*n2*n1])  \
// 				+ FD_D2_O4_A2 * (U[i1 + (i2+2)*n1 + i3*n2*n1] + U[i1 + (i2-2)*n1 + i3*n2*n1])) \
// 				* inv2_d2)

// #define FD_D2_O4_N3(U, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) \
// 		(        (FD_D2_O4_A0 *  U[i1 + i2*n1 + i3    *n2*n1] \
// 				+ FD_D2_O4_A1 * (U[i1 + i2*n1 + (i3+1)*n2*n1] + U[i1 + i2*n1 + (i3-1)*n2*n1])  \
// 				+ FD_D2_O4_A2 * (U[i1 + i2*n1 + (i3+2)*n2*n1] + U[i1 + i2*n1 + (i3-2)*n2*n1])) \
// 				* inv2_d3)
__global__ void kernel_computePressureWithFD_O4(Myfloat *prn, Myfloat *prc, Myfloat *coef, Myfloat inv2_d1, Myfloat inv2_d2, Myfloat inv2_d3, int n1, int n2, int n3, Myint64 i1Start, Myint64 i1End, Myint64 i2Start, Myint64 i2End, Myint64 i3Start, Myint64 i3End)
{
	int size = n1*n2*n3;
	int tid = threadIdx.x + blockIdx.x*blockDim.x;

	const Myfloat FD_D2_O4_A0   =  -5./2. ;
	const Myfloat FD_D2_O4_A1   =  4./3. ;
	const Myfloat FD_D2_O4_A2   =  -1./12. ;
	

	while (tid < size)
	{
		int i3 = tid / (n1*n2);
		int idx = tid-i3*n1*n2;
		int i2 = idx/n1;
		int i1 = idx%n1;

		if (i1 >= i1Start && i1 <= i1End &&
			i2 >= i2Start && i2 <= i2End &&
			i3 >= i3Start && i3 <= i3End   )
		{
			prn[i1+i2*n1+i3*n1*n2] = 2.0 * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] + coef[i1+i2*n1+i3*n1*n2] *
			(((FD_D2_O4_A0 *  prc[i1   + i2*n1 + i3*n2*n1]
			 + FD_D2_O4_A1 * (prc[i1+1 + i2*n1 + i3*n2*n1] + prc[i1-1 + i2*n1 + i3*n2*n1]) 
			 + FD_D2_O4_A2 * (prc[i1+2 + i2*n1 + i3*n2*n1] + prc[i1-2 + i2*n1 + i3*n2*n1]))
			 * inv2_d1)
			+((FD_D2_O4_A0 *  prc[i1 +  i2   *n1 + i3*n2*n1] 
			 + FD_D2_O4_A1 * (prc[i1 + (i2+1)*n1 + i3*n2*n1] + prc[i1 + (i2-1)*n1 + i3*n2*n1])  
			 + FD_D2_O4_A2 * (prc[i1 + (i2+2)*n1 + i3*n2*n1] + prc[i1 + (i2-2)*n1 + i3*n2*n1])) 
			 * inv2_d2)
			+((FD_D2_O4_A0 *  prc[i1 + i2*n1 +  i3   *n2*n1] 
			 + FD_D2_O4_A1 * (prc[i1 + i2*n1 + (i3+1)*n2*n1] + prc[i1 + i2*n1 + (i3-1)*n2*n1])  
			 + FD_D2_O4_A2 * (prc[i1 + i2*n1 + (i3+2)*n2*n1] + prc[i1 + i2*n1 + (i3-2)*n2*n1])) 
			 * inv2_d3));
		}

		tid += blockDim.x * gridDim.x;
	}
}

//-------------------------------------------------------------------------------------------------------

// #define FD_D2_O8_N1(U, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) \
// 		(        (FD_D2_O8_A0 *  U[i1   + i2*n1 + i3*n2*n1] \
// 				+ FD_D2_O8_A1 * (U[i1+1 + i2*n1 + i3*n2*n1] + U[i1-1 + i2*n1 + i3*n2*n1])  \
// 				+ FD_D2_O8_A2 * (U[i1+2 + i2*n1 + i3*n2*n1] + U[i1-2 + i2*n1 + i3*n2*n1])  \
// 				+ FD_D2_O8_A3 * (U[i1+3 + i2*n1 + i3*n2*n1] + U[i1-3 + i2*n1 + i3*n2*n1])  \
// 				+ FD_D2_O8_A4 * (U[i1+4 + i2*n1 + i3*n2*n1] + U[i1-4 + i2*n1 + i3*n2*n1])) \
// 				* inv2_d1)

// #define FD_D2_O8_N2(U, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) \
// 		(        (FD_D2_O8_A0 *  U[i1 + i2    *n1 + i3*n2*n1] \
// 				+ FD_D2_O8_A1 * (U[i1 + (i2+1)*n1 + i3*n2*n1] + U[i1 + (i2-1)*n1 + i3*n2*n1])  \
// 				+ FD_D2_O8_A2 * (U[i1 + (i2+2)*n1 + i3*n2*n1] + U[i1 + (i2-2)*n1 + i3*n2*n1])  \
// 				+ FD_D2_O8_A3 * (U[i1 + (i2+3)*n1 + i3*n2*n1] + U[i1 + (i2-3)*n1 + i3*n2*n1])  \
// 				+ FD_D2_O8_A4 * (U[i1 + (i2+4)*n1 + i3*n2*n1] + U[i1 + (i2-4)*n1 + i3*n2*n1])) \
// 				* inv2_d2)

// #define FD_D2_O8_N3(U, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) \
// 		(        (FD_D2_O8_A0 *  U[i1 + i2*n1 + i3    *n2*n1] \
// 				+ FD_D2_O8_A1 * (U[i1 + i2*n1 + (i3+1)*n2*n1] + U[i1 + i2*n1 + (i3-1)*n2*n1])  \
// 				+ FD_D2_O8_A2 * (U[i1 + i2*n1 + (i3+2)*n2*n1] + U[i1 + i2*n1 + (i3-2)*n2*n1])  \
// 				+ FD_D2_O8_A3 * (U[i1 + i2*n1 + (i3+3)*n2*n1] + U[i1 + i2*n1 + (i3-3)*n2*n1])  \
// 				+ FD_D2_O8_A4 * (U[i1 + i2*n1 + (i3+4)*n2*n1] + U[i1 + i2*n1 + (i3-4)*n2*n1])) \
// 				* inv2_d3)
__global__ void kernel_computePressureWithFD_O8(Myfloat *prn, Myfloat *prc, Myfloat *coef, Myfloat inv2_d1, Myfloat inv2_d2, Myfloat inv2_d3, int n1, int n2, int n3, Myint64 i1Start, Myint64 i1End, Myint64 i2Start, Myint64 i2End, Myint64 i3Start, Myint64 i3End)
{
	int size = n1*n2*n3;
	int tid = threadIdx.x + blockIdx.x*blockDim.x;

	const Myfloat FD_D2_O8_A0   = -205./72. ;
	const Myfloat FD_D2_O8_A1   = 8./5. ;
	const Myfloat FD_D2_O8_A2   = -1/5. ;
	const Myfloat FD_D2_O8_A3   = 8./315. ;
	const Myfloat FD_D2_O8_A4   = -1/560. ;

	while (tid < size)
	{
		int i3 = tid / (n1*n2);
		int idx = tid-i3*n1*n2;
		int i2 = idx/n1;
		int i1 = idx%n1;

		if (i1 >= i1Start && i1 <= i1End &&
			i2 >= i2Start && i2 <= i2End &&
			i3 >= i3Start && i3 <= i3End   )
		{
			prn[i1+i2*n1+i3*n1*n2] = 2.0 * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] + coef[i1+i2*n1+i3*n1*n2] *
				((   (FD_D2_O8_A0 *  prc[i1   + i2*n1 + i3*n2*n1] 
					+ FD_D2_O8_A1 * (prc[i1+1 + i2*n1 + i3*n2*n1] + prc[i1-1 + i2*n1 + i3*n2*n1])  
					+ FD_D2_O8_A2 * (prc[i1+2 + i2*n1 + i3*n2*n1] + prc[i1-2 + i2*n1 + i3*n2*n1])  
					+ FD_D2_O8_A3 * (prc[i1+3 + i2*n1 + i3*n2*n1] + prc[i1-3 + i2*n1 + i3*n2*n1])  
					+ FD_D2_O8_A4 * (prc[i1+4 + i2*n1 + i3*n2*n1] + prc[i1-4 + i2*n1 + i3*n2*n1])) 
					* inv2_d1)
			+(       (FD_D2_O8_A0 *  prc[i1 + i2    *n1 + i3*n2*n1] 
					+ FD_D2_O8_A1 * (prc[i1 + (i2+1)*n1 + i3*n2*n1] + prc[i1 + (i2-1)*n1 + i3*n2*n1])  
					+ FD_D2_O8_A2 * (prc[i1 + (i2+2)*n1 + i3*n2*n1] + prc[i1 + (i2-2)*n1 + i3*n2*n1])  
					+ FD_D2_O8_A3 * (prc[i1 + (i2+3)*n1 + i3*n2*n1] + prc[i1 + (i2-3)*n1 + i3*n2*n1])  
					+ FD_D2_O8_A4 * (prc[i1 + (i2+4)*n1 + i3*n2*n1] + prc[i1 + (i2-4)*n1 + i3*n2*n1])) 
					* inv2_d2)
			+(       (FD_D2_O8_A0 *  prc[i1 + i2*n1 + i3    *n2*n1] 
					+ FD_D2_O8_A1 * (prc[i1 + i2*n1 + (i3+1)*n2*n1] + prc[i1 + i2*n1 + (i3-1)*n2*n1])  
					+ FD_D2_O8_A2 * (prc[i1 + i2*n1 + (i3+2)*n2*n1] + prc[i1 + i2*n1 + (i3-2)*n2*n1])  
					+ FD_D2_O8_A3 * (prc[i1 + i2*n1 + (i3+3)*n2*n1] + prc[i1 + i2*n1 + (i3-3)*n2*n1])  
					+ FD_D2_O8_A4 * (prc[i1 + i2*n1 + (i3+4)*n2*n1] + prc[i1 + i2*n1 + (i3-4)*n2*n1])) 
					* inv2_d3));
		}

		tid += blockDim.x * gridDim.x;
	}
}

//-------------------------------------------------------------------------------------------------------

Grid_Cuda::Grid_Cuda(Grid_type gridTypeIn) : Grid(gridTypeIn)
														{
	printDebug(MID_DEBUG, "IN Grid_Cuda::Grid_Cuda");

	gridMode = "Cuda" ;
	
	d_grid_3d = NULL;
	d_help_3d = NULL;

	printDebug(MID_DEBUG, "OUT Grid_Cuda::Grid_Cuda");
														}

//-------------------------------------------------------------------------------------------------------

Grid_Cuda::Grid_Cuda(Grid_type gridTypeIn, Dim_type dimIn,
		Myint64 n1InnerIn, Myint64 n2InnerIn, Myint64 n3InnerIn) : Grid(gridTypeIn, dimIn,
				n1InnerIn, n2InnerIn, n3InnerIn)
{
	printDebug(MID_DEBUG, "IN Grid_Cuda::Grid_Cuda");

	gridMode = "Cuda" ;
	
	d_grid_3d = NULL;
	d_help_3d = NULL;

	printDebug(MID_DEBUG, "OUT Grid_Cuda::Grid_Cuda");
}

//-------------------------------------------------------------------------------------------------------

Grid_Cuda::~Grid_Cuda(void)
{
	printDebug(MID_DEBUG, "IN Grid_Cuda::~Grid_Cuda");

	//delete[] grid_3d ;
	cudaFree(d_grid_3d);
	cudaFree(d_help_3d);
	cudaCheckError();

	printDebug(MID_DEBUG, "OUT Grid_Cuda::~Grid_Cuda");
}

//-------------------------------------------------------------------------------------------------------

void Grid_Cuda::info(void)
{
	printDebug(FULL_DEBUG, "IN Grid_Cuda::info");

	// parent class info
	Grid::info() ;

	// additional info
	// TO DO

	printDebug(FULL_DEBUG, "IN Grid_Cuda::info");
}


//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_Cuda::FD_LAPLACIAN(Point_type pointType, const Grid& Wgrid, Myint fdOrder)
{
	printDebug(MID_DEBUG, "IN Grid_Cuda::FD_LAPLACIAN");

	// TO DO
	Grid::FD_LAPLACIAN(pointType, Wgrid, fdOrder) ;

	printDebug(MID_DEBUG, "OUT Grid_Cuda::FD_LAPLACIAN");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_Cuda::computePressureWithFD(Grid& prcGridIn, Grid& coefGridIn, Myint fdOrder)
{

	printDebug(FULL_DEBUG, "In Grid_Cuda::computePressureWithFD") ;

	// TO DO
	// Grid::computePressureWithFD(prcGridIn, coefGridIn, fdOrder) ;

	//pointType
	Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
	getGridIndex(INNER_POINTS, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End);

	if (fdOrder != 4 && fdOrder != 8)
	{
		printError("CUDA: only FD order 4/8 implemented");
	}

	const Myfloat inv_d1  = Myfloat(1.0) / d1 ;
	const Myfloat inv_d2  = Myfloat(1.0) / d2 ;
	const Myfloat inv_d3  = Myfloat(1.0) / d3 ;

	const Myfloat inv2_d1 = inv_d1 * inv_d1 ;
	const Myfloat inv2_d2 = inv_d2 * inv_d2 ;
	const Myfloat inv2_d3 = inv_d3 * inv_d3 ;

	Myfloat *prc_d_grid_3d = ((Grid_Cuda&) prcGridIn).d_grid_3d ;
	Myfloat *coef_d_grid_3d = ((Grid_Cuda&) coefGridIn).d_grid_3d ;

	if (fdOrder == 4)
	{
		kernel_computePressureWithFD_O4<<<1024,256>>>(d_grid_3d, prc_d_grid_3d, coef_d_grid_3d,inv2_d1,inv2_d2,inv2_d3,n1,n2,n3,i1Start,i1End,i2Start,i2End,i3Start,i3End);
	}
	else if (fdOrder == 8)
	{
		kernel_computePressureWithFD_O8<<<1024,256>>>(d_grid_3d, prc_d_grid_3d, coef_d_grid_3d,inv2_d1,inv2_d2,inv2_d3,n1,n2,n3,i1Start,i1End,i2Start,i2End,i3Start,i3End);
	}
	
	cudaCheckError();
	cudaDeviceSynchronize();
	
	printDebug(FULL_DEBUG, "Out Grid_Cuda::computePressureWithFD") ;
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_Cuda::initializeGrid(void)
{
	printDebug(FULL_DEBUG, "In Grid_Cuda::initializeGrid") ;

	Grid::initializeGrid() ; // this sets up halos etc.

	if (d_grid_3d == NULL)
	{
		cudaMalloc( (void**)&d_grid_3d, n1*n2*n3*sizeof(Myfloat) );
		cudaCheckError();

		cudaMalloc( (void**)&d_help_3d, 1024*sizeof(Myfloat) );
		cudaCheckError();
	}
	printDebug(FULL_DEBUG, "Out Grid_Cuda::initializeGrid") ;
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------
void Grid_Cuda::fill(Point_type pointType, Myfloat val)
{
	printDebug(FULL_DEBUG, "In Grid_Cuda::fill") ;

	//pointType
	Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
	Grid::getGridIndex(pointType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End);

	kernel_fill_const<<<2048,128>>>(d_grid_3d,val,n1,n2,n3,i1Start,i1End,i2Start,i2End,i3Start,i3End);
	cudaDeviceSynchronize();
	cudaCheckError();

	printDebug(FULL_DEBUG, "Out Grid_Cuda::fill") ;
}

//-------------------------------------------------------------------------------------------------------
void Grid_Cuda::fill(Point_type pointType, Func_type t1,  Func_type t2, Func_type t3,
		Myfloat64 param1, Myfloat64 param2, Myfloat64 param3, Myfloat64 amp)
{
	printDebug(FULL_DEBUG, "In Grid_Cuda::fill") ;

	Grid::fill(pointType, t1,  t2, t3, param1, param2, param3, amp) ;


	//pointType
	Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
	Grid::getGridIndex(pointType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End);

	// we only use sine and linear for now
	int ok = 1;
	if (t1 != FUNC_SINE && t1 != FUNC_LINEAR) ok = 0;
	if (t2 != FUNC_SINE && t2 != FUNC_LINEAR) ok = 0;
	if (t3 != FUNC_SINE && t3 != FUNC_LINEAR) ok = 0;
	if (!ok) printError("CUDA: only FUNC_SINE and FUNC_LINEAR implemented");

	if ((t1==t2 && t2==t3)==false) printError("CUDA: func has to be same in each dimension");

	if (t1 == FUNC_SINE) kernel_fill_sine<<<1024,128>>>  (d_grid_3d,param1,param2,param3,amp,n1,n2,n3,i1Start,i1End,i2Start,i2End,i3Start,i3End,Orig1,Orig2,Orig3,d1,d2,d2);
	else                 kernel_fill_linear<<<1024,128>>>(d_grid_3d,param1,param2,param3,amp,n1,n2,n3,i1Start,i1End,i2Start,i2End,i3Start,i3End,Orig1,Orig2,Orig3,d1,d2,d2);


	printDebug(FULL_DEBUG, "Out Grid_Cuda::fill") ;
}

//-------------------------------------------------------------------------------------------------------

__global__ void kernel_fillArray(Myfloat *gridOut, Myfloat val, Myint64 gridSize)
{
	Myint64 tid = threadIdx.x + blockIdx.x*blockDim.x;
	while (tid < gridSize)
	{
		gridOut[tid] = val ;
		tid += blockDim.x * gridDim.x;
	}
}

void Grid_Cuda::fillArray(Myfloat val)
{
	printDebug(MID_DEBUG, "IN Grid_Cuda::fillArray");

	Myint64 gridSize = n1*n2*n3;
	kernel_fillArray<<<1024,256>>>(d_grid_3d, val, gridSize) ;

	cudaDeviceSynchronize();

	printDebug(MID_DEBUG, "OUT Grid_Cuda::fillArray");
}

//-------------------------------------------------------------------------------------------------------

__global__ void kernel_copyArray(Myfloat *gridOut, Myfloat *gridIn, Myint64 gridSize)
{
	Myint64 tid = threadIdx.x + blockIdx.x*blockDim.x;
	while (tid < gridSize)
	{
		gridOut[tid] = gridIn[tid] ;
		tid += blockDim.x * gridDim.x;
	}
}

void Grid_Cuda::copyArray(const Grid& gridIn)
{
	printDebug(MID_DEBUG, "IN Grid_Cuda::copyArray");

	Myfloat *gridIn_d_grid_3d = ((Grid_Cuda&) gridIn).d_grid_3d ;
	Myint64 gridSize = n1*n2*n3;
	kernel_copyArray<<<1024,256>>>(d_grid_3d, gridIn_d_grid_3d, gridSize) ;

	cudaDeviceSynchronize();

	printDebug(MID_DEBUG, "OUT Grid_Cuda::copyArray");
}

//-------------------------------------------------------------------------------------------------------

__global__ void kernel_addArray(Myfloat *gridOut, Myfloat *gridIn1, Myfloat *gridIn2, Myint64 gridSize)
{
	Myint64 tid = threadIdx.x + blockIdx.x*blockDim.x;
	while (tid < gridSize)
	{
		gridOut[tid] = gridIn1[tid] + gridIn2[tid] ;
		tid += blockDim.x * gridDim.x;
	}
}

void Grid_Cuda::addArray(const Grid& gridIn1, const Grid& gridIn2)
{
	printDebug(MID_DEBUG, "IN Grid_Cuda::addArray");

	Myfloat *gridIn1_d_grid_3d = ((Grid_Cuda&) gridIn1).d_grid_3d ;
	Myfloat *gridIn2_d_grid_3d = ((Grid_Cuda&) gridIn2).d_grid_3d ;
	Myint64 gridSize = n1*n2*n3;
	kernel_addArray<<<1024,256>>>(d_grid_3d, gridIn1_d_grid_3d, gridIn2_d_grid_3d, gridSize) ;

	cudaDeviceSynchronize();

	printDebug(MID_DEBUG, "OUT Grid_Cuda::addArray");
}

//-------------------------------------------------------------------------------------------------------

__global__ void kernel_multiplyArray(Myfloat *gridOut, Myfloat *gridIn1, Myfloat *gridIn2, Myint64 gridSize)
{
	Myint64 tid = threadIdx.x + blockIdx.x*blockDim.x;
	while (tid < gridSize)
	{
		gridOut[tid] = gridIn1[tid] * gridIn2[tid] ;
		tid += blockDim.x * gridDim.x;
	}
}

void Grid_Cuda::multiplyArray(const Grid& gridIn1, const Grid& gridIn2)
{
	printDebug(MID_DEBUG, "IN Grid_Cuda::multiplyArray");

	Myfloat *gridIn1_d_grid_3d = ((Grid_Cuda&) gridIn1).d_grid_3d ;
	Myfloat *gridIn2_d_grid_3d = ((Grid_Cuda&) gridIn2).d_grid_3d ;
	Myint64 gridSize = n1*n2*n3;
	kernel_addArray<<<1024,256>>>(d_grid_3d, gridIn1_d_grid_3d, gridIn2_d_grid_3d, gridSize) ;

	cudaDeviceSynchronize();

	printDebug(MID_DEBUG, "OUT Grid_Cuda::multiplyArray");
}

//-------------------------------------------------------------------------------------------------------

__global__ void kernel_addUpdateArray(Myfloat *gridOut, Myfloat *gridIn, Myint64 gridSize)
{
	Myint64 tid = threadIdx.x + blockIdx.x*blockDim.x;
	while (tid < gridSize)
	{
		gridOut[tid] = gridOut[tid] + gridIn[tid] ;
		tid += blockDim.x * gridDim.x;
	}
}

void Grid_Cuda::addUpdateArray(const Grid& gridIn)
{
	printDebug(MID_DEBUG, "IN Grid_Cuda::addUpdateArray");

	Myfloat *gridIn_d_grid_3d = ((Grid_Cuda&) gridIn).d_grid_3d ;
	Myint64 gridSize = n1*n2*n3;
	kernel_addUpdateArray<<<1024,256>>>(d_grid_3d, gridIn_d_grid_3d, gridSize) ;

	cudaDeviceSynchronize();

	printDebug(MID_DEBUG, "OUT Grid_Cuda::addUpdateArray");
}

//-------------------------------------------------------------------------------------------------------
Myfloat Grid_Cuda::getMin(Point_type pointType)
{
	printDebug(FULL_DEBUG, "In Grid_Cuda::getMin") ;

	if (true)
	{
		//pointType
		Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
		Grid::getGridIndex(pointType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End);
		
		const int numBlocks = 1024;
		kernel_min<<<numBlocks,256>>>(d_grid_3d,d_help_3d,n1,n2,n3,i1Start,i1End,i2Start,i2End,i3Start,i3End);
		cudaCheckError();

		thrust::device_ptr<Myfloat> d_help_3d_ptr = thrust::device_pointer_cast(d_help_3d);
		thrust::device_ptr<hpcscan::Myfloat> vptr = thrust::min_element(thrust::device, d_help_3d_ptr, d_help_3d_ptr + numBlocks);
		float val = *vptr;
		return val;
	}
	else if (true)
	{
		//pointType
		Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
		Grid::getGridIndex(pointType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End);
		kernel_mask<<<1024,256>>>(d_grid_3d,d_help_3d,999,n1,n2,n3,i1Start,i1End,i2Start,i2End,i3Start,i3End);
		thrust::device_ptr<Myfloat> d_help_3d_ptr = thrust::device_pointer_cast(d_help_3d);
		thrust::device_ptr<hpcscan::Myfloat> vptr = thrust::min_element(thrust::device, d_help_3d_ptr, d_help_3d_ptr + n1*n2*n3);
		float val = *vptr;
		return val;
	}
	else
	{
		thrust::device_ptr<Myfloat> d_help_3d_ptr = thrust::device_pointer_cast(d_grid_3d);
		thrust::device_ptr<hpcscan::Myfloat> vptr = thrust::min_element(thrust::device, d_help_3d_ptr, d_help_3d_ptr + n1*n2*n3);
		float val = *vptr;
		return val;
	}
	cudaDeviceSynchronize();
	printDebug(FULL_DEBUG, "Out Grid_Cuda::getMin") ;
}

//-------------------------------------------------------------------------------------------------------
Myfloat Grid_Cuda::getMax(Point_type pointType)
{
	printDebug(FULL_DEBUG, "In Grid_Cuda::getMax") ;

	if (hpcscan::nMpiProc > 1)
	{
		//pointType
		Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
		Grid::getGridIndex(pointType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End);
		kernel_mask<<<1024,256>>>(d_grid_3d,d_help_3d,0,n1,n2,n3,i1Start,i1End,i2Start,i2End,i3Start,i3End);
		thrust::device_ptr<Myfloat> d_help_3d_ptr = thrust::device_pointer_cast(d_help_3d);
		thrust::device_ptr<hpcscan::Myfloat> vptr = thrust::max_element(thrust::device, d_help_3d_ptr, d_help_3d_ptr + n1*n2*n3);
		float val = *vptr;
		return val;
	}
	else
	{
		thrust::device_ptr<Myfloat> d_help_3d_ptr = thrust::device_pointer_cast(d_grid_3d);
		thrust::device_ptr<hpcscan::Myfloat> vptr = thrust::max_element(thrust::device, d_help_3d_ptr, d_help_3d_ptr + n1*n2*n3);
		float val = *vptr;
		return val;
	}
	cudaDeviceSynchronize();
	printDebug(FULL_DEBUG, "Out Grid_Cuda::getMax") ;
}

//-------------------------------------------------------------------------------------------------------
Myfloat Grid_Cuda::L1Err(Point_type pointType, const Grid& gridIn) const
{
	printDebug(FULL_DEBUG, "In Grid_Cuda::L1Err") ;

	// TO DO
	// return(Grid::L1Err(pointType, gridIn)) ;

	//pointType
	Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
	Grid::getGridIndex(pointType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End);

	thrust::device_ptr<Myfloat> d_help_3d_ptr;

	const int numBlocks = 1024;
	
	Myfloat *gridIn_d_grid_3d = ((Grid_Cuda&) gridIn).d_grid_3d ;

	kernel_diff<<<numBlocks,256>>>(d_grid_3d, gridIn_d_grid_3d,d_help_3d,n1,n2,n3,i1Start,i1End,i2Start,i2End,i3Start,i3End);
	cudaCheckError();
	d_help_3d_ptr = thrust::device_pointer_cast(d_help_3d);
	double totErr = thrust::reduce(thrust::device, d_help_3d_ptr, d_help_3d_ptr + numBlocks);

	double totArr;
	if (false)
	{
		kernel_fabsf<<<1024,256>>>(gridIn_d_grid_3d, d_help_3d,n1,n2,n3,i1Start,i1End,i2Start,i2End,i3Start,i3End);
		totArr = thrust::reduce(thrust::device, d_help_3d_ptr, d_help_3d_ptr + n1*n2*n3);
	}
	else // assuming grid values are positive
	{
		d_help_3d_ptr = thrust::device_pointer_cast(gridIn_d_grid_3d);
		totArr = thrust::reduce(thrust::device, d_help_3d_ptr, d_help_3d_ptr + n1*n2*n3);
	}

	cudaDeviceSynchronize();

	if (totArr < MAX_ERR_FLOAT) totArr = 1.0 * npoint ;

	return totErr/totArr;

	

	printDebug(FULL_DEBUG, "Out Grid_Cuda::L1Err") ;
}
//-------------------------------------------------------------------------------------------------------
Rtn_code Grid_Cuda::updatePressure(Point_type pointType, const Grid& prcGrid,
		const Grid& coefGrid, const Grid& laplaGrid)
{
	printDebug(FULL_DEBUG, "In Grid_Cuda::updatePressure") ;

	// TO DO
	// return(Grid::updatePressure(pointType, prcGrid, coefGrid, laplaGrid)) ;

	//pointType
	Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
	Grid::getGridIndex(pointType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End);

	Myfloat *prcGrid_d_grid_3d = ((Grid_Cuda&) prcGrid).d_grid_3d ;
	Myfloat *coefGrid_d_grid_3d = ((Grid_Cuda&) coefGrid).d_grid_3d ;
	Myfloat *laplaGrid_d_grid_3d = ((Grid_Cuda&) laplaGrid).d_grid_3d ;
	kernel_updatePressure<<<1024,256>>>(d_grid_3d, prcGrid_d_grid_3d, coefGrid_d_grid_3d, laplaGrid_d_grid_3d, n1, n2, n3, i1Start, i1End, i2Start, i2End, i3Start, i3End);

	cudaDeviceSynchronize();

	printDebug(FULL_DEBUG, "Out Grid_Cuda::updatePressure") ;
	return(RTN_CODE_OK) ;
}
//-------------------------------------------------------------------------------------------------------
Rtn_code Grid_Cuda::applyBoundaryCondition(BoundCond_type boundCondType)
{
	printDebug(FULL_DEBUG, "In Grid_Cuda::applyBoundaryCondition") ;

	// TO DO
	// return(Grid::applyBoundaryCondition(boundCondType)) ;

	if (boundCondType != BOUND_COND_ANTI_MIRROR)
	{
		printError("CUDA: only BOUND_COND_ANTI_MIRROR boundary condition for now");
	}

	Myint64 	i1halo1_i1Start, i1halo1_i1End, i1halo1_i2Start, i1halo1_i2End, i1halo1_i3Start, i1halo1_i3End,
				i1halo2_i1Start, i1halo2_i1End, i1halo2_i2Start, i1halo2_i2End, i1halo2_i3Start, i1halo2_i3End,
				i2halo1_i1Start, i2halo1_i1End, i2halo1_i2Start, i2halo1_i2End, i2halo1_i3Start, i2halo1_i3End,
				i2halo2_i1Start, i2halo2_i1End, i2halo2_i2Start, i2halo2_i2End, i2halo2_i3Start, i2halo2_i3End,
				i3halo1_i1Start, i3halo1_i1End, i3halo1_i2Start, i3halo1_i2End, i3halo1_i3Start, i3halo1_i3End,
				i3halo2_i1Start, i3halo2_i1End, i3halo2_i2Start, i3halo2_i2End, i3halo2_i3Start, i3halo2_i3End;
	
	Grid::getGridIndex(I1HALO1, &i1halo1_i1Start, &i1halo1_i1End, &i1halo1_i2Start, &i1halo1_i2End, &i1halo1_i3Start, &i1halo1_i3End);
	Grid::getGridIndex(I1HALO2, &i1halo2_i1Start, &i1halo2_i1End, &i1halo2_i2Start, &i1halo2_i2End, &i1halo2_i3Start, &i1halo2_i3End);
	Grid::getGridIndex(I2HALO1, &i2halo1_i1Start, &i2halo1_i1End, &i2halo1_i2Start, &i2halo1_i2End, &i2halo1_i3Start, &i2halo1_i3End);
	Grid::getGridIndex(I2HALO2, &i2halo2_i1Start, &i2halo2_i1End, &i2halo2_i2Start, &i2halo2_i2End, &i2halo2_i3Start, &i2halo2_i3End);
	Grid::getGridIndex(I3HALO1, &i3halo1_i1Start, &i3halo1_i1End, &i3halo1_i2Start, &i3halo1_i2End, &i3halo1_i3Start, &i3halo1_i3End);
	Grid::getGridIndex(I3HALO2, &i3halo2_i1Start, &i3halo2_i1End, &i3halo2_i2Start, &i3halo2_i2End, &i3halo2_i3Start, &i3halo2_i3End);


	kernel_applyBoundaryCondition<<<1024,256>>>(d_grid_3d, n1, n2, n3, 	
		i1halo1_i1Start, i1halo1_i1End, i1halo1_i2Start, i1halo1_i2End, i1halo1_i3Start, i1halo1_i3End,
		i1halo2_i1Start, i1halo2_i1End, i1halo2_i2Start, i1halo2_i2End, i1halo2_i3Start, i1halo2_i3End,
		i2halo1_i1Start, i2halo1_i1End, i2halo1_i2Start, i2halo1_i2End, i2halo1_i3Start, i2halo1_i3End,
		i2halo2_i1Start, i2halo2_i1End, i2halo2_i2Start, i2halo2_i2End, i2halo2_i3Start, i2halo2_i3End,
		i3halo1_i1Start, i3halo1_i1End, i3halo1_i2Start, i3halo1_i2End, i3halo1_i3Start, i3halo1_i3End,
		i3halo2_i1Start, i3halo2_i1End, i3halo2_i2Start, i3halo2_i2End, i3halo2_i3Start, i3halo2_i3End);

	cudaDeviceSynchronize();

	printDebug(FULL_DEBUG, "Out Grid_Cuda::applyBoundaryCondition") ;
	return(RTN_CODE_OK) ;
}

Myfloat Grid_Cuda::getSumAbs(Point_type) const
{
	printWarning("getSumAbs not yet implemented on the GPU, so returning 0.1");
	return 0.1;
}

Myfloat Grid_Cuda::getSumAbsDiff(Point_type, const Grid&) const
{
	printWarning("getSumAbsDiff not yet implemented on the GPU, so returning 0");
	return 0.0;
}
} // namespace hpcscan
