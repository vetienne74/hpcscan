#include "hip/hip_runtime.h"

//-------------------------------------------------------------------------------------------------------
// This grid is activated with command line option -testMode HIP
// Derived class from Grid
// HIP implementation (target GPU)
//-------------------------------------------------------------------------------------------------------

#include "grid_Hip.h"

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

#include <hip/hip_cooperative_groups.h>
// #include <cooperative_groups/reduce.h>

namespace cg = cooperative_groups;

using namespace std;

namespace hpcscan {

static const int blockSize = 512 ;

//-------------------------------------------------------------------------------------------------------

//Macro for checking cuda errors following a cuda launch or api call
#define cudaCheckError() {                                          \
	hipError_t e=hipGetLastError();                                 \
	if(e!=hipSuccess) {                                              \
	  printf("Hip failure %s:%d: '%s'\n",__FILE__,__LINE__,hipGetErrorString(e));           \
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

__global__ void kernel_diff(Myfloat *data1, Myfloat *data2, Myfloat *dataOut,
		int n1, int n2, int n3, Myint64 i1Start, Myint64 i1End, Myint64 i2Start, Myint64 i2End, Myint64 i3Start, Myint64 i3End)
{
	int size = n1*n2*n3;
	int tid = threadIdx.x + blockIdx.x*blockDim.x;

	cg::thread_block cta = cg::this_thread_block();
	__shared__ float sdata[256];

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
	__shared__ float sdata[256];
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

__global__ void kernel_updatePressure(Myfloat *prn, Myfloat *prc, Myfloat *coef, Myfloat *lapla,
		int n1, int n2, int n3, Myint64 i1Start, Myint64 i1End, Myint64 i2Start, Myint64 i2End, Myint64 i3Start, Myint64 i3End)
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

Grid_Hip::Grid_Hip(Grid_type gridTypeIn) : Grid(gridTypeIn)
														{
	printDebug(MID_DEBUG, "IN Grid_Hip::Grid_Hip");

	gridMode = "Hip" ;
	
	d_grid_3d = NULL;
	d_help_3d = NULL;

	printDebug(MID_DEBUG, "OUT Grid_Hip::Grid_Hip");
														}

//-------------------------------------------------------------------------------------------------------

Grid_Hip::Grid_Hip(Grid_type gridTypeIn, Dim_type dimIn,
		Myint64 n1InnerIn, Myint64 n2InnerIn, Myint64 n3InnerIn) : Grid(gridTypeIn, dimIn,
				n1InnerIn, n2InnerIn, n3InnerIn)
{
	printDebug(MID_DEBUG, "IN Grid_Hip::Grid_Hip");

	gridMode = "Hip" ;
	
	d_grid_3d = NULL;
	d_help_3d = NULL;

	printDebug(MID_DEBUG, "OUT Grid_Hip::Grid_Hip");
}

//-------------------------------------------------------------------------------------------------------

Grid_Hip::~Grid_Hip(void)
{
	printDebug(MID_DEBUG, "IN Grid_Hip::~Grid_Hip");

	//delete[] grid_3d ;
	hipFree(d_grid_3d);
	hipFree(d_help_3d);
	cudaCheckError();

	printDebug(MID_DEBUG, "OUT Grid_Hip::~Grid_Hip");
}

//-------------------------------------------------------------------------------------------------------

void Grid_Hip::info(void)
{
	printDebug(FULL_DEBUG, "IN Grid_Hip::info");

	// parent class info
	Grid::info() ;

	// additional info
	printInfo(MASTER, "") ;
	printInfo(MASTER, " * GPU parameters * ") ;
	printInfo(MASTER, " Blocks per grid", 999) ;
	printInfo(MASTER, " Threads per block", 999) ;

	int startDevice = 0;
	int endDevice = 0;
	int deviceCount;
	hipError_t error_id = hipGetDeviceCount(&deviceCount);
	if (error_id != hipSuccess) {
		printError(" In Grid_Hip::info, hipGetDeviceCount", (int) error_id) ;
	}

	if (deviceCount == 0) {
		printError(" No GPU found") ;
	}
	else
	{
		printInfo(MASTER, " Number of GPUs found", deviceCount) ;
	}

	startDevice = 0;
	endDevice = deviceCount - 1;

	for (int currentDevice = startDevice; currentDevice <= endDevice;
			currentDevice++) {
		hipDeviceProp_t deviceProp;
		hipError_t error_id = hipGetDeviceProperties(&deviceProp, currentDevice);

		if (error_id == hipSuccess) {
			string deviceStr = " Device #" + to_string(currentDevice) + "\t";
			printInfo(MASTER, deviceStr, deviceProp.name) ;

			if (deviceProp.computeMode == hipComputeModeProhibited) {
				printError(" Error: device is running in <Compute Mode Prohibited>") ;
			}
		} else {
			printf("hipGetDeviceProperties returned %d\n-> %s\n", (int)error_id,
					hipGetErrorString(error_id));
		}
	}

	printDebug(FULL_DEBUG, "IN Grid_Hip::info");
}


//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_Hip::FD_LAPLACIAN(Point_type pointType, const Grid& Wgrid, Myint fdOrder)
{
	printDebug(MID_DEBUG, "IN Grid_Hip::FD_LAPLACIAN");

	// TO DO

	printDebug(MID_DEBUG, "OUT Grid_Hip::FD_LAPLACIAN");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_Hip::computePressureWithFD(Grid& prcGridIn, Grid& coefGridIn, Myint fdOrder)
{

	printDebug(FULL_DEBUG, "In Grid_Hip::computePressureWithFD") ;

	//pointType
	Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
	getGridIndex(INNER_POINTS, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End);

	if (fdOrder != 4 && fdOrder != 8)
	{
		printError("HIP: only FD order 4/8 implemented");
	}

	const Myfloat inv_d1  = Myfloat(1.0) / d1 ;
	const Myfloat inv_d2  = Myfloat(1.0) / d2 ;
	const Myfloat inv_d3  = Myfloat(1.0) / d3 ;

	const Myfloat inv2_d1 = inv_d1 * inv_d1 ;
	const Myfloat inv2_d2 = inv_d2 * inv_d2 ;
	const Myfloat inv2_d3 = inv_d3 * inv_d3 ;

	Myfloat *prc_d_grid_3d = ((Grid_Hip&) prcGridIn).d_grid_3d ;
	Myfloat *coef_d_grid_3d = ((Grid_Hip&) coefGridIn).d_grid_3d ;

	if (fdOrder == 4)
	{
		hipLaunchKernelGGL(kernel_computePressureWithFD_O4, dim3(1024), dim3(256), 0, 0, d_grid_3d, prc_d_grid_3d, coef_d_grid_3d,inv2_d1,inv2_d2,inv2_d3,n1,n2,n3,i1Start,i1End,i2Start,i2End,i3Start,i3End);
	}
	else if (fdOrder == 8)
	{
		hipLaunchKernelGGL(kernel_computePressureWithFD_O8, dim3(1024), dim3(256), 0, 0, d_grid_3d, prc_d_grid_3d, coef_d_grid_3d,inv2_d1,inv2_d2,inv2_d3,n1,n2,n3,i1Start,i1End,i2Start,i2End,i3Start,i3End);
	}
	
	cudaCheckError();
	hipDeviceSynchronize();
	
	printDebug(FULL_DEBUG, "Out Grid_Hip::computePressureWithFD") ;
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

void Grid_Hip::initializeGrid(void)
{
	printDebug(FULL_DEBUG, "In Grid_Hip::initializeGrid") ;

	Grid::initializeGrid() ; // this sets up halos etc.

	if (d_grid_3d == NULL)
	{
		hipMalloc( (void**)&d_grid_3d, n1*n2*n3*sizeof(Myfloat) );
		cudaCheckError();

		hipMalloc( (void**)&d_help_3d, 1024*sizeof(Myfloat) );
		cudaCheckError();
	}
	printDebug(FULL_DEBUG, "Out Grid_Hip::initializeGrid") ;
}

//-------------------------------------------------------------------------------------------------------
void Grid_Hip::fill(Point_type pointType, Myfloat val)
{
	printDebug(FULL_DEBUG, "In Grid_Hip::fill") ;

	//pointType
	Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
	getGridIndex(pointType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End);

	hipLaunchKernelGGL(kernel_fill_const, dim3(2048), dim3(128), 0, 0, d_grid_3d,val,n1,n2,n3,i1Start,i1End,i2Start,i2End,i3Start,i3End);
	hipDeviceSynchronize();
	cudaCheckError();

	printDebug(FULL_DEBUG, "Out Grid_Hip::fill") ;
}

//-------------------------------------------------------------------------------------------------------
void Grid_Hip::fill(Point_type pointType, Func_type t1,  Func_type t2, Func_type t3,
		Myfloat64 param1, Myfloat64 param2, Myfloat64 param3, Myfloat64 amp)
{
	printDebug(FULL_DEBUG, "In Grid_Hip::fill") ;

	//pointType
	Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
	getGridIndex(pointType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End);

	// we only use sine and linear for now
	int ok = 1;
	if (t1 != FUNC_SINE && t1 != FUNC_LINEAR) ok = 0;
	if (t2 != FUNC_SINE && t2 != FUNC_LINEAR) ok = 0;
	if (t3 != FUNC_SINE && t3 != FUNC_LINEAR) ok = 0;
	if (!ok) printError("HIP: only FUNC_SINE and FUNC_LINEAR implemented");

	if ((t1==t2 && t2==t3)==false) printError("HIP: func has to be same in each dimension");

	if (t1 == FUNC_SINE) hipLaunchKernelGGL(kernel_fill_sine, dim3(1024), dim3(128), 0, 0, d_grid_3d,param1,param2,param3,amp,n1,n2,n3,i1Start,i1End,i2Start,i2End,i3Start,i3End,Orig1,Orig2,Orig3,d1,d2,d2);
	else                 hipLaunchKernelGGL(kernel_fill_linear, dim3(1024), dim3(128), 0, 0, d_grid_3d,param1,param2,param3,amp,n1,n2,n3,i1Start,i1End,i2Start,i2End,i3Start,i3End,Orig1,Orig2,Orig3,d1,d2,d2);


	printDebug(FULL_DEBUG, "Out Grid_Hip::fill") ;
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

void Grid_Hip::fillArray(Myfloat val)
{
	printDebug(MID_DEBUG, "IN Grid_Hip::fillArray");

	Myint64 gridSize = n1*n2*n3;
	hipLaunchKernelGGL(kernel_fillArray, dim3(1024), dim3(256), 0, 0, d_grid_3d, val, gridSize) ;

	hipDeviceSynchronize();

	printDebug(MID_DEBUG, "OUT Grid_Hip::fillArray");
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

void Grid_Hip::copyArray(const Grid& gridIn)
{
	printDebug(MID_DEBUG, "IN Grid_Hip::copyArray");

	Myfloat *gridIn_d_grid_3d = ((Grid_Hip&) gridIn).d_grid_3d ;
	Myint64 gridSize = n1*n2*n3;
	hipLaunchKernelGGL(kernel_copyArray, dim3(1024), dim3(256), 0, 0, d_grid_3d, gridIn_d_grid_3d, gridSize) ;

	hipDeviceSynchronize();

	printDebug(MID_DEBUG, "OUT Grid_Hip::copyArray");
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

void Grid_Hip::addArray(const Grid& gridIn1, const Grid& gridIn2)
{
	printDebug(MID_DEBUG, "IN Grid_Hip::addArray");

	Myfloat *gridIn1_d_grid_3d = ((Grid_Hip&) gridIn1).d_grid_3d ;
	Myfloat *gridIn2_d_grid_3d = ((Grid_Hip&) gridIn2).d_grid_3d ;
	Myint64 gridSize = n1*n2*n3;
	hipLaunchKernelGGL(kernel_addArray, dim3(1024), dim3(256), 0, 0, d_grid_3d, gridIn1_d_grid_3d, gridIn2_d_grid_3d, gridSize) ;

	hipDeviceSynchronize();

	printDebug(MID_DEBUG, "OUT Grid_Hip::addArray");
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

void Grid_Hip::multiplyArray(const Grid& gridIn1, const Grid& gridIn2)
{
	printDebug(MID_DEBUG, "IN Grid_Hip::multiplyArray");

	Myfloat *gridIn1_d_grid_3d = ((Grid_Hip&) gridIn1).d_grid_3d ;
	Myfloat *gridIn2_d_grid_3d = ((Grid_Hip&) gridIn2).d_grid_3d ;
	Myint64 gridSize = n1*n2*n3;
	hipLaunchKernelGGL(kernel_multiplyArray, dim3(1024), dim3(256), 0, 0, d_grid_3d, gridIn1_d_grid_3d, gridIn2_d_grid_3d, gridSize) ;

	hipDeviceSynchronize();

	printDebug(MID_DEBUG, "OUT Grid_Hip::multiplyArray");
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

void Grid_Hip::addUpdateArray(const Grid& gridIn)
{
	printDebug(MID_DEBUG, "IN Grid_Hip::addUpdateArray");

	Myfloat *gridIn_d_grid_3d = ((Grid_Hip&) gridIn).d_grid_3d ;
	Myint64 gridSize = n1*n2*n3;
	hipLaunchKernelGGL(kernel_addUpdateArray, dim3(1024), dim3(256), 0, 0, d_grid_3d, gridIn_d_grid_3d, gridSize) ;

	hipDeviceSynchronize();

	printDebug(MID_DEBUG, "OUT Grid_Hip::addUpdateArray");
}

//-------------------------------------------------------------------------------------------------------

__global__ void minCommSingleBlock(const Myfloat *a, Myfloat *out,
		int n1, int n2, int n3, Myint64 i1Start, Myint64 i1End, Myint64 i2Start, Myint64 i2End, Myint64 i3Start, Myint64 i3End)
{
	int idx = threadIdx.x;
	Myfloat minval = +FLT_MAX;
	Myint64 arraySize = n1*n2*n3;

	for (Myint64 i = idx; i < arraySize; i += blockSize)
	{
		int i3 = i / (n1*n2);
		int idx2 = i -i3*n1*n2;
		int i2 = idx2/n1;
		int i1 = idx2%n1;

		if (i1 >= i1Start && i1 <= i1End &&
				i2 >= i2Start && i2 <= i2End &&
				i3 >= i3Start && i3 <= i3End   )
		{
			if (a[i] < minval) minval = a[i] ;
		}
	}

	__shared__ Myfloat r[blockSize];
	r[idx] = minval;
	__syncthreads();
	for (int size = blockSize/2; size>0; size/=2) { //uniform
		if (idx<size)
			if (r[idx+size] < r[idx]) r[idx] = r[idx+size];
		__syncthreads();
	}
	if (idx == 0)
		*out = r[0];
}

Myfloat Grid_Hip::getMin(Point_type pointType)
{
	printDebug(FULL_DEBUG, "In Grid_Hip::getMin") ;

	Myfloat val = 0 ;
	Myfloat *d_val ;
	hipMalloc((void**)&d_val, sizeof(Myfloat) * 1);

	Myfloat *grid_d_grid_3d = d_grid_3d ;
	Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
	getGridIndex(pointType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End);
	hipLaunchKernelGGL(minCommSingleBlock, dim3(1), dim3(blockSize), 0, 0, grid_d_grid_3d, d_val,
			n1, n2, n3, i1Start, i1End, i2Start, i2End, i3Start, i3End);
	hipDeviceSynchronize();
	hipMemcpy(&val, d_val, sizeof(Myfloat), hipMemcpyDeviceToHost);
	cudaCheckError();
	hipFree(d_val);

	printDebug(FULL_DEBUG, "Out Grid_Hip::getMin") ;

	return val ;
}

//-------------------------------------------------------------------------------------------------------

__global__ void maxCommSingleBlock(const Myfloat *a, Myfloat *out,
		int n1, int n2, int n3, Myint64 i1Start, Myint64 i1End, Myint64 i2Start, Myint64 i2End, Myint64 i3Start, Myint64 i3End)
{
	int idx = threadIdx.x;
	Myfloat maxval = -FLT_MAX;
	Myint64 arraySize = n1*n2*n3;

	for (Myint64 i = idx; i < arraySize; i += blockSize)
	{
		int i3 = i / (n1*n2);
		int idx2 = i -i3*n1*n2;
		int i2 = idx2/n1;
		int i1 = idx2%n1;

		if (i1 >= i1Start && i1 <= i1End &&
				i2 >= i2Start && i2 <= i2End &&
				i3 >= i3Start && i3 <= i3End   )
		{
			if (a[i] > maxval) maxval = a[i] ;
		}
	}

	__shared__ Myfloat r[blockSize];
	r[idx] = maxval;
	__syncthreads();
	for (int size = blockSize/2; size>0; size/=2) { //uniform
		if (idx<size)
			if (r[idx+size] > r[idx]) r[idx] = r[idx+size];
		__syncthreads();
	}
	if (idx == 0)
		*out = r[0];
}

Myfloat Grid_Hip::getMax(Point_type pointType)
{
	printDebug(FULL_DEBUG, "In Grid_Hip::getMax") ;

	Myfloat val = 0 ;
	Myfloat *d_val ;
	hipMalloc((void**)&d_val, sizeof(Myfloat) * 1);

	Myfloat *grid_d_grid_3d = d_grid_3d ;
	Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
	getGridIndex(pointType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End);
	hipLaunchKernelGGL(maxCommSingleBlock, dim3(1), dim3(blockSize), 0, 0, grid_d_grid_3d, d_val,
			n1, n2, n3, i1Start, i1End, i2Start, i2End, i3Start, i3End);
	hipDeviceSynchronize();
	hipMemcpy(&val, d_val, sizeof(Myfloat), hipMemcpyDeviceToHost);
	cudaCheckError();
	hipFree(d_val);

	printDebug(FULL_DEBUG, "Out Grid_Hip::getMax") ;

	return val ;
}


//-------------------------------------------------------------------------------------------------------

Myfloat Grid_Hip::L1Err(Point_type pointType, const Grid& gridIn) const
{
	printDebug(FULL_DEBUG, "In Grid_Hip::L1Err") ;

	Myfloat64 sum1 , sum2 ;
	sum1 = getSumAbsDiff(pointType, gridIn) ;
	sum2 = gridIn.getSumAbs(pointType) ;

	// prevent divide by zero
	if (sum2 < MAX_ERR_FLOAT) sum2 = 1.0 * npoint ;
	Myfloat err = sum1 / sum2 ;

	printDebug(LIGHT_DEBUG, "sum1", sum1) ;
	printDebug(LIGHT_DEBUG, "sum2", sum2) ;
	printDebug(LIGHT_DEBUG, "err", err) ;

	if (std::isnan(err))
	{
		printError("In Grid_Hip::L1Err, std::isnan(err)") ;
	}

	printDebug(FULL_DEBUG, "Out Grid_Hip::L1Err") ;
	return(err) ;
}

//-------------------------------------------------------------------------------------------------------

Myfloat Grid_Hip::allProcL1Err(Point_type pointType, const Grid& gridIn) const
{
	printDebug(LIGHT_DEBUG, "IN Grid_Hip::allProcL1Err");

	Myfloat64 sum1Loc = 0.0 ;
	Myfloat64 sum2Loc = 0.0 ;
	Myfloat64 sum1 = 0.0 ;
	Myfloat64 sum2 = 0.0 ;

	sum1Loc = getSumAbsDiff(pointType, gridIn) ;
	sum2Loc = gridIn.getSumAbs(pointType) ;

	// reduction
	MPI_Reduce(&sum1Loc, &sum1, 1, MPI_MYFLOAT64, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&sum2Loc, &sum2, 1, MPI_MYFLOAT64, MPI_SUM, 0, MPI_COMM_WORLD);

	// prevent divide by zero
	if (sum2 == 0.0) sum2 = 1.0 * npoint ;
	Myfloat err = sum1 / sum2 ;

	printDebug(LIGHT_DEBUG, "sum1", sum1) ;
	printDebug(LIGHT_DEBUG, "sum2", sum2) ;
	printDebug(LIGHT_DEBUG, "err", err) ;

	if (std::isnan(err))
	{
		printError("In Grid_Hip::allProcL1Err, std::isnan(err)") ;
	}

	printDebug(LIGHT_DEBUG, "OUT Grid_Hip::allProcL1Err");
	return(err) ;
}

//-------------------------------------------------------------------------------------------------------
Rtn_code Grid_Hip::updatePressure(Point_type pointType, const Grid& prcGrid,
		const Grid& coefGrid, const Grid& laplaGrid)
{
	printDebug(FULL_DEBUG, "In Grid_Hip::updatePressure") ;

	//pointType
	Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
	getGridIndex(pointType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End);

	Myfloat *prcGrid_d_grid_3d = ((Grid_Hip&) prcGrid).d_grid_3d ;
	Myfloat *coefGrid_d_grid_3d = ((Grid_Hip&) coefGrid).d_grid_3d ;
	Myfloat *laplaGrid_d_grid_3d = ((Grid_Hip&) laplaGrid).d_grid_3d ;
	hipLaunchKernelGGL(kernel_updatePressure, dim3(1024), dim3(256), 0, 0, d_grid_3d, prcGrid_d_grid_3d, coefGrid_d_grid_3d, laplaGrid_d_grid_3d,
			n1, n2, n3, i1Start, i1End, i2Start, i2End, i3Start, i3End);

	hipDeviceSynchronize();

	printDebug(FULL_DEBUG, "Out Grid_Hip::updatePressure") ;
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_Hip::exchangeHalo(MPI_comm_mode_type commMode, Point_type pointType)
{
	printDebug(FULL_DEBUG, "IN Grid_Hip::exchangeHalo");

	if (commMode == MPI_COMM_MODE_SENDRECV)
	{
		printDebug(FULL_DEBUG, "MPI_COMM_MODE_SENDRECV") ;

		MPI_Status status ;
		Myint rankSend, rankRecv ;
		MPI_Datatype typeSend, typeRecv ;
		Myfloat *bufSend, *bufRecv ;

		Point_type bufSendPointType, bufRecvPointType ;

		if (pointType == I1HALO1)
		{
			printDebug(FULL_DEBUG, "I1HALO1") ;
			Myint64 i1Send = i1InnerEnd - haloWidth + 1 ;
			Myint64 i2Send = i2InnerStart ;
			Myint64 i3Send = i3InnerStart ;
			Myint64 i1Recv = i1Halo1Start ;
			Myint64 i2Recv = i2InnerStart ;
			Myint64 i3Recv = i3InnerStart ;
			Myint64 idxSend = i1Send+i2Send*n1+i3Send*n1*n2 ;
			Myint64 idxRecv = i1Recv+i2Recv*n1+i3Recv*n1*n2 ;
			bufSend = &(grid_3d[idxSend]) ;
			bufRecv = &(grid_3d[idxRecv]) ;
			typeSend = i1HaloDataType ;
			typeRecv = i1HaloDataType ;
			rankSend = i1ProcIdEnd ;
			rankRecv = i1ProcIdStart ;

			bufRecvPointType = I1HALO1 ;
			bufSendPointType = I1INNERHALO2 ;
		}

		else if (pointType == I1HALO2)
		{
			printDebug(FULL_DEBUG, "I1HALO2") ;
			Myint64 i1Send = i1InnerStart ;
			Myint64 i2Send = i2InnerStart ;
			Myint64 i3Send = i3InnerStart ;
			Myint64 i1Recv = i1Halo2Start ;
			Myint64 i2Recv = i2InnerStart ;
			Myint64 i3Recv = i3InnerStart ;
			Myint64 idxSend = i1Send+i2Send*n1+i3Send*n1*n2 ;
			Myint64 idxRecv = i1Recv+i2Recv*n1+i3Recv*n1*n2 ;
			bufSend = &(grid_3d[idxSend]) ;
			bufRecv = &(grid_3d[idxRecv]) ;
			typeSend = i1HaloDataType ;
			typeRecv = i1HaloDataType ;
			rankSend = i1ProcIdStart ;
			rankRecv = i1ProcIdEnd ;

			bufRecvPointType = I1HALO2 ;
			bufSendPointType = I1INNERHALO1 ;
		}

		else if (pointType == I2HALO1)
		{
			printDebug(FULL_DEBUG, "I2HALO1") ;
			Myint64 i1Send = i1InnerStart ;
			Myint64 i2Send = i2InnerEnd - haloWidth + 1 ;
			Myint64 i3Send = i3InnerStart ;
			Myint64 i1Recv = i1InnerStart ;
			Myint64 i2Recv = i2Halo1Start ;
			Myint64 i3Recv = i3InnerStart ;
			Myint64 idxSend = i1Send+i2Send*n1+i3Send*n1*n2 ;
			Myint64 idxRecv = i1Recv+i2Recv*n1+i3Recv*n1*n2 ;
			bufSend = &(grid_3d[idxSend]) ;
			bufRecv = &(grid_3d[idxRecv]) ;
			typeSend = i2HaloDataType ;
			typeRecv = i2HaloDataType ;
			rankSend = i2ProcIdEnd ;
			rankRecv = i2ProcIdStart ;

			bufRecvPointType = I2HALO1 ;
			bufSendPointType = I2INNERHALO2 ;
		}

		else if (pointType == I2HALO2)
		{
			printDebug(FULL_DEBUG, "I2HALO2") ;
			Myint64 i1Send = i1InnerStart ;
			Myint64 i2Send = i2InnerStart ;
			Myint64 i3Send = i3InnerStart ;
			Myint64 i1Recv = i1InnerStart ;
			Myint64 i2Recv = i2Halo2Start ;
			Myint64 i3Recv = i3InnerStart ;
			Myint64 idxSend = i1Send+i2Send*n1+i3Send*n1*n2 ;
			Myint64 idxRecv = i1Recv+i2Recv*n1+i3Recv*n1*n2 ;
			bufSend = &(grid_3d[idxSend]) ;
			bufRecv = &(grid_3d[idxRecv]) ;
			typeSend = i2HaloDataType ;
			typeRecv = i2HaloDataType ;
			rankSend = i2ProcIdStart ;
			rankRecv = i2ProcIdEnd ;

			bufRecvPointType = I2HALO2 ;
			bufSendPointType = I2INNERHALO1 ;
		}

		else if (pointType == I3HALO1)
		{
			printDebug(FULL_DEBUG, "I3HALO1") ;
			Myint64 i1Send = i1InnerStart ;
			Myint64 i2Send = i2InnerStart ;
			Myint64 i3Send = i3InnerEnd - haloWidth + 1 ;
			Myint64 i1Recv = i1InnerStart ;
			Myint64 i2Recv = i2InnerStart ;
			Myint64 i3Recv = i3Halo1Start ;
			Myint64 idxSend = i1Send+i2Send*n1+i3Send*n1*n2 ;
			Myint64 idxRecv = i1Recv+i2Recv*n1+i3Recv*n1*n2 ;
			bufSend = &(grid_3d[idxSend]) ;
			bufRecv = &(grid_3d[idxRecv]) ;
			typeSend = i3HaloDataType ;
			typeRecv = i3HaloDataType ;
			rankSend = i3ProcIdEnd ;
			rankRecv = i3ProcIdStart ;

			bufRecvPointType = I3HALO1 ;
			bufSendPointType = I3INNERHALO2 ;
		}

		else if (pointType == I3HALO2)
		{
			printDebug(FULL_DEBUG, "I3HALO2") ;
			Myint64 i1Send = i1InnerStart ;
			Myint64 i2Send = i2InnerStart ;
			Myint64 i3Send = i3InnerStart ;
			Myint64 i1Recv = i1InnerStart ;
			Myint64 i2Recv = i2InnerStart ;
			Myint64 i3Recv = i3Halo2Start ;
			Myint64 idxSend = i1Send+i2Send*n1+i3Send*n1*n2 ;
			Myint64 idxRecv = i1Recv+i2Recv*n1+i3Recv*n1*n2 ;
			bufSend = &(grid_3d[idxSend]) ;
			bufRecv = &(grid_3d[idxRecv]) ;
			typeSend = i3HaloDataType ;
			typeRecv = i3HaloDataType ;
			rankSend = i3ProcIdStart ;
			rankRecv = i3ProcIdEnd ;

			bufRecvPointType = I3HALO2 ;
			bufSendPointType = I3INNERHALO1 ;
		}

		else
		{
			printError("IN Grid_Hip::exchangeHalo, invalid pointType", pointType) ;
			return(RTN_CODE_KO) ;
		}

		// copy halo to send from device to host
		if (rankSend != MPI_PROC_NULL)
			copyGridDeviceToHost(bufSendPointType) ;

		// call MPI_Sendrecv
		printDebug(FULL_DEBUG, "MPI_Sendrecv", rankSend, rankRecv) ;
		MPI_Sendrecv(bufSend, 1, typeSend, rankSend, 0,
				bufRecv, 1, typeRecv, rankRecv, 0,
				MPI_COMM_WORLD, &status);

		// copy halo received from host to device
		if (rankRecv != MPI_PROC_NULL)
			copyGridHostToDevice(bufRecvPointType) ;

	}
	else
	{
		printError("IN Grid_Hip::exchangeHalo, invalid commMode", commMode) ;
		return(RTN_CODE_KO) ;
	}

	printDebug(FULL_DEBUG, "OUT Grid_Hip::exchangeHalo");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

void Grid_Hip::copyGridDeviceToHost(Point_type pointType)
{
	printDebug(FULL_DEBUG, "In Grid_Hip::copyGridDeviceToHost") ;

	Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
	getGridIndex(pointType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End);

	//if (pointType == I1INNERHALO1)
	for (Myint64 i3 = i3Start; i3<= i3End; i3++)
	{
		for (Myint64 i2 = i2Start; i2<= i2End; i2++)
		{
			for (Myint64 i1 = i1Start; i1<= i1End; i1++)
			{
				// TMP copy point to point
				Myint64 idx = i1+i2*n1+i3*n1*n2 ;
				hipMemcpy(&(grid_3d[idx]), &(d_grid_3d[idx]), sizeof(Myfloat), hipMemcpyDeviceToHost) ;
			}
		}
	}

	hipDeviceSynchronize();

	printDebug(FULL_DEBUG, "In Grid_Hip::copyGridHostToDevice") ;
}

//-------------------------------------------------------------------------------------------------------

void Grid_Hip::copyGridHostToDevice(Point_type pointType)
{
	printDebug(FULL_DEBUG, "In Grid_Hip::copyGridDeviceToHost") ;

	Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
	getGridIndex(pointType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End);

	//if (pointType == I1HALO1)
	for (Myint64 i3 = i3Start; i3<= i3End; i3++)
	{
		for (Myint64 i2 = i2Start; i2<= i2End; i2++)
		{
			for (Myint64 i1 = i1Start; i1<= i1End; i1++)
			{
				// TMP copy point to point
				Myint64 idx = i1+i2*n1+i3*n1*n2 ;
				hipMemcpy(&(d_grid_3d[idx]), &(grid_3d[idx]), sizeof(Myfloat), hipMemcpyHostToDevice) ;
			}
		}
	}

	hipDeviceSynchronize();

	printDebug(FULL_DEBUG, "In Grid_Hip::copyGridHostToDevice") ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_Hip::applyBoundaryCondition(BoundCond_type boundCondType)
{
	printDebug(FULL_DEBUG, "In Grid_Hip::applyBoundaryCondition") ;

	if (boundCondType != BOUND_COND_ANTI_MIRROR)
	{
		printError("HIP: only BOUND_COND_ANTI_MIRROR boundary condition for now");
	}

	Myint64 	i1halo1_i1Start, i1halo1_i1End, i1halo1_i2Start, i1halo1_i2End, i1halo1_i3Start, i1halo1_i3End,
				i1halo2_i1Start, i1halo2_i1End, i1halo2_i2Start, i1halo2_i2End, i1halo2_i3Start, i1halo2_i3End,
				i2halo1_i1Start, i2halo1_i1End, i2halo1_i2Start, i2halo1_i2End, i2halo1_i3Start, i2halo1_i3End,
				i2halo2_i1Start, i2halo2_i1End, i2halo2_i2Start, i2halo2_i2End, i2halo2_i3Start, i2halo2_i3End,
				i3halo1_i1Start, i3halo1_i1End, i3halo1_i2Start, i3halo1_i2End, i3halo1_i3Start, i3halo1_i3End,
				i3halo2_i1Start, i3halo2_i1End, i3halo2_i2Start, i3halo2_i2End, i3halo2_i3Start, i3halo2_i3End;
	
	getGridIndex(I1HALO1, &i1halo1_i1Start, &i1halo1_i1End, &i1halo1_i2Start, &i1halo1_i2End, &i1halo1_i3Start, &i1halo1_i3End);
	getGridIndex(I1HALO2, &i1halo2_i1Start, &i1halo2_i1End, &i1halo2_i2Start, &i1halo2_i2End, &i1halo2_i3Start, &i1halo2_i3End);
	getGridIndex(I2HALO1, &i2halo1_i1Start, &i2halo1_i1End, &i2halo1_i2Start, &i2halo1_i2End, &i2halo1_i3Start, &i2halo1_i3End);
	getGridIndex(I2HALO2, &i2halo2_i1Start, &i2halo2_i1End, &i2halo2_i2Start, &i2halo2_i2End, &i2halo2_i3Start, &i2halo2_i3End);
	getGridIndex(I3HALO1, &i3halo1_i1Start, &i3halo1_i1End, &i3halo1_i2Start, &i3halo1_i2End, &i3halo1_i3Start, &i3halo1_i3End);
	getGridIndex(I3HALO2, &i3halo2_i1Start, &i3halo2_i1End, &i3halo2_i2Start, &i3halo2_i2End, &i3halo2_i3Start, &i3halo2_i3End);


	hipLaunchKernelGGL(kernel_applyBoundaryCondition, dim3(1024), dim3(256), 0, 0, d_grid_3d, n1, n2, n3, 	
		i1halo1_i1Start, i1halo1_i1End, i1halo1_i2Start, i1halo1_i2End, i1halo1_i3Start, i1halo1_i3End,
		i1halo2_i1Start, i1halo2_i1End, i1halo2_i2Start, i1halo2_i2End, i1halo2_i3Start, i1halo2_i3End,
		i2halo1_i1Start, i2halo1_i1End, i2halo1_i2Start, i2halo1_i2End, i2halo1_i3Start, i2halo1_i3End,
		i2halo2_i1Start, i2halo2_i1End, i2halo2_i2Start, i2halo2_i2End, i2halo2_i3Start, i2halo2_i3End,
		i3halo1_i1Start, i3halo1_i1End, i3halo1_i2Start, i3halo1_i2End, i3halo1_i3Start, i3halo1_i3End,
		i3halo2_i1Start, i3halo2_i1End, i3halo2_i2Start, i3halo2_i2End, i3halo2_i3Start, i3halo2_i3End);

	hipDeviceSynchronize();

	printDebug(FULL_DEBUG, "Out Grid_Hip::applyBoundaryCondition") ;
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

__global__ void sumAbsCommSingleBlock(const Myfloat *a, Myfloat *out,
		int n1, int n2, int n3, Myint64 i1Start, Myint64 i1End, Myint64 i2Start, Myint64 i2End, Myint64 i3Start, Myint64 i3End)
{
	int idx = threadIdx.x;
	Myfloat64 sum = 0;
	Myint64 arraySize = n1*n2*n3;

	for (Myint64 i = idx; i < arraySize; i += blockSize)
	{
		int i3 = i / (n1*n2);
		int idx2 = i -i3*n1*n2;
		int i2 = idx2/n1;
		int i1 = i - i3*(n1*n2) - i2*n1;

		if (i1 >= i1Start && i1 <= i1End &&
				i2 >= i2Start && i2 <= i2End &&
				i3 >= i3Start && i3 <= i3End   )
		{
			sum += fabs(a[i]) ;
		}
	}

	__shared__ Myfloat r[blockSize];
	r[idx] = sum;
	__syncthreads();
	for (int size = blockSize/2; size>0; size/=2) { //uniform
		if (idx<size)
			r[idx] += r[idx+size];
		__syncthreads();
	}
	if (idx == 0)
		*out = r[0];
}

Myfloat Grid_Hip::getSumAbs(Point_type pointType) const
{
	printDebug(LIGHT_DEBUG, "IN Grid_Hip::getSumAbs");

	Myfloat sum = 0 ;
	Myfloat *d_sum ;
	hipMalloc((void**)&d_sum, sizeof(Myfloat) * 1);

	Myfloat *grid_d_grid_3d = d_grid_3d ;
	Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
	getGridIndex(pointType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End);
	hipLaunchKernelGGL(sumAbsCommSingleBlock, dim3(1), dim3(blockSize), 0, 0, grid_d_grid_3d, d_sum,
			n1, n2, n3, i1Start, i1End, i2Start, i2End, i3Start, i3End);
	hipDeviceSynchronize();
	hipMemcpy(&sum, d_sum, sizeof(Myfloat), hipMemcpyDeviceToHost);
	cudaCheckError();
	hipFree(d_sum);

	// reduction
	Myfloat64 sum2Loc = sum ;
	Myfloat64 sum2 = 0.0 ;
	if (gridType == GRID_LOCAL)
	{
		MPI_Reduce(&sum2Loc, &sum2, 1, MPI_MYFLOAT64, MPI_SUM, 0, MPI_COMM_WORLD);
	}
	else
	{
		sum2 = sum2Loc ;
	}

	printDebug(LIGHT_DEBUG, "sum2", sum2) ;

	if (std::isnan(sum2))
	{
		printError("In Grid_Hip::getSumAbs, std::isnan(sum2)") ;
	}

	printDebug(LIGHT_DEBUG, "OUT Grid_Hip::getSumAbs");

	return(sum2) ;
}

//-------------------------------------------------------------------------------------------------------

__global__ void sumAbsDiffCommSingleBlock(const Myfloat *a, const Myfloat *b, Myfloat *out,
		int n1, int n2, int n3, Myint64 i1Start, Myint64 i1End, Myint64 i2Start, Myint64 i2End, Myint64 i3Start, Myint64 i3End)
{
	int idx = threadIdx.x;
	Myfloat64 sum = 0;
	Myint64 arraySize = n1*n2*n3;

	for (Myint64 i = idx; i < arraySize; i += blockSize)
	{
		int i3 = i / (n1*n2);
		int idx2 = i -i3*n1*n2;
		int i2 = idx2/n1;
		int i1 = idx2%n1;

		if (i1 >= i1Start && i1 <= i1End &&
				i2 >= i2Start && i2 <= i2End &&
				i3 >= i3Start && i3 <= i3End   )
		{

			sum += fabs(a[i] - b[i]) ;
		}
	}

	__shared__ Myfloat r[blockSize];
	r[idx] = sum;
	__syncthreads();
	for (int size = blockSize/2; size>0; size/=2) { //uniform
		if (idx<size)
			r[idx] += r[idx+size];
		__syncthreads();
	}
	if (idx == 0)
		*out = r[0];
}

Myfloat Grid_Hip::getSumAbsDiff(Point_type pointType, const Grid& gridIn) const
{

	printDebug(LIGHT_DEBUG, "IN Grid_Hip::getSumAbsDiff");

	Myfloat sum = 0 ;
	Myfloat *d_sum ;
	hipMalloc((void**)&d_sum, sizeof(Myfloat) * 1);

	Myfloat *grid_d_grid_3d = d_grid_3d ;
	Myfloat *gridIn_d_grid_3d = ((Grid_Hip&) gridIn).d_grid_3d ;
	Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
	getGridIndex(pointType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End);
	hipLaunchKernelGGL(sumAbsDiffCommSingleBlock, dim3(1), dim3(blockSize), 0, 0, grid_d_grid_3d, gridIn_d_grid_3d, d_sum,
			n1, n2, n3, i1Start, i1End, i2Start, i2End, i3Start, i3End);
	hipDeviceSynchronize();
	hipMemcpy(&sum, d_sum, sizeof(Myfloat), hipMemcpyDeviceToHost);
	cudaCheckError();
	hipFree(d_sum);

	// reduction
	Myfloat64 sum1Loc = sum ;
	Myfloat64 sum1 = 0.0 ;
	if (gridType == GRID_LOCAL)
	{
		MPI_Reduce(&sum1Loc, &sum1, 1, MPI_MYFLOAT64, MPI_SUM, 0, MPI_COMM_WORLD);
	}
	else
	{
		sum1 = sum1Loc ;
	}

	printDebug(LIGHT_DEBUG, "sum1", sum1) ;

	if (std::isnan(sum1))
	{
		printError("In Grid_Hip::getSumAbsDiff, std::isnan(sum1)") ;
	}

	printDebug(LIGHT_DEBUG, "OUT Grid_Hip::getSumAbsDiff");

	return(sum1) ;
}

//-------------------------------------------------------------------------------------------------------

__global__ void maxErrCommSingleBlock(const Myfloat *a, const Myfloat *b, Myfloat *out,
		int n1, int n2, int n3, Myint64 i1Start, Myint64 i1End, Myint64 i2Start, Myint64 i2End, Myint64 i3Start, Myint64 i3End) {
	int idx = threadIdx.x;
	Myfloat err = -FLT_MAX, err2 = 0.0 ;
	Myint64 arraySize = n1*n2*n3;

	for (Myint64 i = idx; i < arraySize; i += blockSize)
	{
		int i3 = i / (n1*n2);
		int idx2 = i -i3*n1*n2;
		int i2 = idx2/n1;
		int i1 = idx2%n1;

		if (i1 >= i1Start && i1 <= i1End &&
				i2 >= i2Start && i2 <= i2End &&
				i3 >= i3Start && i3 <= i3End   )
		{
			if (fabs(b[i]) < MAX_ERR_FLOAT)
			{
				err2 = fabs(a[i] - b[i]) ;
			}
			else
			{
				err2 = fabs(a[i] - b[i]) / b[i] ;
			}

			if (err2 > err)
			{
				err = err2 ;
			}
		}
	}

	__shared__ Myfloat r[blockSize];
	r[idx] = err;
	__syncthreads();
	for (int size = blockSize/2; size>0; size/=2) { //uniform
		if (idx<size)
			if (r[idx+size] > r[idx]) r[idx] = r[idx+size];
		__syncthreads();
	}
	if (idx == 0)
		*out = r[0];
}

Myfloat Grid_Hip::maxErr(Point_type pointType, const Grid& gridIn) const
{
	printDebug(FULL_DEBUG, "IN Grid_Hip::maxErr");

	Myfloat err = 0 ;
	Myfloat *d_err ;
	hipMalloc((void**)&d_err, sizeof(Myfloat) * 1);

	Myfloat *grid_d_grid_3d = d_grid_3d ;
	Myfloat *gridIn_d_grid_3d = ((Grid_Hip&) gridIn).d_grid_3d ;
	Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
	getGridIndex(pointType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End);
	hipLaunchKernelGGL(maxErrCommSingleBlock, dim3(1), dim3(blockSize), 0, 0, grid_d_grid_3d, gridIn_d_grid_3d, d_err,
			n1, n2, n3, i1Start, i1End, i2Start, i2End, i3Start, i3End);
	hipDeviceSynchronize();
	hipMemcpy(&err, d_err, sizeof(Myfloat), hipMemcpyDeviceToHost);
	hipFree(d_err);
	cudaCheckError();

	printDebug(FULL_DEBUG, "OUT Grid_Hip::maxErr");
	return(err) ;
}

} // namespace hpcscan
