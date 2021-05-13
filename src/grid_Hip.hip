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

#include <hip/hip_cooperative_groups.h>

namespace cg = cooperative_groups;

using namespace std;

namespace hpcscan {

static const int blockSize = 256 ;

//-------------------------------------------------------------------------------------------------------

//Macro for checking cuda errors following a cuda launch or api call
#define cudaCheckError() {                                          \
	hipError_t e=hipGetLastError();                                 \
	if(e!=hipSuccess) {                                              \
	  printf("Hip failure %s:%d: '%s'\n",__FILE__,__LINE__,hipGetErrorString(e));           \
	  exit(0); \
	}                                                                 \
   }

//*******************************************************************************************************
// BEGINING OF HIP KERNELS
//*******************************************************************************************************

//-------------------------------------------------------------------------------------------------------
// retrieve minimum value (1st step)
// multi-block reduction on the input array dataIn
// each block finds its minimum and stores into the array dataOut at entry dataOut[blockIdx.x]

__global__ void kernel_multiBlk_minval(Myfloat *dataIn, Myfloat *dataOut,
		int n1, int n2, int n3, Myint64 i1Start, Myint64 i1End, Myint64 i2Start, Myint64 i2End, Myint64 i3Start, Myint64 i3End)
{
	Myint64 size = n1*n2*n3;
	Myint64 tid = threadIdx.x + blockIdx.x*blockDim.x;

	// dynamic shared memory
	HIP_DYNAMIC_SHARED( float, sdata)

	// set to max float value
	sdata[threadIdx.x] = +FLT_MAX ;

	// each thread find its minimum
	while (tid < size)
	{
		// convert 1d index to 3d indexes
		unsigned int i3 = tid / (n1*n2);
		unsigned int idx = tid-i3*n1*n2;
		unsigned int i2 = idx/n1;
		unsigned int i1 = idx%n1;

		// check if point fall into target area
		if (i1 >= i1Start && i1 <= i1End &&
			i2 >= i2Start && i2 <= i2End &&
			i3 >= i3Start && i3 <= i3End   )
		{
			// update min value
			Myfloat val = dataIn[tid];
			if (val < sdata[threadIdx.x]) sdata[threadIdx.x] = val;
		}

		tid += blockDim.x * gridDim.x;
	}

	__syncthreads();

	// find minimum between all threads
	for (unsigned int s = blockDim.x / 2; s > 0; s/=2)
	{
		if (threadIdx.x < s)
		{
			Myfloat val = sdata[threadIdx.x + s];
			if (val < sdata[threadIdx.x]) sdata[threadIdx.x] = val;
		}
		__syncthreads();
	}

	// write result for the block into global array
	if (threadIdx.x == 0) dataOut[blockIdx.x] = sdata[0];
}

//-------------------------------------------------------------------------------------------------------
// retrieve minimum value (2nd step)
// single block reduction on the input array dataInOut
// the minimum is stored at first entry dataInOut[0]

__global__ void kernel_singleBlk_minval(Myfloat *dataInOut, int dataInOutSize)
{
	int idx = threadIdx.x;
	for (int size = dataInOutSize/2; size>0; size/=2) {
		if (idx<size)
			if (dataInOut[idx+size] < dataInOut[idx]) dataInOut[idx] = dataInOut[idx+size];
		__syncthreads();
	}
}

//-------------------------------------------------------------------------------------------------------
// retrieve maximum value (1st step)
// multi-block reduction on the input array dataIn
// each block finds its maximum and stores into the array dataOut at entry dataOut[blockIdx.x]

__global__ void kernel_multiBlk_maxval(Myfloat *dataIn, Myfloat *dataOut,
		int n1, int n2, int n3, Myint64 i1Start, Myint64 i1End, Myint64 i2Start, Myint64 i2End, Myint64 i3Start, Myint64 i3End)
{
	Myint64 size = n1*n2*n3;
	Myint64 tid = threadIdx.x + blockIdx.x*blockDim.x;

	// dynamic shared memory
	HIP_DYNAMIC_SHARED( float, sdata)

	// set to min float value
	sdata[threadIdx.x] = -FLT_MAX ;

	// each thread find its minimum
	while (tid < size)
	{
		// convert 1d index to 3d indexes
		unsigned int i3 = tid / (n1*n2);
		unsigned int idx = tid-i3*n1*n2;
		unsigned int i2 = idx/n1;
		unsigned int i1 = idx%n1;

		// check if point fall into target area
		if (i1 >= i1Start && i1 <= i1End &&
			i2 >= i2Start && i2 <= i2End &&
			i3 >= i3Start && i3 <= i3End   )
		{
			// update max value
			Myfloat val = dataIn[tid];
			if (val > sdata[threadIdx.x]) sdata[threadIdx.x] = val;
		}

		tid += blockDim.x * gridDim.x;
	}

	__syncthreads();

	// find maximum between all threads
	for (unsigned int s = blockDim.x / 2; s > 0; s/=2)
	{
		if (threadIdx.x < s)
		{
			Myfloat val = sdata[threadIdx.x + s];
			if (val > sdata[threadIdx.x]) sdata[threadIdx.x] = val;
		}
		__syncthreads();
	}

	// write result for the block into global array
	if (threadIdx.x == 0) dataOut[blockIdx.x] = sdata[0];
}

//-------------------------------------------------------------------------------------------------------
// retrieve maximum value (2nd step)
// single block reduction on the input array dataInOut
// the maximum is stored at first entry dataInOut[0]

__global__ void kernel_singleBlk_maxval(Myfloat *dataInOut, int dataInOutSize)
{
	int idx = threadIdx.x;
	for (int size = dataInOutSize/2; size>0; size/=2) {
		if (idx<size)
			if (dataInOut[idx+size] > dataInOut[idx]) dataInOut[idx] = dataInOut[idx+size];
		__syncthreads();
	}
}

//-------------------------------------------------------------------------------------------------------
// sum abs values (1st step)
// multi-block reduction on the input array dataIn
// each block does the sum and stores into the array dataOut at entry dataOut[blockIdx.x]

__global__ void kernel_multiBlk_sumAbs(Myfloat *dataIn, Myfloat *dataOut,
		int n1, int n2, int n3, Myint64 i1Start, Myint64 i1End, Myint64 i2Start, Myint64 i2End, Myint64 i3Start, Myint64 i3End)
{
	Myint64 size = n1*n2*n3;
	Myint64 tid = threadIdx.x + blockIdx.x*blockDim.x;

	// dynamic shared memory
	HIP_DYNAMIC_SHARED( float, sdata)

	// set to zero
	sdata[threadIdx.x] = 0.0 ;

	// each thread sums
	while (tid < size)
	{
		// convert 1d index to 3d indexes
		unsigned int i3 = tid / (n1*n2);
		unsigned int idx = tid-i3*n1*n2;
		unsigned int i2 = idx/n1;
		unsigned int i1 = idx%n1;

		// check if point fall into target area
		if (i1 >= i1Start && i1 <= i1End &&
			i2 >= i2Start && i2 <= i2End &&
			i3 >= i3Start && i3 <= i3End   )
		{
			// update sum
			Myfloat val = fabs(dataIn[tid]) ;
			sdata[threadIdx.x] += val ;
		}

		tid += blockDim.x * gridDim.x;
	}

	__syncthreads();

	// sum between all threads
	for (unsigned int s = blockDim.x / 2; s > 0; s/=2)
	{
		if (threadIdx.x < s)
		{
			Myfloat val = sdata[threadIdx.x + s];
			sdata[threadIdx.x] += val;
		}
		__syncthreads();
	}

	// write result for the block into global array
	if (threadIdx.x == 0) dataOut[blockIdx.x] = sdata[0];
}

//-------------------------------------------------------------------------------------------------------
// sum values (absolute values) (2nd step)
// single block reduction on the input array dataInOut
// the maximum is stored at first entry dataInOut[0]

__global__ void kernel_singleBlk_sum(Myfloat *dataInOut, int dataInOutSize)
{
	int idx = threadIdx.x;
	for (int size = dataInOutSize/2; size>0; size/=2) {
		if (idx<size)
			dataInOut[idx] += dataInOut[idx+size];
		__syncthreads();
	}
}


//-------------------------------------------------------------------------------------------------------

// there is probably an easier way to implement this (3d blocks?)
__global__ void kernel_fill_const(Myfloat *data, Myfloat val, const int n1, const int n2, const int n3, Myint64 i1Start, Myint64 i1End, Myint64 i2Start, Myint64 i2End, Myint64 i3Start, Myint64 i3End)
{
	int size = n1*n2*n3;
	int tid = threadIdx.x + blockIdx.x*blockDim.x;

	// 1D to 3D index
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

__global__ void kernel_applyBoundaryCondition(Dim_type dim, Myfloat *data, int n1, int n2, int n3,
			Myint I1HALO1_neigh, Myint64 i1halo1_i1Start, Myint64 i1halo1_i1End, Myint64 i1halo1_i2Start, Myint64 i1halo1_i2End, Myint64 i1halo1_i3Start, Myint64 i1halo1_i3End,
			Myint I1HALO2_neigh, Myint64 i1halo2_i1Start, Myint64 i1halo2_i1End, Myint64 i1halo2_i2Start, Myint64 i1halo2_i2End, Myint64 i1halo2_i3Start, Myint64 i1halo2_i3End,
			Myint I2HALO1_neigh, Myint64 i2halo1_i1Start, Myint64 i2halo1_i1End, Myint64 i2halo1_i2Start, Myint64 i2halo1_i2End, Myint64 i2halo1_i3Start, Myint64 i2halo1_i3End,
			Myint I2HALO2_neigh, Myint64 i2halo2_i1Start, Myint64 i2halo2_i1End, Myint64 i2halo2_i2Start, Myint64 i2halo2_i2End, Myint64 i2halo2_i3Start, Myint64 i2halo2_i3End,
			Myint I3HALO1_neigh, Myint64 i3halo1_i1Start, Myint64 i3halo1_i1End, Myint64 i3halo1_i2Start, Myint64 i3halo1_i2End, Myint64 i3halo1_i3Start, Myint64 i3halo1_i3End,
			Myint I3HALO2_neigh, Myint64 i3halo2_i1Start, Myint64 i3halo2_i1End, Myint64 i3halo2_i2Start, Myint64 i3halo2_i2End, Myint64 i3halo2_i3Start, Myint64 i3halo2_i3End)
			
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
		if (I1HALO1_neigh == MPI_PROC_NULL)
		{
			Myint64 iInner = i1halo1_i1End+1;
			if (tid == iInner+i2*n1+i3*n1*n2) data[tid] = 0.0 ;

			if (i1 >= i1halo1_i1Start && i1 <= i1halo1_i1End &&
					i2 >= i1halo1_i2Start && i2 <= i1halo1_i2End &&
					i3 >= i1halo1_i3Start && i3 <= i1halo1_i3End   )
			{
				data[tid] = -data[(iInner+iInner-i1)+i2*n1+i3*n1*n2];
			}
		}

		// I1HALO2
		if (I1HALO2_neigh == MPI_PROC_NULL)
		{
			Myint64 iInner = i1halo2_i1Start-1;
			if (tid == iInner+i2*n1+i3*n1*n2) data[tid] = 0.0 ;

			if (i1 >= i1halo2_i1Start && i1 <= i1halo2_i1End &&
					i2 >= i1halo2_i2Start && i2 <= i1halo2_i2End &&
					i3 >= i1halo2_i3Start && i3 <= i1halo2_i3End   )
			{
				data[tid] = -data[(iInner-(i1-iInner))+i2*n1+i3*n1*n2];
			}
		}

		if (dim >= DIM2)
		{
			// I2HALO1
			if (I2HALO1_neigh == MPI_PROC_NULL)
			{
				Myint64 iInner = i2halo1_i2End+1;
				if (tid == i1+iInner*n1+i3*n1*n2) data[tid] = 0.0 ;

				if (i1 >= i2halo1_i1Start && i1 <= i2halo1_i1End &&
						i2 >= i2halo1_i2Start && i2 <= i2halo1_i2End &&
						i3 >= i2halo1_i3Start && i3 <= i2halo1_i3End   )
				{
					data[tid] = -data[i1+(iInner+iInner-i2)*n1+i3*n1*n2];
				}
			}

			// I2HALO2
			if (I2HALO2_neigh == MPI_PROC_NULL)
			{
				Myint64 iInner = i2halo2_i2Start-1;
				if (tid == i1+iInner*n1+i3*n1*n2) data[tid] = 0.0 ;

				if (i1 >= i2halo2_i1Start && i1 <= i2halo2_i1End &&
						i2 >= i2halo2_i2Start && i2 <= i2halo2_i2End &&
						i3 >= i2halo2_i3Start && i3 <= i2halo2_i3End   )
				{
					data[tid] = -data[i1+(iInner-(i2-iInner))*n1+i3*n1*n2];
				}
			}
		}

		if (dim >= DIM3)
		{
			// I3HALO1
			if (I3HALO1_neigh == MPI_PROC_NULL)
			{
				Myint64 iInner = i3halo1_i3End+1;
				if (tid == i1+i2*n1+iInner*n1*n2) data[tid] = 0.0 ;

				if (i1 >= i3halo1_i1Start && i1 <= i3halo1_i1End &&
						i2 >= i3halo1_i2Start && i2 <= i3halo1_i2End &&
						i3 >= i3halo1_i3Start && i3 <= i3halo1_i3End   )
				{
					data[tid] = -data[i1+i2*n1+(iInner+iInner-i3)*n1*n2];
				}
			}

			// I3HALO2
			if (I3HALO2_neigh == MPI_PROC_NULL)
			{
				Myint64 iInner = i3halo2_i3Start-1;
				if (tid == i1+i2*n1+iInner*n1*n2) data[tid] = 0.0 ;

				if (i1 >= i3halo2_i1Start && i1 <= i3halo2_i1End &&
						i2 >= i3halo2_i2Start && i2 <= i3halo2_i2End &&
						i3 >= i3halo2_i3Start && i3 <= i3halo2_i3End   )
				{
					data[tid] = -data[i1+i2*n1+(iInner-(i3-iInner))*n1*n2];
				}
			}
		}

		tid += blockDim.x * gridDim.x;
	}
}

//-------------------------------------------------------------------------------------------------------

__global__ void kernel_computePressureWithFD(Dim_type dim, Myint fdOrder, Myfloat *prn, Myfloat *prc, Myfloat *coef, Myfloat inv2_d1, Myfloat inv2_d2, Myfloat inv2_d3,
		int n1, int n2, int n3, Myint64 i1Start, Myint64 i1End, Myint64 i2Start, Myint64 i2End, Myint64 i3Start, Myint64 i3End)
{
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
			// compute FD for 1D
			if (dim == DIM1)
			{
				if (fdOrder == 2)
				{
					prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
							coef[i1+i2*n1+i3*n1*n2] *
							FD_D2_O2_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
				}
				else if (fdOrder == 4)
				{
					prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
							coef[i1+i2*n1+i3*n1*n2] *
							FD_D2_O4_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
				}
				else if (fdOrder == 8)
				{
					prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
							coef[i1+i2*n1+i3*n1*n2] *
							FD_D2_O8_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
				}
				else if (fdOrder == 12)
				{
					prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
							coef[i1+i2*n1+i3*n1*n2] *
							FD_D2_O12_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
				}
				else if (fdOrder == 16)
				{
					prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
							coef[i1+i2*n1+i3*n1*n2] *
							FD_D2_O16_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
				}
			}

			// compute FD for 2D
			else if (dim == DIM2)
			{
				if (fdOrder == 2)
				{
					prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
							coef[i1+i2*n1+i3*n1*n2] *
							(FD_D2_O2_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O2_N2(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)) ;
				}
				else if (fdOrder == 4)
				{
					prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
							coef[i1+i2*n1+i3*n1*n2] *
							(FD_D2_O4_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O4_N2(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)) ;
				}
				else if (fdOrder == 8)
				{
					prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
							coef[i1+i2*n1+i3*n1*n2] *
							(FD_D2_O8_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O8_N2(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)) ;
				}
				else if (fdOrder == 12)
				{
					prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
							coef[i1+i2*n1+i3*n1*n2] *
							(FD_D2_O12_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O12_N2(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)) ;
				}
				else if (fdOrder == 16)
				{
					prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
							coef[i1+i2*n1+i3*n1*n2] *
							(FD_D2_O16_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O16_N2(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)) ;
				}
			}

			// compute FD for 3D
			else if (dim == DIM3)
			{
				if (fdOrder == 2)
				{
					prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
							coef[i1+i2*n1+i3*n1*n2] *
							(FD_D2_O2_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O2_N2(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O2_N3(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)) ;
				}
				else if (fdOrder == 4)
				{
					prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
							coef[i1+i2*n1+i3*n1*n2] *
							(FD_D2_O4_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O4_N2(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O4_N3(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)) ;
				}
				else if (fdOrder == 8)
				{
					prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
							coef[i1+i2*n1+i3*n1*n2] *
							(FD_D2_O8_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O8_N2(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O8_N3(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)) ;
				}
				else if (fdOrder == 12)
				{
					prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
							coef[i1+i2*n1+i3*n1*n2] *
							(FD_D2_O12_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O12_N2(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O12_N3(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)) ;
				}
				else if (fdOrder == 16)
				{
					prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
							coef[i1+i2*n1+i3*n1*n2] *
							(FD_D2_O16_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O16_N2(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O16_N3(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)) ;
				}
			}
		}

		tid += blockDim.x * gridDim.x;
	}
}

//*******************************************************************************************************
// END OF HIP KERNELS
//*******************************************************************************************************

//-------------------------------------------------------------------------------------------------------

Grid_Hip::Grid_Hip(Grid_type gridTypeIn) : Grid(gridTypeIn)
														{
	printDebug(MID_DEBUG, "IN Grid_Hip::Grid_Hip");

	gridMode = GRID_MODE_HIP ;

	d_grid_3d = NULL;
	d_help_3d = NULL;

	gpuBlkSize  = Config::Instance()->gpuBlkSize ;
	gpuGridSize = Config::Instance()->gpuGridSize ;

	printDebug(MID_DEBUG, "OUT Grid_Hip::Grid_Hip");
														}

//-------------------------------------------------------------------------------------------------------

Grid_Hip::Grid_Hip(Grid_type gridTypeIn, Dim_type dimIn,
		Myint64 n1InnerIn, Myint64 n2InnerIn, Myint64 n3InnerIn) : Grid(gridTypeIn, dimIn,
				n1InnerIn, n2InnerIn, n3InnerIn)
{
	printDebug(MID_DEBUG, "IN Grid_Hip::Grid_Hip");

	gridMode = GRID_MODE_HIP ;

	d_grid_3d = NULL;
	d_help_3d = NULL;

	gpuBlkSize  = Config::Instance()->gpuBlkSize ;
	gpuGridSize = Config::Instance()->gpuGridSize ;


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
	printInfo(MASTER, " Blocks per grid", gpuGridSize) ;
	printInfo(MASTER, " Threads per block", gpuBlkSize) ;

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

void Grid_Hip::write(string file_name)
{
	printDebug(LIGHT_DEBUG, "IN Grid_Hip::write");

	// each proc write is own file

	if (Config::Instance()->writeGrid)
	{
		// copy grid from device to host
		copyGridDeviceToHost(ALL_POINTS) ;

		Grid::write(file_name) ;
	}

	printDebug(LIGHT_DEBUG, "OUT Grid_Hip::write");
}


//-------------------------------------------------------------------------------------------------------

__global__ void kernel_FD_D2_N1(Myint fdOrder, Myfloat *w, Myfloat *u, Myfloat inv2_d1, Myfloat inv2_d2, Myfloat inv2_d3,
		int n1, int n2, int n3, Myint64 i1Start, Myint64 i1End, Myint64 i2Start, Myint64 i2End, Myint64 i3Start, Myint64 i3End)
{
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
			// compute FD along N1
			if (fdOrder == 2)
			{
				w[i1+i2*n1+i3*n1*n2] =
						FD_D2_O2_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
			}
			else if (fdOrder == 4)
			{
				w[i1+i2*n1+i3*n1*n2] =
						FD_D2_O4_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
			}
			else if (fdOrder == 8)
			{
				w[i1+i2*n1+i3*n1*n2] =
						FD_D2_O8_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
			}
			else if (fdOrder == 12)
			{
				w[i1+i2*n1+i3*n1*n2] =
						FD_D2_O12_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
			}
			else if (fdOrder == 16)
			{
				w[i1+i2*n1+i3*n1*n2] =
						FD_D2_O16_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
			}
		}

		tid += blockDim.x * gridDim.x;
	}
}

Rtn_code Grid_Hip::FD_D2_N1(Point_type pointType, const Grid& Wgrid, Myint fdOrder)
{
	printDebug(MID_DEBUG, "IN Grid_Hip::FD_D2_N1");

	// check grids are same size
	if (this->sameSize(Wgrid) != true)
	{
		printError("Grid::FD_D2_N1, grids have not same size") ;
		return(RTN_CODE_KO) ;
	}

	Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
	getGridIndex(INNER_POINTS, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End) ;

	const Myfloat inv_d1  = Myfloat(1.0) / d1 ;
	const Myfloat inv_d2  = Myfloat(1.0) / d2 ;
	const Myfloat inv_d3  = Myfloat(1.0) / d3 ;

	const Myfloat inv2_d1 = inv_d1 * inv_d1 ;
	const Myfloat inv2_d2 = inv_d2 * inv_d2 ;
	const Myfloat inv2_d3 = inv_d3 * inv_d3 ;

	Myfloat * d_w = ((Grid_Hip&) Wgrid).d_grid_3d ;
	Myfloat * d_u = this->d_grid_3d ;

	hipLaunchKernelGGL(kernel_FD_D2_N1, dim3(1024), dim3(256), 0, 0, fdOrder, d_w, d_u,inv2_d1,inv2_d2,inv2_d3,
			n1,n2,n3,i1Start,i1End,i2Start,i2End,i3Start,i3End);

	cudaCheckError();
	hipDeviceSynchronize();

	printDebug(MID_DEBUG, "OUT Grid_Hip::FD_D2_N1");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

__global__ void kernel_FD_D2_N2(Myint fdOrder, Myfloat *w, Myfloat *u, Myfloat inv2_d1, Myfloat inv2_d2, Myfloat inv2_d3,
		int n1, int n2, int n3, Myint64 i1Start, Myint64 i1End, Myint64 i2Start, Myint64 i2End, Myint64 i3Start, Myint64 i3End)
{
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
			// compute FD along N2
			if (fdOrder == 2)
			{
				w[i1+i2*n1+i3*n1*n2] =
						FD_D2_O2_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
			}
			else if (fdOrder == 4)
			{
				w[i1+i2*n1+i3*n1*n2] =
						FD_D2_O4_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
			}
			else if (fdOrder == 8)
			{
				w[i1+i2*n1+i3*n1*n2] =
						FD_D2_O8_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
			}
			else if (fdOrder == 12)
			{
				w[i1+i2*n1+i3*n1*n2] =
						FD_D2_O12_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
			}
			else if (fdOrder == 16)
			{
				w[i1+i2*n1+i3*n1*n2] =
						FD_D2_O16_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
			}
		}

		tid += blockDim.x * gridDim.x;
	}
}

Rtn_code Grid_Hip::FD_D2_N2(Point_type pointType, const Grid& Wgrid, Myint fdOrder)
{
	printDebug(MID_DEBUG, "IN Grid_Hip::FD_D2_N2");

	// check grids are same size
	if (this->sameSize(Wgrid) != true)
	{
		printError("Grid::FD_D2_N2, grids have not same size") ;
		return(RTN_CODE_KO) ;
	}

	Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
	getGridIndex(INNER_POINTS, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End) ;

	const Myfloat inv_d1  = Myfloat(1.0) / d1 ;
	const Myfloat inv_d2  = Myfloat(1.0) / d2 ;
	const Myfloat inv_d3  = Myfloat(1.0) / d3 ;

	const Myfloat inv2_d1 = inv_d1 * inv_d1 ;
	const Myfloat inv2_d2 = inv_d2 * inv_d2 ;
	const Myfloat inv2_d3 = inv_d3 * inv_d3 ;

	Myfloat * d_w = ((Grid_Hip&) Wgrid).d_grid_3d ;
	Myfloat * d_u = this->d_grid_3d ;

	hipLaunchKernelGGL(kernel_FD_D2_N2, dim3(1024), dim3(256), 0, 0, fdOrder, d_w, d_u,inv2_d1,inv2_d2,inv2_d3,
			n1,n2,n3,i1Start,i1End,i2Start,i2End,i3Start,i3End);

	cudaCheckError();
	hipDeviceSynchronize();

	printDebug(MID_DEBUG, "OUT Grid_Hip::FD_D2_N2");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

__global__ void kernel_FD_D2_N3(Myint fdOrder, Myfloat *w, Myfloat *u, Myfloat inv2_d1, Myfloat inv2_d2, Myfloat inv2_d3,
		int n1, int n2, int n3, Myint64 i1Start, Myint64 i1End, Myint64 i2Start, Myint64 i2End, Myint64 i3Start, Myint64 i3End)
{
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
			// compute FD along N3
			if (fdOrder == 2)
			{
				w[i1+i2*n1+i3*n1*n2] =
						FD_D2_O2_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
			}
			else if (fdOrder == 4)
			{
				w[i1+i2*n1+i3*n1*n2] =
						FD_D2_O4_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
			}
			else if (fdOrder == 8)
			{
				w[i1+i2*n1+i3*n1*n2] =
						FD_D2_O8_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
			}
			else if (fdOrder == 12)
			{
				w[i1+i2*n1+i3*n1*n2] =
						FD_D2_O12_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
			}
			else if (fdOrder == 16)
			{
				w[i1+i2*n1+i3*n1*n2] =
						FD_D2_O16_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
			}
		}

		tid += blockDim.x * gridDim.x;
	}
}

Rtn_code Grid_Hip::FD_D2_N3(Point_type pointType, const Grid& Wgrid, Myint fdOrder)
{
	printDebug(MID_DEBUG, "IN Grid_Hip::FD_D2_N3");

	// check grids are same size
	if (this->sameSize(Wgrid) != true)
	{
		printError("Grid::FD_D2_N3, grids have not same size") ;
		return(RTN_CODE_KO) ;
	}

	Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
	getGridIndex(INNER_POINTS, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End) ;

	const Myfloat inv_d1  = Myfloat(1.0) / d1 ;
	const Myfloat inv_d2  = Myfloat(1.0) / d2 ;
	const Myfloat inv_d3  = Myfloat(1.0) / d3 ;

	const Myfloat inv2_d1 = inv_d1 * inv_d1 ;
	const Myfloat inv2_d2 = inv_d2 * inv_d2 ;
	const Myfloat inv2_d3 = inv_d3 * inv_d3 ;

	Myfloat * d_w = ((Grid_Hip&) Wgrid).d_grid_3d ;
	Myfloat * d_u = this->d_grid_3d ;

	hipLaunchKernelGGL(kernel_FD_D2_N3, dim3(1024), dim3(256), 0, 0, fdOrder, d_w, d_u,inv2_d1,inv2_d2,inv2_d3,
			n1,n2,n3,i1Start,i1End,i2Start,i2End,i3Start,i3End);

	cudaCheckError();
	hipDeviceSynchronize();

	printDebug(MID_DEBUG, "OUT Grid_Hip::FD_D2_N3");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

__global__ void kernel_FD_LAPLACIAN(Dim_type dim, Myint fdOrder, Myfloat *w, Myfloat *u, Myfloat inv2_d1, Myfloat inv2_d2, Myfloat inv2_d3,
		int n1, int n2, int n3, Myint64 i1Start, Myint64 i1End, Myint64 i2Start, Myint64 i2End, Myint64 i3Start, Myint64 i3End)
{
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
			// compute FD Laplacian for 1D
			if (dim == DIM1)
			{
				// same as FD_D2_N1
				if (fdOrder == 2)
				{
					w[i1+i2*n1+i3*n1*n2] =
							FD_D2_O2_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
				}
				else if (fdOrder == 4)
				{
					w[i1+i2*n1+i3*n1*n2] =
							FD_D2_O4_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
				}
				else if (fdOrder == 8)
				{
					w[i1+i2*n1+i3*n1*n2] =
							FD_D2_O8_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
				}
				else if (fdOrder == 12)
				{
					w[i1+i2*n1+i3*n1*n2] =
							FD_D2_O12_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
				}
				else if (fdOrder == 16)
				{
					w[i1+i2*n1+i3*n1*n2] =
							FD_D2_O16_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
				}
			}

			// compute FD Laplacian for 2D
			else if (dim == DIM2)
			{
				if (fdOrder == 2)
				{
					w[i1+i2*n1+i3*n1*n2] =
							FD_D2_O2_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
							+ FD_D2_O2_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
				}
				else if (fdOrder == 4)
				{
					w[i1+i2*n1+i3*n1*n2] =
							FD_D2_O4_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
							+ FD_D2_O4_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
				}
				else if (fdOrder == 8)
				{
					w[i1+i2*n1+i3*n1*n2] =
							FD_D2_O8_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
							+ FD_D2_O8_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
				}
				else if (fdOrder == 12)
				{
					w[i1+i2*n1+i3*n1*n2] =
							FD_D2_O12_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
							+ FD_D2_O12_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
				}
				else if (fdOrder == 16)
				{
					w[i1+i2*n1+i3*n1*n2] =
							FD_D2_O16_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
							+ FD_D2_O16_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
				}
			}

			// compute FD Laplacian for 3D
			else if (dim == DIM3)
			{
				if (fdOrder == 2)
				{
					w[i1+i2*n1+i3*n1*n2] =
							FD_D2_O2_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
							+ FD_D2_O2_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
							+ FD_D2_O2_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
				}
				else if (fdOrder == 4)
				{
					w[i1+i2*n1+i3*n1*n2] =
							FD_D2_O4_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
							+ FD_D2_O4_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
							+ FD_D2_O4_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
				}
				else if (fdOrder == 8)
				{
					w[i1+i2*n1+i3*n1*n2] =
							FD_D2_O8_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
							+ FD_D2_O8_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
							+ FD_D2_O8_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
				}
				else if (fdOrder == 12)
				{
					w[i1+i2*n1+i3*n1*n2] =
							FD_D2_O12_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
							+ FD_D2_O12_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
							+ FD_D2_O12_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
				}
				else if (fdOrder == 16)
				{
					w[i1+i2*n1+i3*n1*n2] =
							FD_D2_O16_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
							+ FD_D2_O16_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
							+ FD_D2_O16_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
				}
			}
		}

		tid += blockDim.x * gridDim.x;
	}
}

Rtn_code Grid_Hip::FD_LAPLACIAN(Point_type pointType, const Grid& Wgrid, Myint fdOrder)
{
	printDebug(MID_DEBUG, "IN Grid_Hip::FD_LAPLACIAN");

	// check grids are same size
	if (this->sameSize(Wgrid) != true)
	{
		printError("Grid::FD_LAPLACIAN, grids have not same size") ;
		return(RTN_CODE_KO) ;
	}

	Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
	getGridIndex(INNER_POINTS, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End) ;

	const Myfloat inv_d1  = Myfloat(1.0) / d1 ;
	const Myfloat inv_d2  = Myfloat(1.0) / d2 ;
	const Myfloat inv_d3  = Myfloat(1.0) / d3 ;

	const Myfloat inv2_d1 = inv_d1 * inv_d1 ;
	const Myfloat inv2_d2 = inv_d2 * inv_d2 ;
	const Myfloat inv2_d3 = inv_d3 * inv_d3 ;

	Myfloat * d_w = ((Grid_Hip&) Wgrid).d_grid_3d ;
	Myfloat * d_u = this->d_grid_3d ;

	hipLaunchKernelGGL(kernel_FD_LAPLACIAN, dim3(1024), dim3(256), 0, 0, dim, fdOrder, d_w, d_u,inv2_d1,inv2_d2,inv2_d3,
				n1,n2,n3,i1Start,i1End,i2Start,i2End,i3Start,i3End);

	cudaCheckError();
	hipDeviceSynchronize();

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

	const Myfloat inv_d1  = Myfloat(1.0) / d1 ;
	const Myfloat inv_d2  = Myfloat(1.0) / d2 ;
	const Myfloat inv_d3  = Myfloat(1.0) / d3 ;

	const Myfloat inv2_d1 = inv_d1 * inv_d1 ;
	const Myfloat inv2_d2 = inv_d2 * inv_d2 ;
	const Myfloat inv2_d3 = inv_d3 * inv_d3 ;

	Myfloat *prc_d_grid_3d = ((Grid_Hip&) prcGridIn).d_grid_3d ;
	Myfloat *coef_d_grid_3d = ((Grid_Hip&) coefGridIn).d_grid_3d ;
	
	hipLaunchKernelGGL(kernel_computePressureWithFD, dim3(1024), dim3(256), 0, 0, dim, fdOrder, d_grid_3d, prc_d_grid_3d, coef_d_grid_3d,inv2_d1,inv2_d2,inv2_d3,
					n1,n2,n3,i1Start,i1End,i2Start,i2End,i3Start,i3End);

	cudaCheckError();
	hipDeviceSynchronize();
	
	printDebug(FULL_DEBUG, "Out Grid_Hip::computePressureWithFD") ;
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

void Grid_Hip::initializeGrid(void)
{
	printDebug(FULL_DEBUG, "In Grid_Hip::initializeGrid") ;

	Grid::initializeGrid() ;

	if (d_grid_3d == NULL)
	{
		// allocate the grid on the device
		hipMalloc( (void**)&d_grid_3d, npoint * sizeof(Myfloat) );
		cudaCheckError();

		// allocate 1d array of the device used to perform reduction operation
		hipMalloc( (void**)&d_help_3d, (gpuGridSize) * sizeof(Myfloat) );
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

	// this function is critical for validation purpose of most of functions
	// however, it is not included in the performance benchmark
	// it is therefore convenient to build the grid on the CPU
	// and then copy the grid to the GPU

	Grid::fill(pointType, t1, t2, t3, param1, param2, param3, amp) ;
	copyGridHostToDevice(pointType) ;

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

Myfloat Grid_Hip::getMin(Point_type pointType)
{
	printDebug(FULL_DEBUG, "In Grid_Hip::getMin") ;

	Myfloat val = 0 ;

	Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
	getGridIndex(pointType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End);

	hipLaunchKernelGGL(kernel_multiBlk_minval, dim3(gpuGridSize), dim3(gpuBlkSize), gpuBlkSize * sizeof(Myfloat), 0, d_grid_3d, d_help_3d,
			n1, n2, n3, i1Start, i1End, i2Start, i2End, i3Start, i3End) ;
	hipDeviceSynchronize();

	hipLaunchKernelGGL(kernel_singleBlk_minval, dim3(1), dim3(gpuBlkSize), 0, 0, d_help_3d, gpuGridSize) ;
	hipDeviceSynchronize();

	hipMemcpy(&val, &(d_help_3d[0]), sizeof(Myfloat), hipMemcpyDeviceToHost);
	cudaCheckError();

	printDebug(FULL_DEBUG, "Out Grid_Hip::getMin") ;

	return val ;
}

//-------------------------------------------------------------------------------------------------------

Myfloat Grid_Hip::getMax(Point_type pointType)
{
	printDebug(FULL_DEBUG, "In Grid_Hip::getMax") ;

	Myfloat val = 0 ;

	Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
	getGridIndex(pointType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End);

	hipLaunchKernelGGL(kernel_multiBlk_maxval, dim3(gpuGridSize), dim3(gpuBlkSize), gpuBlkSize * sizeof(Myfloat), 0, d_grid_3d, d_help_3d,
			n1, n2, n3, i1Start, i1End, i2Start, i2End, i3Start, i3End) ;
	hipDeviceSynchronize();

	hipLaunchKernelGGL(kernel_singleBlk_maxval, dim3(1), dim3(gpuBlkSize), 0, 0, d_help_3d, gpuGridSize) ;
	hipDeviceSynchronize();

	hipMemcpy(&val, &(d_help_3d[0]), sizeof(Myfloat), hipMemcpyDeviceToHost);
	cudaCheckError();

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

	if (pointType == ALL_POINTS)
	{
		// copy all points with one call to hipMemcpy
		Myint64 idx = 0 ;
		hipMemcpy(&(grid_3d[idx]), &(d_grid_3d[idx]), npoint * sizeof(Myfloat), hipMemcpyDeviceToHost) ;
	}
	else
	{
		for (Myint64 i3 = i3Start; i3<= i3End; i3++)
		{
			for (Myint64 i2 = i2Start; i2<= i2End; i2++)
			{
				// copy 1d segment from i1Start to i1End
				Myint64 i1 = i1Start ;
				Myint64 idx = i1+i2*n1+i3*n1*n2 ;
				Myint64 nn = i1End - i1Start + 1 ;
				hipMemcpy(&(grid_3d[idx]), &(d_grid_3d[idx]), nn * sizeof(Myfloat), hipMemcpyDeviceToHost) ;
			}
		}
	}

	hipDeviceSynchronize();

	printDebug(FULL_DEBUG, "Out Grid_Hip::copyGridDeviceToHost") ;
}

//-------------------------------------------------------------------------------------------------------

void Grid_Hip::copyGridHostToDevice(Point_type pointType)
{
	printDebug(FULL_DEBUG, "In Grid_Hip::copyGridHostToDevice") ;

	Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
	getGridIndex(pointType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End);

	if (pointType == ALL_POINTS)
	{
		// copy all points with one call to hipMemcpy
		Myint64 idx = 0 ;
		hipMemcpy(&(d_grid_3d[idx]), &(grid_3d[idx]), npoint * sizeof(Myfloat), hipMemcpyHostToDevice) ;
	}
	else
	{
		for (Myint64 i3 = i3Start; i3<= i3End; i3++)
		{
			for (Myint64 i2 = i2Start; i2<= i2End; i2++)
			{
				// copy 1d segment
				Myint64 i1 = i1Start ;
				Myint64 idx = i1+i2*n1+i3*n1*n2 ;
				Myint64 nn = i1End - i1Start + 1 ;
				hipMemcpy(&(d_grid_3d[idx]), &(grid_3d[idx]), nn * sizeof(Myfloat), hipMemcpyHostToDevice) ;
			}
		}
	}

	hipDeviceSynchronize();

	printDebug(FULL_DEBUG, "Out Grid_Hip::copyGridHostToDevice") ;
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


	hipLaunchKernelGGL(kernel_applyBoundaryCondition, dim3(1024), dim3(256), 0, 0, dim, d_grid_3d, n1, n2, n3,
			getNeighbourProc(I1HALO1), i1halo1_i1Start, i1halo1_i1End, i1halo1_i2Start, i1halo1_i2End, i1halo1_i3Start, i1halo1_i3End,
			getNeighbourProc(I1HALO2), i1halo2_i1Start, i1halo2_i1End, i1halo2_i2Start, i1halo2_i2End, i1halo2_i3Start, i1halo2_i3End,
			getNeighbourProc(I2HALO1), i2halo1_i1Start, i2halo1_i1End, i2halo1_i2Start, i2halo1_i2End, i2halo1_i3Start, i2halo1_i3End,
			getNeighbourProc(I2HALO2), i2halo2_i1Start, i2halo2_i1End, i2halo2_i2Start, i2halo2_i2End, i2halo2_i3Start, i2halo2_i3End,
			getNeighbourProc(I3HALO1), i3halo1_i1Start, i3halo1_i1End, i3halo1_i2Start, i3halo1_i2End, i3halo1_i3Start, i3halo1_i3End,
			getNeighbourProc(I3HALO2), i3halo2_i1Start, i3halo2_i1End, i3halo2_i2Start, i3halo2_i2End, i3halo2_i3Start, i3halo2_i3End);

	hipDeviceSynchronize();

	printDebug(FULL_DEBUG, "Out Grid_Hip::applyBoundaryCondition") ;
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Myfloat Grid_Hip::getSumAbs(Point_type pointType) const
{
	printDebug(LIGHT_DEBUG, "IN Grid_Hip::getSumAbs");

	Myfloat sum = 0 ;

	Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
	getGridIndex(pointType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End);

	hipLaunchKernelGGL(kernel_multiBlk_sumAbs, dim3(gpuGridSize), dim3(gpuBlkSize), gpuBlkSize * sizeof(Myfloat), 0, d_grid_3d, d_help_3d,
			n1, n2, n3, i1Start, i1End, i2Start, i2End, i3Start, i3End) ;
	hipDeviceSynchronize();

	hipLaunchKernelGGL(kernel_singleBlk_sum, dim3(1), dim3(gpuBlkSize), 0, 0, d_help_3d, gpuGridSize) ;
	hipDeviceSynchronize();

	hipMemcpy(&sum, &(d_help_3d[0]), sizeof(Myfloat), hipMemcpyDeviceToHost);
	cudaCheckError();

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

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_Hip::sendWithMPI(Myint64 nGridPoint, Myint procDestId)
{

	printDebug(FULL_DEBUG, "In Grid_Hip::sendWithMPI") ;

	// copy from device to host
	Myint64 idx = 0 ;
	hipMemcpy(&(grid_3d[idx]), &(d_grid_3d[idx]), nGridPoint * sizeof(Myfloat), hipMemcpyDeviceToHost) ;

	MPI_Send(grid_3d, nGridPoint, MPI_MYFLOAT, procDestId, 0, MPI_COMM_WORLD) ;

	printDebug(FULL_DEBUG, "Out Grid_Hip::sendWithMPI") ;
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_Hip::recvWithMPI(Myint64 nGridPoint, Myint procSrcId)
{

	printDebug(FULL_DEBUG, "In Grid_Hip::recvWithMPI") ;

	MPI_Status status ;
	MPI_Recv(grid_3d, nGridPoint, MPI_MYFLOAT, procSrcId, 0, MPI_COMM_WORLD, &status) ;

	// copy from host to device
	Myint64 idx = 0 ;
	hipMemcpy(&(d_grid_3d[idx]), &(grid_3d[idx]), npoint * sizeof(Myfloat), hipMemcpyHostToDevice) ;

	if (status.MPI_ERROR != MPI_SUCCESS)
	{
		//printError("MPI ERROR", status.MPI_ERROR) ;
	}

	printDebug(FULL_DEBUG, "Out Grid_Hip::recvWithMPI") ;
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_Hip::sendRecvWithMPI(const Grid& gridDest, Myint idSend, Myint idRecv, Myint64 nGridPoint)
{

	printDebug(FULL_DEBUG, "In Grid_Hip::sendRecvWithMPI") ;

	Myfloat *bufSend = grid_3d ;
	Myfloat *bufRecv = gridDest.grid_3d ;

	// copy from device to host
	Myint64 idx = 0 ;
	hipMemcpy(&(grid_3d[idx]), &(d_grid_3d[idx]), nGridPoint * sizeof(Myfloat), hipMemcpyDeviceToHost) ;

	MPI_Status status ;
	MPI_Sendrecv(bufSend, nGridPoint, MPI_MYFLOAT, idSend, 0,
			bufRecv, nGridPoint, MPI_MYFLOAT, idRecv, 0,
			MPI_COMM_WORLD, &status) ;
	if (status.MPI_ERROR != MPI_SUCCESS)
	{
		//printError("MPI ERROR", status.MPI_ERROR) ;
	}

	// copy from host to device
	Myfloat *gridDest_d_grid_3d = ((Grid_Hip&) gridDest).d_grid_3d ;
	hipMemcpy(&(gridDest_d_grid_3d[idx]), &(gridDest.grid_3d[idx]), npoint * sizeof(Myfloat), hipMemcpyHostToDevice) ;

	printDebug(FULL_DEBUG, "Out Grid_Hip::sendRecvWithMPI") ;
	return(RTN_CODE_OK) ;
}


} // namespace hpcscan
