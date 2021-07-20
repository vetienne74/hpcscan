
//-------------------------------------------------------------------------------------------------------
// This grid is activated with command line option -testMode CUDA_Optim
// Derived class from Grid_Cuda
// CUDA implementation (target GPU)
//-------------------------------------------------------------------------------------------------------

#include "grid_Cuda_Optim.h"

#include <algorithm> // for min and max
#include <cassert>
#include <cfloat>  // for FLT_MAX ;
#include <cmath>   // for fabs
#include <cstddef> // for NULL
#include <fstream>
#include <stdio.h>

#include <cooperative_groups.h>
#include "mpi.h"

#include "config.h"
#include "constant.h"
#include "fdm.h"
#include "global.h"
#include "output_report.h"

using namespace std;

namespace cg = cooperative_groups;

namespace hpcscan {

#define MAX_BLOCK_DIM_X 64
#define MAX_BLOCK_DIM_Y 16

//-------------------------------------------------------------------------------------------------------

// Macro for checking cuda errors following a cuda launch or api call
#define cudaCheckError() { \
		cudaError_t e=cudaGetLastError(); \
		if(e!=cudaSuccess) { \
			printf("Cuda failure %s:%d: '%s'\n",__FILE__,__LINE__,cudaGetErrorString(e)); \
			printError(" CUDA ERROR") ; \
			exit(0); \
		} \
}

//*******************************************************************************************************
// BEGINING OF CUDA KERNELS
//*******************************************************************************************************

// With max FD order = 16, there are at most 9 FD coefficients in the stencil
#define MAX_FD_COEF 9
__constant__ Myfloat stencil[MAX_FD_COEF];

// Set maximum number of threads per block for all kernels
//__launch_bounds__(MAX_BLOCK_DIM_Y*MAX_BLOCK_DIM_X)

//-------------------------------------------------------------------------------------------------------
// compute Laplacian
//   u
// output w

__global__ void kernelOpt_FD_LAPLACIAN_O4(const Myint fdOrder, Myfloat *output, Myfloat *input,
		const Myfloat inv2_d1, const Myfloat inv2_d2, const Myfloat inv2_d3,
		const Myint n1, const Myint n2, const Myint n3,
		const Myint64 i1Start, const Myint64 i1End, const Myint64 i2Start, const Myint64 i2End, const Myint64 i3Start, const Myint64 i3End,
		const int dimx, const int dimy, const int dimz)
{
#define RADIUS_O4 2
	bool validr = true;
	bool validw = true;

	const int gtidx = blockIdx.x * blockDim.x + threadIdx.x;
	const int gtidy = blockIdx.y * blockDim.y + threadIdx.y;
	const int gtidz = blockIdx.z * dimz ;
	const int ltidx = threadIdx.x;
	const int ltidy = threadIdx.y;
	const int workx = blockDim.x;
	const int worky = blockDim.y;

	// Handle to thread block group
	cg::thread_block cta = cg::this_thread_block();
	__shared__ float tile[MAX_BLOCK_DIM_Y + 2 * RADIUS_O4][MAX_BLOCK_DIM_X + 2 * RADIUS_O4];

	const int stride_y = n1 ;
	const int stride_z = stride_y * n2 ;

	int inputIndex  = gtidz * stride_z ;
	int outputIndex = gtidz * stride_z ;

	// Advance inputIndex to start of inner volume
	inputIndex += i2Start * stride_y + i1Start ;

	// Advance inputIndex to target element
	inputIndex += gtidy * stride_y + gtidx;

	Myfloat infront[RADIUS_O4];
	Myfloat behind[RADIUS_O4];
	Myfloat current;
	const Myfloat fdCoef0 = stencil[0] * (inv2_d1 + inv2_d2 + inv2_d3) ;

	const int tx = ltidx + RADIUS_O4;
	const int ty = ltidy + RADIUS_O4;

	// Check in bounds
	if ((gtidx > i1End) || (gtidy > i2End))
		validr = false;

	if ((gtidx >= dimx) || (gtidy >= dimy))
		validw = false;

	// Preload the "infront" and "behind" data
	for (int i = RADIUS_O4 - 2 ; i >= 0 ; i--)
	{
		if (validr)
			behind[i] = input[inputIndex];

		inputIndex += stride_z;
	}

	if (validr)
		current = input[inputIndex];

	outputIndex = inputIndex;
	inputIndex += stride_z;

	for (int i = 0 ; i < RADIUS_O4 ; i++)
	{
		if (validr)
			infront[i] = input[inputIndex];

		inputIndex += stride_z;
	}

	// set max element to process along z (n3)
	int maxZ = gtidz + dimz ;
	if (maxZ > (i3End-RADIUS_O4+1))
	{
		maxZ = i3End-RADIUS_O4+1 ;
	}

	// Step through the xy-planes
#pragma unroll 5
	for (int iz = gtidz ; iz < maxZ ; iz++)
	{
		// Advance the slice (move the thread-front)
		for (int i = RADIUS_O4 - 1 ; i > 0 ; i--)
			behind[i] = behind[i - 1];

		behind[0] = current;
		current = infront[0];
#pragma unroll 2

		for (int i = 0 ; i < RADIUS_O4 - 1 ; i++)
			infront[i] = infront[i + 1];

		if (validr)
			infront[RADIUS_O4 - 1] = input[inputIndex];

		inputIndex  += stride_z;
		outputIndex += stride_z;
		cg::sync(cta);

		// Note that for the work items on the boundary of the problem, the
		// supplied index when reading the halo (below) may wrap to the
		// previous/next row or even the previous/next xy-plane. This is
		// acceptable since a) we disable the output write for these work
		// items and b) there is at least one xy-plane before/after the
		// current plane, so the access will be within bounds.

		// Update the data slice in the local tile
		// Halo above & below
		if (ltidy < RADIUS_O4)
		{
			tile[ltidy][tx]                  = input[outputIndex - RADIUS_O4 * stride_y];
			tile[ltidy + worky + RADIUS_O4][tx] = input[outputIndex + worky * stride_y];
		}

		// Halo left & right
		if (ltidx < RADIUS_O4)
		{
			tile[ty][ltidx]                  = input[outputIndex - RADIUS_O4];
			tile[ty][ltidx + workx + RADIUS_O4] = input[outputIndex + workx];
		}

		tile[ty][tx] = current;
		cg::sync(cta);

		// Compute the output value
		Myfloat value = fdCoef0 * current ;
#pragma unroll 2

		for (int i = 1 ; i <= RADIUS_O4 ; i++)
		{
			value += stencil[i] * (inv2_d3 * (infront[i-1] + behind[i-1]) +
					inv2_d2 * (tile[ty - i][tx] + tile[ty + i][tx]) +
					inv2_d1 * (tile[ty][tx - i] + tile[ty][tx + i]));
		}

		// Store the output value
		if (validw)
			output[outputIndex] = value;
	}
}

__global__ void kernelOpt_FD_LAPLACIAN_O6(const Myint fdOrder, Myfloat *output, Myfloat *input,
		const Myfloat inv2_d1, const Myfloat inv2_d2, const Myfloat inv2_d3,
		const Myint n1, const Myint n2, const Myint n3,
		const Myint64 i1Start, const Myint64 i1End, const Myint64 i2Start, const Myint64 i2End, const Myint64 i3Start, const Myint64 i3End,
		const int dimx, const int dimy, const int dimz)
{
#define RADIUS_O6 3
	bool validr = true;
	bool validw = true;

	const int gtidx = blockIdx.x * blockDim.x + threadIdx.x;
	const int gtidy = blockIdx.y * blockDim.y + threadIdx.y;
	const int gtidz = blockIdx.z * dimz ;
	const int ltidx = threadIdx.x;
	const int ltidy = threadIdx.y;
	const int workx = blockDim.x;
	const int worky = blockDim.y;

	// Handle to thread block group
	cg::thread_block cta = cg::this_thread_block();
	__shared__ float tile[MAX_BLOCK_DIM_Y + 2 * RADIUS_O6][MAX_BLOCK_DIM_X + 2 * RADIUS_O6];

	const int stride_y = n1 ;
	const int stride_z = stride_y * n2 ;

	int inputIndex  = gtidz * stride_z ;
	int outputIndex = gtidz * stride_z ;

	// Advance inputIndex to start of inner volume
	inputIndex += i2Start * stride_y + i1Start ;

	// Advance inputIndex to target element
	inputIndex += gtidy * stride_y + gtidx;

	Myfloat infront[RADIUS_O6];
	Myfloat behind[RADIUS_O6];
	Myfloat current;
	const Myfloat fdCoef0 = stencil[0] * (inv2_d1 + inv2_d2 + inv2_d3) ;

	const int tx = ltidx + RADIUS_O6;
	const int ty = ltidy + RADIUS_O6;

	// Check in bounds
	if ((gtidx > i1End) || (gtidy > i2End))
		validr = false;

	if ((gtidx >= dimx) || (gtidy >= dimy))
		validw = false;

	// Preload the "infront" and "behind" data
	for (int i = RADIUS_O6 - 2 ; i >= 0 ; i--)
	{
		if (validr)
			behind[i] = input[inputIndex];

		inputIndex += stride_z;
	}

	if (validr)
		current = input[inputIndex];

	outputIndex = inputIndex;
	inputIndex += stride_z;

	for (int i = 0 ; i < RADIUS_O6 ; i++)
	{
		if (validr)
			infront[i] = input[inputIndex];

		inputIndex += stride_z;
	}

	// set max element to process along z (n3)
	int maxZ = gtidz + dimz ;
	if (maxZ > (i3End-RADIUS_O6+1))
	{
		maxZ = i3End-RADIUS_O6+1 ;
	}

	// Step through the xy-planes
#pragma unroll 7
	for (int iz = gtidz ; iz < maxZ ; iz++)
	{
		// Advance the slice (move the thread-front)
		for (int i = RADIUS_O6 - 1 ; i > 0 ; i--)
			behind[i] = behind[i - 1];

		behind[0] = current;
		current = infront[0];
#pragma unroll 3

		for (int i = 0 ; i < RADIUS_O6 - 1 ; i++)
			infront[i] = infront[i + 1];

		if (validr)
			infront[RADIUS_O6 - 1] = input[inputIndex];

		inputIndex  += stride_z;
		outputIndex += stride_z;
		cg::sync(cta);

		// Note that for the work items on the boundary of the problem, the
		// supplied index when reading the halo (below) may wrap to the
		// previous/next row or even the previous/next xy-plane. This is
		// acceptable since a) we disable the output write for these work
		// items and b) there is at least one xy-plane before/after the
		// current plane, so the access will be within bounds.

		// Update the data slice in the local tile
		// Halo above & below
		if (ltidy < RADIUS_O6)
		{
			tile[ltidy][tx]                  = input[outputIndex - RADIUS_O6 * stride_y];
			tile[ltidy + worky + RADIUS_O6][tx] = input[outputIndex + worky * stride_y];
		}

		// Halo left & right
		if (ltidx < RADIUS_O6)
		{
			tile[ty][ltidx]                  = input[outputIndex - RADIUS_O6];
			tile[ty][ltidx + workx + RADIUS_O6] = input[outputIndex + workx];
		}

		tile[ty][tx] = current;
		cg::sync(cta);

		// Compute the output value
		Myfloat value = fdCoef0 * current ;
#pragma unroll 3

		for (int i = 1 ; i <= RADIUS_O6 ; i++)
		{
			value += stencil[i] * (inv2_d3 * (infront[i-1] + behind[i-1]) +
					inv2_d2 * (tile[ty - i][tx] + tile[ty + i][tx]) +
					inv2_d1 * (tile[ty][tx - i] + tile[ty][tx + i]));
		}

		// Store the output value
		if (validw)
			output[outputIndex] = value;
	}
}

__global__ void kernelOpt_FD_LAPLACIAN_O8(const Myint fdOrder, Myfloat *output, Myfloat *input,
		const Myfloat inv2_d1, const Myfloat inv2_d2, const Myfloat inv2_d3,
		const Myint n1, const Myint n2, const Myint n3,
		const Myint64 i1Start, const Myint64 i1End, const Myint64 i2Start, const Myint64 i2End, const Myint64 i3Start, const Myint64 i3End,
		const int dimx, const int dimy, const int dimz)
{
#define RADIUS_O8 4
	bool validr = true;
	bool validw = true;

	const int gtidx = blockIdx.x * blockDim.x + threadIdx.x;
	const int gtidy = blockIdx.y * blockDim.y + threadIdx.y;
	const int gtidz = blockIdx.z * dimz ;
	const int ltidx = threadIdx.x;
	const int ltidy = threadIdx.y;
	const int workx = blockDim.x;
	const int worky = blockDim.y;

	// Handle to thread block group
	cg::thread_block cta = cg::this_thread_block();
	__shared__ float tile[MAX_BLOCK_DIM_Y + 2 * RADIUS_O8][MAX_BLOCK_DIM_X + 2 * RADIUS_O8];

	const int stride_y = n1 ;
	const int stride_z = stride_y * n2 ;

	int inputIndex  = gtidz * stride_z ;
	int outputIndex = gtidz * stride_z ;

	// Advance inputIndex to start of inner volume
	inputIndex += i2Start * stride_y + i1Start ;

	// Advance inputIndex to target element
	inputIndex += gtidy * stride_y + gtidx;

	Myfloat infront[RADIUS_O8];
	Myfloat behind[RADIUS_O8];
	Myfloat current;
	const Myfloat fdCoef0 = stencil[0] * (inv2_d1 + inv2_d2 + inv2_d3) ;

	const int tx = ltidx + RADIUS_O8;
	const int ty = ltidy + RADIUS_O8;

	// Check in bounds
	if ((gtidx > i1End) || (gtidy > i2End))
		validr = false;

	if ((gtidx >= dimx) || (gtidy >= dimy))
		validw = false;

	// Preload the "infront" and "behind" data
	for (int i = RADIUS_O8 - 2 ; i >= 0 ; i--)
	{
		if (validr)
			behind[i] = input[inputIndex];

		inputIndex += stride_z;
	}

	if (validr)
		current = input[inputIndex];

	outputIndex = inputIndex;
	inputIndex += stride_z;

	for (int i = 0 ; i < RADIUS_O8 ; i++)
	{
		if (validr)
			infront[i] = input[inputIndex];

		inputIndex += stride_z;
	}

	// set max element to process along z (n3)
	int maxZ = gtidz + dimz ;
	if (maxZ > (i3End-RADIUS_O8+1))
	{
		maxZ = i3End-RADIUS_O8+1 ;
	}

	// Step through the xy-planes
#pragma unroll 9	
	for (int iz = gtidz ; iz < maxZ ; iz++)
	{			
		// Advance the slice (move the thread-front)
		for (int i = RADIUS_O8 - 1 ; i > 0 ; i--)
			behind[i] = behind[i - 1];

		behind[0] = current;
		current = infront[0];
#pragma unroll 4

		for (int i = 0 ; i < RADIUS_O8 - 1 ; i++)
			infront[i] = infront[i + 1];

		if (validr)
			infront[RADIUS_O8 - 1] = input[inputIndex];

		inputIndex  += stride_z;
		outputIndex += stride_z;
		cg::sync(cta);

		// Note that for the work items on the boundary of the problem, the
		// supplied index when reading the halo (below) may wrap to the
		// previous/next row or even the previous/next xy-plane. This is
		// acceptable since a) we disable the output write for these work
		// items and b) there is at least one xy-plane before/after the
		// current plane, so the access will be within bounds.

		// Update the data slice in the local tile
		// Halo above & below
		if (ltidy < RADIUS_O8)
		{
			tile[ltidy][tx]                  = input[outputIndex - RADIUS_O8 * stride_y];
			tile[ltidy + worky + RADIUS_O8][tx] = input[outputIndex + worky * stride_y];
		}

		// Halo left & right
		if (ltidx < RADIUS_O8)
		{
			tile[ty][ltidx]                  = input[outputIndex - RADIUS_O8];
			tile[ty][ltidx + workx + RADIUS_O8] = input[outputIndex + workx];
		}

		tile[ty][tx] = current;
		cg::sync(cta);

		// Compute the output value
		Myfloat value = fdCoef0 * current ;
#pragma unroll 4

		for (int i = 1 ; i <= RADIUS_O8 ; i++)
		{
			value += stencil[i] * (inv2_d3 * (infront[i-1] + behind[i-1]) + 
					inv2_d2 * (tile[ty - i][tx] + tile[ty + i][tx]) +
					inv2_d1 * (tile[ty][tx - i] + tile[ty][tx + i]));
		}

		// Store the output value
		if (validw)
			output[outputIndex] = value;
	}
}

__global__ void kernelOpt_FD_LAPLACIAN_O10(const Myint fdOrder, Myfloat *output, Myfloat *input,
		const Myfloat inv2_d1, const Myfloat inv2_d2, const Myfloat inv2_d3,
		const Myint n1, const Myint n2, const Myint n3,
		const Myint64 i1Start, const Myint64 i1End, const Myint64 i2Start, const Myint64 i2End, const Myint64 i3Start, const Myint64 i3End,
		const int dimx, const int dimy, const int dimz)
{
#define RADIUS_O10 5
	bool validr = true;
	bool validw = true;

	const int gtidx = blockIdx.x * blockDim.x + threadIdx.x;
	const int gtidy = blockIdx.y * blockDim.y + threadIdx.y;
	const int gtidz = blockIdx.z * dimz ;
	const int ltidx = threadIdx.x;
	const int ltidy = threadIdx.y;
	const int workx = blockDim.x;
	const int worky = blockDim.y;

	// Handle to thread block group
	cg::thread_block cta = cg::this_thread_block();
	__shared__ float tile[MAX_BLOCK_DIM_Y + 2 * RADIUS_O10][MAX_BLOCK_DIM_X + 2 * RADIUS_O10];

	const int stride_y = n1 ;
	const int stride_z = stride_y * n2 ;

	int inputIndex  = gtidz * stride_z ;
	int outputIndex = gtidz * stride_z ;

	// Advance inputIndex to start of inner volume
	inputIndex += i2Start * stride_y + i1Start ;

	// Advance inputIndex to target element
	inputIndex += gtidy * stride_y + gtidx;

	Myfloat infront[RADIUS_O10];
	Myfloat behind[RADIUS_O10];
	Myfloat current;
	const Myfloat fdCoef0 = stencil[0] * (inv2_d1 + inv2_d2 + inv2_d3) ;

	const int tx = ltidx + RADIUS_O10;
	const int ty = ltidy + RADIUS_O10;

	// Check in bounds
	if ((gtidx > i1End) || (gtidy > i2End))
		validr = false;

	if ((gtidx >= dimx) || (gtidy >= dimy))
		validw = false;

	// Preload the "infront" and "behind" data
	for (int i = RADIUS_O10 - 2 ; i >= 0 ; i--)
	{
		if (validr)
			behind[i] = input[inputIndex];

		inputIndex += stride_z;
	}

	if (validr)
		current = input[inputIndex];

	outputIndex = inputIndex;
	inputIndex += stride_z;

	for (int i = 0 ; i < RADIUS_O10 ; i++)
	{
		if (validr)
			infront[i] = input[inputIndex];

		inputIndex += stride_z;
	}

	// set max element to process along z (n3)
	int maxZ = gtidz + dimz ;
	if (maxZ > (i3End-RADIUS_O10+1))
	{
		maxZ = i3End-RADIUS_O10+1 ;
	}

	// Step through the xy-planes
#pragma unroll 11
	for (int iz = gtidz ; iz < maxZ ; iz++)
	{
		// Advance the slice (move the thread-front)
		for (int i = RADIUS_O10 - 1 ; i > 0 ; i--)
			behind[i] = behind[i - 1];

		behind[0] = current;
		current = infront[0];
#pragma unroll 5

		for (int i = 0 ; i < RADIUS_O10 - 1 ; i++)
			infront[i] = infront[i + 1];

		if (validr)
			infront[RADIUS_O10 - 1] = input[inputIndex];

		inputIndex  += stride_z;
		outputIndex += stride_z;
		cg::sync(cta);

		// Note that for the work items on the boundary of the problem, the
		// supplied index when reading the halo (below) may wrap to the
		// previous/next row or even the previous/next xy-plane. This is
		// acceptable since a) we disable the output write for these work
		// items and b) there is at least one xy-plane before/after the
		// current plane, so the access will be within bounds.

		// Update the data slice in the local tile
		// Halo above & below
		if (ltidy < RADIUS_O10)
		{
			tile[ltidy][tx]                  = input[outputIndex - RADIUS_O10 * stride_y];
			tile[ltidy + worky + RADIUS_O10][tx] = input[outputIndex + worky * stride_y];
		}

		// Halo left & right
		if (ltidx < RADIUS_O10)
		{
			tile[ty][ltidx]                  = input[outputIndex - RADIUS_O10];
			tile[ty][ltidx + workx + RADIUS_O10] = input[outputIndex + workx];
		}

		tile[ty][tx] = current;
		cg::sync(cta);

		// Compute the output value
		Myfloat value = fdCoef0 * current ;
#pragma unroll 5

		for (int i = 1 ; i <= RADIUS_O10 ; i++)
		{
			value += stencil[i] * (inv2_d3 * (infront[i-1] + behind[i-1]) +
					inv2_d2 * (tile[ty - i][tx] + tile[ty + i][tx]) +
					inv2_d1 * (tile[ty][tx - i] + tile[ty][tx + i]));
		}

		// Store the output value
		if (validw)
			output[outputIndex] = value;
	}
}

__global__ void kernelOpt_FD_LAPLACIAN_O12(const Myint fdOrder, Myfloat *output, Myfloat *input,
		const Myfloat inv2_d1, const Myfloat inv2_d2, const Myfloat inv2_d3,
		const Myint n1, const Myint n2, const Myint n3,
		const Myint64 i1Start, const Myint64 i1End, const Myint64 i2Start, const Myint64 i2End, const Myint64 i3Start, const Myint64 i3End,
		const int dimx, const int dimy, const int dimz)
{
#define RADIUS_O12 6
	bool validr = true;
	bool validw = true;

	const int gtidx = blockIdx.x * blockDim.x + threadIdx.x;
	const int gtidy = blockIdx.y * blockDim.y + threadIdx.y;
	const int gtidz = blockIdx.z * dimz ;
	const int ltidx = threadIdx.x;
	const int ltidy = threadIdx.y;
	const int workx = blockDim.x;
	const int worky = blockDim.y;

	// Handle to thread block group
	cg::thread_block cta = cg::this_thread_block();
	__shared__ float tile[MAX_BLOCK_DIM_Y + 2 * RADIUS_O12][MAX_BLOCK_DIM_X + 2 * RADIUS_O12];

	const int stride_y = n1 ;
	const int stride_z = stride_y * n2 ;

	int inputIndex  = gtidz * stride_z ;
	int outputIndex = gtidz * stride_z ;

	// Advance inputIndex to start of inner volume
	inputIndex += i2Start * stride_y + i1Start ;

	// Advance inputIndex to target element
	inputIndex += gtidy * stride_y + gtidx;

	Myfloat infront[RADIUS_O12];
	Myfloat behind[RADIUS_O12];
	Myfloat current;
	const Myfloat fdCoef0 = stencil[0] * (inv2_d1 + inv2_d2 + inv2_d3) ;

	const int tx = ltidx + RADIUS_O12;
	const int ty = ltidy + RADIUS_O12;

	// Check in bounds
	if ((gtidx > i1End) || (gtidy > i2End))
		validr = false;

	if ((gtidx >= dimx) || (gtidy >= dimy))
		validw = false;

	// Preload the "infront" and "behind" data
	for (int i = RADIUS_O12 - 2 ; i >= 0 ; i--)
	{
		if (validr)
			behind[i] = input[inputIndex];

		inputIndex += stride_z;
	}

	if (validr)
		current = input[inputIndex];

	outputIndex = inputIndex;
	inputIndex += stride_z;

	for (int i = 0 ; i < RADIUS_O12 ; i++)
	{
		if (validr)
			infront[i] = input[inputIndex];

		inputIndex += stride_z;
	}

	// set max element to process along z (n3)
	int maxZ = gtidz + dimz ;
	if (maxZ > (i3End-RADIUS_O12+1))
	{
		maxZ = i3End-RADIUS_O12+1 ;
	}

	// Step through the xy-planes
#pragma unroll 13
	for (int iz = gtidz ; iz < maxZ ; iz++)
	{
		// Advance the slice (move the thread-front)
		for (int i = RADIUS_O12 - 1 ; i > 0 ; i--)
			behind[i] = behind[i - 1];

		behind[0] = current;
		current = infront[0];
#pragma unroll 6

		for (int i = 0 ; i < RADIUS_O12 - 1 ; i++)
			infront[i] = infront[i + 1];

		if (validr)
			infront[RADIUS_O12 - 1] = input[inputIndex];

		inputIndex  += stride_z;
		outputIndex += stride_z;
		cg::sync(cta);

		// Note that for the work items on the boundary of the problem, the
		// supplied index when reading the halo (below) may wrap to the
		// previous/next row or even the previous/next xy-plane. This is
		// acceptable since a) we disable the output write for these work
		// items and b) there is at least one xy-plane before/after the
		// current plane, so the access will be within bounds.

		// Update the data slice in the local tile
		// Halo above & below
		if (ltidy < RADIUS_O12)
		{
			tile[ltidy][tx]                  = input[outputIndex - RADIUS_O12 * stride_y];
			tile[ltidy + worky + RADIUS_O12][tx] = input[outputIndex + worky * stride_y];
		}

		// Halo left & right
		if (ltidx < RADIUS_O12)
		{
			tile[ty][ltidx]                  = input[outputIndex - RADIUS_O12];
			tile[ty][ltidx + workx + RADIUS_O12] = input[outputIndex + workx];
		}

		tile[ty][tx] = current;
		cg::sync(cta);

		// Compute the output value
		Myfloat value = fdCoef0 * current ;
#pragma unroll 6

		for (int i = 1 ; i <= RADIUS_O12 ; i++)
		{
			value += stencil[i] * (inv2_d3 * (infront[i-1] + behind[i-1]) +
					inv2_d2 * (tile[ty - i][tx] + tile[ty + i][tx]) +
					inv2_d1 * (tile[ty][tx - i] + tile[ty][tx + i]));
		}

		// Store the output value
		if (validw)
			output[outputIndex] = value;
	}
}

__global__ void kernelOpt_FD_LAPLACIAN_O14(const Myint fdOrder, Myfloat *output, Myfloat *input,
		const Myfloat inv2_d1, const Myfloat inv2_d2, const Myfloat inv2_d3,
		const Myint n1, const Myint n2, const Myint n3,
		const Myint64 i1Start, const Myint64 i1End, const Myint64 i2Start, const Myint64 i2End, const Myint64 i3Start, const Myint64 i3End,
		const int dimx, const int dimy, const int dimz)
{
#define RADIUS_O14 7
	bool validr = true;
	bool validw = true;

	const int gtidx = blockIdx.x * blockDim.x + threadIdx.x;
	const int gtidy = blockIdx.y * blockDim.y + threadIdx.y;
	const int gtidz = blockIdx.z * dimz ;
	const int ltidx = threadIdx.x;
	const int ltidy = threadIdx.y;
	const int workx = blockDim.x;
	const int worky = blockDim.y;

	// Handle to thread block group
	cg::thread_block cta = cg::this_thread_block();
	__shared__ float tile[MAX_BLOCK_DIM_Y + 2 * RADIUS_O14][MAX_BLOCK_DIM_X + 2 * RADIUS_O14];

	const int stride_y = n1 ;
	const int stride_z = stride_y * n2 ;

	int inputIndex  = gtidz * stride_z ;
	int outputIndex = gtidz * stride_z ;

	// Advance inputIndex to start of inner volume
	inputIndex += i2Start * stride_y + i1Start ;

	// Advance inputIndex to target element
	inputIndex += gtidy * stride_y + gtidx;

	Myfloat infront[RADIUS_O14];
	Myfloat behind[RADIUS_O14];
	Myfloat current;
	const Myfloat fdCoef0 = stencil[0] * (inv2_d1 + inv2_d2 + inv2_d3) ;

	const int tx = ltidx + RADIUS_O14;
	const int ty = ltidy + RADIUS_O14;

	// Check in bounds
	if ((gtidx > i1End) || (gtidy > i2End))
		validr = false;

	if ((gtidx >= dimx) || (gtidy >= dimy))
		validw = false;

	// Preload the "infront" and "behind" data
	for (int i = RADIUS_O14 - 2 ; i >= 0 ; i--)
	{
		if (validr)
			behind[i] = input[inputIndex];

		inputIndex += stride_z;
	}

	if (validr)
		current = input[inputIndex];

	outputIndex = inputIndex;
	inputIndex += stride_z;

	for (int i = 0 ; i < RADIUS_O14 ; i++)
	{
		if (validr)
			infront[i] = input[inputIndex];

		inputIndex += stride_z;
	}

	// set max element to process along z (n3)
	int maxZ = gtidz + dimz ;
	if (maxZ > (i3End-RADIUS_O14+1))
	{
		maxZ = i3End-RADIUS_O14+1 ;
	}

	// Step through the xy-planes
#pragma unroll 15
	for (int iz = gtidz ; iz < maxZ ; iz++)
	{
		// Advance the slice (move the thread-front)
		for (int i = RADIUS_O14 - 1 ; i > 0 ; i--)
			behind[i] = behind[i - 1];

		behind[0] = current;
		current = infront[0];
#pragma unroll 7

		for (int i = 0 ; i < RADIUS_O14 - 1 ; i++)
			infront[i] = infront[i + 1];

		if (validr)
			infront[RADIUS_O14 - 1] = input[inputIndex];

		inputIndex  += stride_z;
		outputIndex += stride_z;
		cg::sync(cta);

		// Note that for the work items on the boundary of the problem, the
		// supplied index when reading the halo (below) may wrap to the
		// previous/next row or even the previous/next xy-plane. This is
		// acceptable since a) we disable the output write for these work
		// items and b) there is at least one xy-plane before/after the
		// current plane, so the access will be within bounds.

		// Update the data slice in the local tile
		// Halo above & below
		if (ltidy < RADIUS_O14)
		{
			tile[ltidy][tx]                  = input[outputIndex - RADIUS_O14 * stride_y];
			tile[ltidy + worky + RADIUS_O14][tx] = input[outputIndex + worky * stride_y];
		}

		// Halo left & right
		if (ltidx < RADIUS_O14)
		{
			tile[ty][ltidx]                  = input[outputIndex - RADIUS_O14];
			tile[ty][ltidx + workx + RADIUS_O14] = input[outputIndex + workx];
		}

		tile[ty][tx] = current;
		cg::sync(cta);

		// Compute the output value
		Myfloat value = fdCoef0 * current ;
#pragma unroll 7

		for (int i = 1 ; i <= RADIUS_O14 ; i++)
		{
			value += stencil[i] * (inv2_d3 * (infront[i-1] + behind[i-1]) +
					inv2_d2 * (tile[ty - i][tx] + tile[ty + i][tx]) +
					inv2_d1 * (tile[ty][tx - i] + tile[ty][tx + i]));
		}

		// Store the output value
		if (validw)
			output[outputIndex] = value;
	}
}

__global__ void kernelOpt_FD_LAPLACIAN_O16(const Myint fdOrder, Myfloat *output, Myfloat *input,
		const Myfloat inv2_d1, const Myfloat inv2_d2, const Myfloat inv2_d3,
		const Myint n1, const Myint n2, const Myint n3,
		const Myint64 i1Start, const Myint64 i1End, const Myint64 i2Start, const Myint64 i2End, const Myint64 i3Start, const Myint64 i3End,
		const int dimx, const int dimy, const int dimz)
{
#define RADIUS_O16 8
	bool validr = true;
	bool validw = true;

	const int gtidx = blockIdx.x * blockDim.x + threadIdx.x;
	const int gtidy = blockIdx.y * blockDim.y + threadIdx.y;
	const int gtidz = blockIdx.z * dimz ;
	const int ltidx = threadIdx.x;
	const int ltidy = threadIdx.y;
	const int workx = blockDim.x;
	const int worky = blockDim.y;

	// Handle to thread block group
	cg::thread_block cta = cg::this_thread_block();
	__shared__ float tile[MAX_BLOCK_DIM_Y + 2 * RADIUS_O16][MAX_BLOCK_DIM_X + 2 * RADIUS_O16];

	const int stride_y = n1 ;
	const int stride_z = stride_y * n2 ;

	int inputIndex  = gtidz * stride_z ;
	int outputIndex = gtidz * stride_z ;

	// Advance inputIndex to start of inner volume
	inputIndex += i2Start * stride_y + i1Start ;

	// Advance inputIndex to target element
	inputIndex += gtidy * stride_y + gtidx;

	Myfloat infront[RADIUS_O16];
	Myfloat behind[RADIUS_O16];
	Myfloat current;
	const Myfloat fdCoef0 = stencil[0] * (inv2_d1 + inv2_d2 + inv2_d3) ;

	const int tx = ltidx + RADIUS_O16;
	const int ty = ltidy + RADIUS_O16;

	// Check in bounds
	if ((gtidx > i1End) || (gtidy > i2End))
		validr = false;

	if ((gtidx >= dimx) || (gtidy >= dimy))
		validw = false;

	// Preload the "infront" and "behind" data
	for (int i = RADIUS_O16 - 2 ; i >= 0 ; i--)
	{
		if (validr)
			behind[i] = input[inputIndex];

		inputIndex += stride_z;
	}

	if (validr)
		current = input[inputIndex];

	outputIndex = inputIndex;
	inputIndex += stride_z;

	for (int i = 0 ; i < RADIUS_O16 ; i++)
	{
		if (validr)
			infront[i] = input[inputIndex];

		inputIndex += stride_z;
	}

	// set max element to process along z (n3)
	int maxZ = gtidz + dimz ;
	if (maxZ > (i3End-RADIUS_O16+1))
	{
		maxZ = i3End-RADIUS_O16+1 ;
	}

	// Step through the xy-planes
#pragma unroll 17
	for (int iz = gtidz ; iz < maxZ ; iz++)
	{
		// Advance the slice (move the thread-front)
		for (int i = RADIUS_O16 - 1 ; i > 0 ; i--)
			behind[i] = behind[i - 1];

		behind[0] = current;
		current = infront[0];
#pragma unroll 8

		for (int i = 0 ; i < RADIUS_O16 - 1 ; i++)
			infront[i] = infront[i + 1];

		if (validr)
			infront[RADIUS_O16 - 1] = input[inputIndex];

		inputIndex  += stride_z;
		outputIndex += stride_z;
		cg::sync(cta);

		// Note that for the work items on the boundary of the problem, the
		// supplied index when reading the halo (below) may wrap to the
		// previous/next row or even the previous/next xy-plane. This is
		// acceptable since a) we disable the output write for these work
		// items and b) there is at least one xy-plane before/after the
		// current plane, so the access will be within bounds.

		// Update the data slice in the local tile
		// Halo above & below
		if (ltidy < RADIUS_O16)
		{
			tile[ltidy][tx]                  = input[outputIndex - RADIUS_O16 * stride_y];
			tile[ltidy + worky + RADIUS_O16][tx] = input[outputIndex + worky * stride_y];
		}

		// Halo left & right
		if (ltidx < RADIUS_O16)
		{
			tile[ty][ltidx]                  = input[outputIndex - RADIUS_O16];
			tile[ty][ltidx + workx + RADIUS_O16] = input[outputIndex + workx];
		}

		tile[ty][tx] = current;
		cg::sync(cta);

		// Compute the output value
		Myfloat value = fdCoef0 * current ;
#pragma unroll 8

		for (int i = 1 ; i <= RADIUS_O16 ; i++)
		{
			value += stencil[i] * (inv2_d3 * (infront[i-1] + behind[i-1]) +
					inv2_d2 * (tile[ty - i][tx] + tile[ty + i][tx]) +
					inv2_d1 * (tile[ty][tx - i] + tile[ty][tx + i]));
		}

		// Store the output value
		if (validw)
			output[outputIndex] = value;
	}
}


//-------------------------------------------------------------------------------------------------------
// update pressure wavefield (used in propagator) - 3D
// input/output prn
// input prc

__global__ void kernelOpt_computePressureWithFD_3D_O8(const Dim_type dim, const Myint fdOrder, Myfloat *prn, Myfloat *prc, Myfloat *coef,
		const Myfloat inv2_d1, const Myfloat inv2_d2, const Myfloat inv2_d3,
		const Myint n1, const Myint n2, const Myint n3,
		const Myint64 i1Start, const Myint64 i1End, const Myint64 i2Start, const Myint64 i2End, const Myint64 i3Start, const Myint64 i3End,
		const int dimx, const int dimy, const int dimz)
{
#define RADIUS 4
	bool validr = true;
	bool validw = true;

	const int gtidx = blockIdx.x * blockDim.x + threadIdx.x;
	const int gtidy = blockIdx.y * blockDim.y + threadIdx.y;
	const int gtidz = blockIdx.z * dimz ;
	const int ltidx = threadIdx.x;
	const int ltidy = threadIdx.y;
	const int workx = blockDim.x;
	const int worky = blockDim.y;

	// Handle to thread block group
	cg::thread_block cta = cg::this_thread_block();
	__shared__ float tile[MAX_BLOCK_DIM_Y + 2 * RADIUS][MAX_BLOCK_DIM_X + 2 * RADIUS];

	const int stride_y = n1 ;
	const int stride_z = stride_y * n2 ;

	int inputIndex  = gtidz * stride_z ;
	int outputIndex = gtidz * stride_z ;

	// Advance inputIndex to start of inner volume
	inputIndex += i2Start * stride_y + i1Start ;

	// Advance inputIndex to target element
	inputIndex += gtidy * stride_y + gtidx;

	float infront[RADIUS];
	float behind[RADIUS];
	float current;
	const Myfloat fdCoef0 = stencil[0] * (inv2_d1 + inv2_d2 + inv2_d3) ;

	const int tx = ltidx + RADIUS;
	const int ty = ltidy + RADIUS;

	// Check in bounds
	if ((gtidx > i1End) || (gtidy > i2End))
		validr = false;

	if ((gtidx >= dimx) || (gtidy >= dimy))
		validw = false;

	// Preload the "infront" and "behind" data
	for (int i = RADIUS - 2 ; i >= 0 ; i--)
	{
		if (validr)
			behind[i] = prc[inputIndex];

		inputIndex += stride_z;
	}

	if (validr)
		current = prc[inputIndex];

	outputIndex = inputIndex;
	inputIndex += stride_z;

	for (int i = 0 ; i < RADIUS ; i++)
	{
		if (validr)
			infront[i] = prc[inputIndex];

		inputIndex += stride_z;
	}

	// set max element to process along z (n3)
	int maxZ = gtidz + dimz ;
	if (maxZ > (i3End-RADIUS+1)) 
	{
		maxZ = i3End-RADIUS+1 ;	
	}

	// Step through the xy-planes
#pragma unroll 9
	for (int iz = gtidz ; iz < maxZ ; iz++)	
	{
		// Advance the slice (move the thread-front)
		for (int i = RADIUS - 1 ; i > 0 ; i--)
			behind[i] = behind[i - 1];

		behind[0] = current;
		current = infront[0];
#pragma unroll 4

		for (int i = 0 ; i < RADIUS - 1 ; i++)
			infront[i] = infront[i + 1];

		if (validr)
			infront[RADIUS - 1] = prc[inputIndex];

		inputIndex  += stride_z;
		outputIndex += stride_z;
		cg::sync(cta);

		// Note that for the work items on the boundary of the problem, the
		// supplied index when reading the halo (below) may wrap to the
		// previous/next row or even the previous/next xy-plane. This is
		// acceptable since a) we disable the output write for these work
		// items and b) there is at least one xy-plane before/after the
		// current plane, so the access will be within bounds.

		// Update the data slice in the local tile
		// Halo above & below
		if (ltidy < RADIUS)
		{
			tile[ltidy][tx]                  = prc[outputIndex - RADIUS * stride_y];
			tile[ltidy + worky + RADIUS][tx] = prc[outputIndex + worky * stride_y];
		}

		// Halo left & right
		if (ltidx < RADIUS)
		{
			tile[ty][ltidx]                  = prc[outputIndex - RADIUS];
			tile[ty][ltidx + workx + RADIUS] = prc[outputIndex + workx];
		}

		tile[ty][tx] = current;
		cg::sync(cta);

		// Compute the output value
		float value = fdCoef0 * current ;
#pragma unroll 4

		for (int i = 1 ; i <= RADIUS ; i++)
		{
			value += stencil[i] * (inv2_d3 * (infront[i-1] + behind[i-1]) +
					inv2_d2 * (tile[ty - i][tx] + tile[ty + i][tx]) +
					inv2_d1 * (tile[ty][tx - i] + tile[ty][tx + i]));
		}


		// Store the output value
		if (validw)
			prn[outputIndex] = TWO * current - prn[outputIndex] +
			coef[outputIndex] * value ;
	}

}

//-------------------------------------------------------------------------------------------------------

Grid_Cuda_Optim::Grid_Cuda_Optim(Grid_type gridTypeIn) : Grid_Cuda(gridTypeIn)
{
	printDebug(MID_DEBUG, "IN Grid_Cuda_Optim::Grid_Cuda_Optim");

	gridMode = GRID_MODE_CUDA_OPTIM ;

	gpuFDBlkSize1  = UNSPECIFIED ;
	gpuFDBlkSize2  = UNSPECIFIED ;
	gpuFDBlkSize3  = UNSPECIFIED ;
	gpuFDGridSize1 = UNSPECIFIED ;
	gpuFDGridSize2 = UNSPECIFIED ;
	gpuFDGridSize3 = UNSPECIFIED ;

	printDebug(MID_DEBUG, "OUT Grid_Cuda_Optim::Grid_Cuda_Optim");
}

//-------------------------------------------------------------------------------------------------------

Grid_Cuda_Optim::Grid_Cuda_Optim(Grid_type gridTypeIn, Dim_type dimIn,
		Myint64 n1InnerIn, Myint64 n2InnerIn, Myint64 n3InnerIn) : Grid_Cuda(gridTypeIn, dimIn,
				n1InnerIn, n2InnerIn, n3InnerIn)
{
	printDebug(MID_DEBUG, "IN Grid_Cuda_Optim::Grid_Cuda_Optim");

	gridMode = GRID_MODE_CUDA_OPTIM ;

	gpuFDBlkSize1  = UNSPECIFIED ;
	gpuFDBlkSize2  = UNSPECIFIED ;
	gpuFDBlkSize3  = UNSPECIFIED ;
	gpuFDGridSize1 = UNSPECIFIED ;
	gpuFDGridSize2 = UNSPECIFIED ;
	gpuFDGridSize3 = UNSPECIFIED ;

	printDebug(MID_DEBUG, "OUT Grid_Cuda_Optim::Grid_Cuda_Optim");
}

//-------------------------------------------------------------------------------------------------------

Grid_Cuda_Optim::~Grid_Cuda_Optim(void)
{
	printDebug(MID_DEBUG, "IN Grid_Cuda_Optim::~Grid_Cuda_Optim");

	Grid_Cuda::~Grid_Cuda() ;

	printDebug(MID_DEBUG, "OUT Grid_Cuda_Optim::~Grid_Cuda_Optim");
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_Cuda_Optim::initializeGrid(void)
{
	printDebug(FULL_DEBUG, "In Grid_Cuda_Optim::initializeGrid") ;

	Grid_Cuda::initializeGrid() ;

	// Copy the FD coefficients to the device coefficient buffer (constant memory)
	{
		auto FD_coef = getFD_D2coefVector(Config::Instance()->fdOrder) ;
		cudaMemcpyToSymbol(stencil, (void *) FD_coef.data(), FD_coef.size() * sizeof(Myfloat)) ;
		cudaCheckError();
	}

	printDebug(FULL_DEBUG, "Out Grid_Cuda_Optim::initializeGrid") ;
	return(RTN_CODE_OK) ;
}


//-------------------------------------------------------------------------------------------------------

void Grid_Cuda_Optim::info(void)
{
	printDebug(FULL_DEBUG, "IN Grid_Cuda_Optim::info");

	// parent class info
	Grid_Cuda::info() ;

	// determine block size
	gpuFDBlkSize1 = Config::Instance()->cb1 ;
	if (gpuFDBlkSize1 < haloWidth) gpuFDBlkSize1 = haloWidth ;
	if (gpuFDBlkSize1 > MAX_BLOCK_DIM_X) gpuFDBlkSize1 = MAX_BLOCK_DIM_X ;

	gpuFDBlkSize2 = Config::Instance()->cb2 ;
	if (gpuFDBlkSize2 < haloWidth) gpuFDBlkSize2 = haloWidth ;
	if (gpuFDBlkSize2 > MAX_BLOCK_DIM_Y) gpuFDBlkSize2 = MAX_BLOCK_DIM_Y ;

	gpuFDBlkSize3 = 1 ;

	// determine grid size
	gpuFDGridSize1 = n1Inner / gpuFDBlkSize1 ;
	if (n1Inner % gpuFDBlkSize1) gpuFDGridSize1++ ;

	gpuFDGridSize2 = n2Inner / gpuFDBlkSize2 ;
	if (n2Inner % gpuFDBlkSize2) gpuFDGridSize2++ ;

	maxFDBlkSize3 = min(Config::Instance()->cb3, (int) n3Inner) ;
	gpuFDGridSize3 = n3Inner / maxFDBlkSize3 ;
	if (n3Inner % Config::Instance()->cb3) gpuFDGridSize3++ ;

	printInfo(ALL, " - FD computations with 3D blocks") ;

	printInfo(ALL, " FD Blk Size1\t", gpuFDBlkSize1) ;
	printInfo(ALL, " FD Blk Size2\t", gpuFDBlkSize2) ;
	printInfo(ALL, " FD Blk Size3\t", maxFDBlkSize3) ;

	printInfo(ALL, " FD Grid Size1\t", gpuFDGridSize1) ;
	printInfo(ALL, " FD Grid Size2\t", gpuFDGridSize2) ;
	printInfo(ALL, " FD Grid Size3\t", gpuFDGridSize3) ;


	// max number of threads
	struct cudaFuncAttributes funcAttrib ;		
	cudaFuncGetAttributes(&funcAttrib, kernelOpt_FD_LAPLACIAN_O8) ;
	printInfo(MASTER, " Max threads/blk Lapla", funcAttrib.maxThreadsPerBlock) ;
	if ((gpuFDBlkSize1 * gpuFDBlkSize2) > funcAttrib.maxThreadsPerBlock)
	{
		printError("Grid_Cuda_Optim::info, FD block are too large") ;
		return ;
	}

	cudaFuncGetAttributes(&funcAttrib, kernelOpt_computePressureWithFD_3D_O8) ;
	printInfo(MASTER, " Max threads/blk Propa", funcAttrib.maxThreadsPerBlock) ;
	if ((gpuFDBlkSize1 * gpuFDBlkSize2) > funcAttrib.maxThreadsPerBlock)
	{
		printError("Grid_Cuda_Optim::info, FD block are too large") ;
		return ;
	}	

	printDebug(FULL_DEBUG, "OUT Grid_Cuda_Optim::info");
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_Cuda_Optim::FD_LAPLACIAN(Point_type pointType, const Grid& Wgrid, Myint fdOrder)
{
	printDebug(MID_DEBUG, "IN Grid_Cuda_Optim::FD_LAPLACIAN");

	// check grids are same size
	if (this->sameSize(Wgrid) != true)
	{
		printError("Grid_Cuda_Optim::FD_LAPLACIAN, grids have not same size") ;
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

	Myfloat * d_w = ((Grid_Cuda_Optim&) Wgrid).d_grid_3d ;
	Myfloat * d_u = this->d_grid_3d ;		

	int dimx = n1Inner ;
	int dimy = n2Inner ;
	int dimz = maxFDBlkSize3 ;

	dim3 BlkSize(gpuFDBlkSize1, gpuFDBlkSize2, gpuFDBlkSize3) ;
	dim3 GridSize(gpuFDGridSize1, gpuFDGridSize2, gpuFDGridSize3) ;

	if ((dim == DIM3) && (fdOrder == 4))
	{

		kernelOpt_FD_LAPLACIAN_O4<<<GridSize, BlkSize>>>(fdOrder, d_w, d_u, inv2_d1, inv2_d2, inv2_d3,
				n1, n2, n3, i1Start, i1End, i2Start, i2End, i3Start, i3End, dimx, dimy, dimz);
	}
	else if ((dim == DIM3) && (fdOrder == 6))
	{

		kernelOpt_FD_LAPLACIAN_O6<<<GridSize, BlkSize>>>(fdOrder, d_w, d_u,inv2_d1,inv2_d2,inv2_d3,
				n1,n2,n3,i1Start,i1End,i2Start,i2End,i3Start,i3End, dimx, dimy, dimz);
	}
	else if ((dim == DIM3) && (fdOrder == 8))
	{

		kernelOpt_FD_LAPLACIAN_O8<<<GridSize, BlkSize>>>(fdOrder, d_w, d_u,inv2_d1,inv2_d2,inv2_d3,
				n1,n2,n3,i1Start,i1End,i2Start,i2End,i3Start,i3End, dimx, dimy, dimz);
	}
	else if ((dim == DIM3) && (fdOrder == 10))
	{

		kernelOpt_FD_LAPLACIAN_O10<<<GridSize, BlkSize>>>(fdOrder, d_w, d_u,inv2_d1,inv2_d2,inv2_d3,
				n1,n2,n3,i1Start,i1End,i2Start,i2End,i3Start,i3End, dimx, dimy, dimz);
	}
	else if ((dim == DIM3) && (fdOrder == 12))
	{

		kernelOpt_FD_LAPLACIAN_O12<<<GridSize, BlkSize>>>(fdOrder, d_w, d_u,inv2_d1,inv2_d2,inv2_d3,
				n1,n2,n3,i1Start,i1End,i2Start,i2End,i3Start,i3End, dimx, dimy, dimz);
	}
	else if ((dim == DIM3) && (fdOrder == 14))
	{

		kernelOpt_FD_LAPLACIAN_O14<<<GridSize, BlkSize>>>(fdOrder, d_w, d_u,inv2_d1,inv2_d2,inv2_d3,
				n1,n2,n3,i1Start,i1End,i2Start,i2End,i3Start,i3End, dimx, dimy, dimz);
	}
	else if ((dim == DIM3) && (fdOrder == 16))
	{

		kernelOpt_FD_LAPLACIAN_O16<<<GridSize, BlkSize>>>(fdOrder, d_w, d_u,inv2_d1,inv2_d2,inv2_d3,
				n1,n2,n3,i1Start,i1End,i2Start,i2End,i3Start,i3End, dimx, dimy, dimz);
	}
	else
	{
		Grid_Cuda::FD_LAPLACIAN(pointType, Wgrid, fdOrder) ;
	}

	cudaCheckError();
	cudaDeviceSynchronize();

	printDebug(MID_DEBUG, "OUT Grid_Cuda_Optim::FD_LAPLACIAN");
	return(RTN_CODE_OK) ;
}

Rtn_code Grid_Cuda_Optim::computePressureWithFD(Grid& prcGridIn, Grid& coefGridIn, Myint fdOrder)
{

	printDebug(FULL_DEBUG, "In Grid_Cuda_Optim::computePressureWithFD") ;

	// check grids are same size
	if (this->sameSize(prcGridIn) != true)
	{
		printError("In Grid_Cuda_Optim::computePressureWithFD, grids have not same size") ;
		return(RTN_CODE_KO) ;
	}

	//pointType
	Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
	getGridIndex(INNER_POINTS, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End);

	const Myfloat inv_d1  = Myfloat(1.0) / d1 ;
	const Myfloat inv_d2  = Myfloat(1.0) / d2 ;
	const Myfloat inv_d3  = Myfloat(1.0) / d3 ;

	const Myfloat inv2_d1 = inv_d1 * inv_d1 ;
	const Myfloat inv2_d2 = inv_d2 * inv_d2 ;
	const Myfloat inv2_d3 = inv_d3 * inv_d3 ;

	Myfloat *prc_d_grid_3d = ((Grid_Cuda_Optim&) prcGridIn).d_grid_3d ;
	Myfloat *coef_d_grid_3d = ((Grid_Cuda_Optim&) coefGridIn).d_grid_3d ;

	int dimx = n1Inner ;
	int dimy = n2Inner ;
	int dimz = maxFDBlkSize3 ;

	dim3 BlkSize(gpuFDBlkSize1, gpuFDBlkSize2, gpuFDBlkSize3) ;
	dim3 GridSize(gpuFDGridSize1, gpuFDGridSize2, gpuFDGridSize3) ;

	if ((dim == DIM3) && (fdOrder == 8))
	{			   
		kernelOpt_computePressureWithFD_3D_O8<<<GridSize, BlkSize>>>(dim, fdOrder, d_grid_3d, prc_d_grid_3d, coef_d_grid_3d,inv2_d1,inv2_d2,inv2_d3,
				n1,n2,n3,i1Start,i1End,i2Start,i2End,i3Start,i3End, dimx, dimy, dimz);
	}
	else
	{
		Grid_Cuda::computePressureWithFD(prcGridIn, coefGridIn, fdOrder) ;
	}

	cudaCheckError();
	cudaDeviceSynchronize();

	printDebug(FULL_DEBUG, "Out Grid_Cuda_Optim::computePressureWithFD") ;
	return(RTN_CODE_OK) ;
}


} // namespace hpcscan
