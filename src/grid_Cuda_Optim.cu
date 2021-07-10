
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

#include "mpi.h"

#include "config.h"
#include "constant.h"
#include "fdm.h"
#include "global.h"
#include "output_report.h"

using namespace std;

// start code from NVIDIA
#include <cooperative_groups.h>

namespace cg = cooperative_groups;

// Note: If you change the RADIUS, you should also change the unrolling below
#define RADIUS 4

#define k_blockDimX    32
#define k_blockDimMaxY 16
#define k_blockSizeMin 128
#define k_blockSizeMax (k_blockDimX * k_blockDimMaxY)
// end code from NVIDIA

namespace hpcscan {

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

__constant__ float stencil[RADIUS + 1];

//-------------------------------------------------------------------------------------------------------
// compute derivative along axis 1
// input u
// output w

__global__ void kernelOpt_FD_N1_O8(const Myint fdOrder, Myfloat *output, Myfloat *input,
		const Myfloat inv2_d1, const Myfloat inv2_d2, const Myfloat inv2_d3,
		const Myint n1, const Myint n2, const Myint n3,
		const Myint64 i1Start, const Myint64 i1End, const Myint64 i2Start, const Myint64 i2End, const Myint64 i3Start, const Myint64 i3End,
		const int dimx, const int dimy, const int dimz)
{
	bool validr = true;
	bool validw = true;

	const int gtidx = blockIdx.x * blockDim.x + threadIdx.x;
	const int gtidy = blockIdx.y * blockDim.y + threadIdx.y;
	const int ltidx = threadIdx.x;
	const int ltidy = threadIdx.y;
	const int workx = blockDim.x;
	const int worky = blockDim.y;

	// Handle to thread block group
	cg::thread_block cta = cg::this_thread_block();
	__shared__ float tile[k_blockDimMaxY + 2 * RADIUS][k_blockDimX + 2 * RADIUS];

	const int stride_y = dimx + 2 * RADIUS;
	const int stride_z = stride_y * (dimy + 2 * RADIUS);

	int inputIndex  = 0;
	int outputIndex = 0;

	// Advance inputIndex to start of inner volume
	inputIndex += RADIUS * stride_y + RADIUS;

	// Advance inputIndex to target element
	inputIndex += gtidy * stride_y + gtidx;

	float infront[RADIUS];
	float behind[RADIUS];
	float current;

	const int tx = ltidx + RADIUS;
	const int ty = ltidy + RADIUS;

	// Check in bounds
	if ((gtidx >= dimx + RADIUS) || (gtidy >= dimy + RADIUS))
		validr = false;

	if ((gtidx >= dimx) || (gtidy >= dimy))
		validw = false;

	// Preload the "infront" and "behind" data
	for (int i = RADIUS - 2 ; i >= 0 ; i--)
	{
		if (validr)
			behind[i] = input[inputIndex];

		inputIndex += stride_z;
	}

	if (validr)
		current = input[inputIndex];

	outputIndex = inputIndex;
	inputIndex += stride_z;

	for (int i = 0 ; i < RADIUS ; i++)
	{
		if (validr)
			infront[i] = input[inputIndex];

		inputIndex += stride_z;
	}

	// Step through the xy-planes
#pragma unroll 9

	for (int iz = 0 ; iz < dimz ; iz++)
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
			infront[RADIUS - 1] = input[inputIndex];

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
			tile[ltidy][tx]                  = input[outputIndex - RADIUS * stride_y];
			tile[ltidy + worky + RADIUS][tx] = input[outputIndex + worky * stride_y];
		}

		// Halo left & right
		if (ltidx < RADIUS)
		{
			tile[ty][ltidx]                  = input[outputIndex - RADIUS];
			tile[ty][ltidx + workx + RADIUS] = input[outputIndex + workx];
		}

		tile[ty][tx] = current;
		cg::sync(cta);

		// Compute the output value
		float value = stencil[0] * current ;
#pragma unroll 4

		for (int i = 1 ; i <= RADIUS ; i++)
		{
			//value += stencil[i] * (infront[i-1] + behind[i-1] + tile[ty - i][tx] + tile[ty + i][tx] + tile[ty][tx - i] + tile[ty][tx + i]);
			value += stencil[i] * (tile[ty][tx - i] + tile[ty][tx + i]) ; // d1
		}

		// Store the output value
		if (validw)
			output[outputIndex] = value;
	}
}

//-------------------------------------------------------------------------------------------------------
// compute derivative along axis 2
// input u
// output w

__global__ void kernelOpt_FD_N2_O8(const Myint fdOrder, Myfloat *output, Myfloat *input,
		const Myfloat inv2_d1, const Myfloat inv2_d2, const Myfloat inv2_d3,
		const Myint n1, const Myint n2, const Myint n3,
		const Myint64 i1Start, const Myint64 i1End, const Myint64 i2Start, const Myint64 i2End, const Myint64 i3Start, const Myint64 i3End,
		const int dimx, const int dimy, const int dimz)
{
	bool validr = true;
	bool validw = true;

	const int gtidx = blockIdx.x * blockDim.x + threadIdx.x;
	const int gtidy = blockIdx.y * blockDim.y + threadIdx.y;
	const int ltidx = threadIdx.x;
	const int ltidy = threadIdx.y;
	const int workx = blockDim.x;
	const int worky = blockDim.y;

	// Handle to thread block group
	cg::thread_block cta = cg::this_thread_block();
	__shared__ float tile[k_blockDimMaxY + 2 * RADIUS][k_blockDimX + 2 * RADIUS];

	const int stride_y = dimx + 2 * RADIUS;
	const int stride_z = stride_y * (dimy + 2 * RADIUS);

	int inputIndex  = 0;
	int outputIndex = 0;

	// Advance inputIndex to start of inner volume
	inputIndex += RADIUS * stride_y + RADIUS;

	// Advance inputIndex to target element
	inputIndex += gtidy * stride_y + gtidx;

	float infront[RADIUS];
	float behind[RADIUS];
	float current;

	const int tx = ltidx + RADIUS;
	const int ty = ltidy + RADIUS;

	// Check in bounds
	if ((gtidx >= dimx + RADIUS) || (gtidy >= dimy + RADIUS))
		validr = false;

	if ((gtidx >= dimx) || (gtidy >= dimy))
		validw = false;

	// Preload the "infront" and "behind" data
	for (int i = RADIUS - 2 ; i >= 0 ; i--)
	{
		if (validr)
			behind[i] = input[inputIndex];

		inputIndex += stride_z;
	}

	if (validr)
		current = input[inputIndex];

	outputIndex = inputIndex;
	inputIndex += stride_z;

	for (int i = 0 ; i < RADIUS ; i++)
	{
		if (validr)
			infront[i] = input[inputIndex];

		inputIndex += stride_z;
	}

	// Step through the xy-planes
#pragma unroll 9

	for (int iz = 0 ; iz < dimz ; iz++)
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
			infront[RADIUS - 1] = input[inputIndex];

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
			tile[ltidy][tx]                  = input[outputIndex - RADIUS * stride_y];
			tile[ltidy + worky + RADIUS][tx] = input[outputIndex + worky * stride_y];
		}

		// Halo left & right
		if (ltidx < RADIUS)
		{
			tile[ty][ltidx]                  = input[outputIndex - RADIUS];
			tile[ty][ltidx + workx + RADIUS] = input[outputIndex + workx];
		}

		tile[ty][tx] = current;
		cg::sync(cta);

		// Compute the output value
		float value = stencil[0] * current ;
#pragma unroll 4

		for (int i = 1 ; i <= RADIUS ; i++)
		{
			//value += stencil[i] * (infront[i-1] + behind[i-1] + tile[ty - i][tx] + tile[ty + i][tx] + tile[ty][tx - i] + tile[ty][tx + i]);
			value += stencil[i] * (tile[ty - i][tx] + tile[ty + i][tx]); // d2
		}

		// Store the output value
		if (validw)
			output[outputIndex] = value;
	}
}

//-------------------------------------------------------------------------------------------------------
// compute derivative along axis 3
// input u
// output w

__global__ void kernelOpt_FD_N3_O8(const Myint fdOrder, Myfloat *output, Myfloat *input,
		const Myfloat inv2_d1, const Myfloat inv2_d2, const Myfloat inv2_d3,
		const Myint n1, const Myint n2, const Myint n3,
		const Myint64 i1Start, const Myint64 i1End, const Myint64 i2Start, const Myint64 i2End, const Myint64 i3Start, const Myint64 i3End,
		const int dimx, const int dimy, const int dimz)
{
	bool validr = true;
	bool validw = true;

	const int gtidx = blockIdx.x * blockDim.x + threadIdx.x;
	const int gtidy = blockIdx.y * blockDim.y + threadIdx.y;
	const int ltidx = threadIdx.x;
	const int ltidy = threadIdx.y;
	const int workx = blockDim.x;
	const int worky = blockDim.y;

	// Handle to thread block group
	cg::thread_block cta = cg::this_thread_block();
	__shared__ float tile[k_blockDimMaxY + 2 * RADIUS][k_blockDimX + 2 * RADIUS];

	const int stride_y = dimx + 2 * RADIUS;
	const int stride_z = stride_y * (dimy + 2 * RADIUS);

	int inputIndex  = 0;
	int outputIndex = 0;

	// Advance inputIndex to start of inner volume
	inputIndex += RADIUS * stride_y + RADIUS;

	// Advance inputIndex to target element
	inputIndex += gtidy * stride_y + gtidx;

	float infront[RADIUS];
	float behind[RADIUS];
	float current;

	const int tx = ltidx + RADIUS;
	const int ty = ltidy + RADIUS;

	// Check in bounds
	if ((gtidx >= dimx + RADIUS) || (gtidy >= dimy + RADIUS))
		validr = false;

	if ((gtidx >= dimx) || (gtidy >= dimy))
		validw = false;

	// Preload the "infront" and "behind" data
	for (int i = RADIUS - 2 ; i >= 0 ; i--)
	{
		if (validr)
			behind[i] = input[inputIndex];

		inputIndex += stride_z;
	}

	if (validr)
		current = input[inputIndex];

	outputIndex = inputIndex;
	inputIndex += stride_z;

	for (int i = 0 ; i < RADIUS ; i++)
	{
		if (validr)
			infront[i] = input[inputIndex];

		inputIndex += stride_z;
	}

	// Step through the xy-planes
#pragma unroll 9

	for (int iz = 0 ; iz < dimz ; iz++)
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
			infront[RADIUS - 1] = input[inputIndex];

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
		//if (ltidy < RADIUS)
		//{
		//	tile[ltidy][tx]                  = input[outputIndex - RADIUS * stride_y];
		//	tile[ltidy + worky + RADIUS][tx] = input[outputIndex + worky * stride_y];
		//}

		// Halo left & right
		//if (ltidx < RADIUS)
		//{
		//	tile[ty][ltidx]                  = input[outputIndex - RADIUS];
		//	tile[ty][ltidx + workx + RADIUS] = input[outputIndex + workx];
		//}

		tile[ty][tx] = current;
		cg::sync(cta);

		// Compute the output value
		float value = stencil[0] * current ;
#pragma unroll 4

		for (int i = 1 ; i <= RADIUS ; i++)
		{
			//value += stencil[i] * (infront[i-1] + behind[i-1] + tile[ty - i][tx] + tile[ty + i][tx] + tile[ty][tx - i] + tile[ty][tx + i]);
			value += stencil[i] * (infront[i-1] + behind[i-1]); // d3
		}

		// Store the output value
		if (validw)
			output[outputIndex] = value;
	}
}


//-------------------------------------------------------------------------------------------------------
// compute Laplacian
// input u
// output w

__global__ void kernelOpt_FD_LAPLACIAN_O8(const Myint fdOrder, Myfloat *output, Myfloat *input,
		const Myfloat inv2_d1, const Myfloat inv2_d2, const Myfloat inv2_d3,
		const Myint n1, const Myint n2, const Myint n3,
		const Myint64 i1Start, const Myint64 i1End, const Myint64 i2Start, const Myint64 i2End, const Myint64 i3Start, const Myint64 i3End,
		const int dimx, const int dimy, const int dimz)
{
	bool validr = true;
	bool validw = true;

	const int gtidx = blockIdx.x * blockDim.x + threadIdx.x;
	const int gtidy = blockIdx.y * blockDim.y + threadIdx.y;
	const int ltidx = threadIdx.x;
	const int ltidy = threadIdx.y;
	const int workx = blockDim.x;
	const int worky = blockDim.y;

	// Handle to thread block group
	cg::thread_block cta = cg::this_thread_block();
	__shared__ float tile[k_blockDimMaxY + 2 * RADIUS][k_blockDimX + 2 * RADIUS];

	const int stride_y = dimx + 2 * RADIUS;
	const int stride_z = stride_y * (dimy + 2 * RADIUS);

	int inputIndex  = 0;
	int outputIndex = 0;

	// Advance inputIndex to start of inner volume
	inputIndex += RADIUS * stride_y + RADIUS;

	// Advance inputIndex to target element
	inputIndex += gtidy * stride_y + gtidx;

	float infront[RADIUS];
	float behind[RADIUS];
	float current;

	const int tx = ltidx + RADIUS;
	const int ty = ltidy + RADIUS;

	// Check in bounds
	if ((gtidx >= dimx + RADIUS) || (gtidy >= dimy + RADIUS))
		validr = false;

	if ((gtidx >= dimx) || (gtidy >= dimy))
		validw = false;

	// Preload the "infront" and "behind" data
	for (int i = RADIUS - 2 ; i >= 0 ; i--)
	{
		if (validr)
			behind[i] = input[inputIndex];

		inputIndex += stride_z;
	}

	if (validr)
		current = input[inputIndex];

	outputIndex = inputIndex;
	inputIndex += stride_z;

	for (int i = 0 ; i < RADIUS ; i++)
	{
		if (validr)
			infront[i] = input[inputIndex];

		inputIndex += stride_z;
	}

	// Step through the xy-planes
#pragma unroll 9

	for (int iz = 0 ; iz < dimz ; iz++)
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
			infront[RADIUS - 1] = input[inputIndex];

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
			tile[ltidy][tx]                  = input[outputIndex - RADIUS * stride_y];
			tile[ltidy + worky + RADIUS][tx] = input[outputIndex + worky * stride_y];
		}

		// Halo left & right
		if (ltidx < RADIUS)
		{
			tile[ty][ltidx]                  = input[outputIndex - RADIUS];
			tile[ty][ltidx + workx + RADIUS] = input[outputIndex + workx];
		}

		tile[ty][tx] = current;
		cg::sync(cta);

		// Compute the output value
		float value = stencil[0] * current * 3 ;
#pragma unroll 4

		for (int i = 1 ; i <= RADIUS ; i++)
		{
			value += stencil[i] * (infront[i-1] + behind[i-1] + tile[ty - i][tx] + tile[ty + i][tx] + tile[ty][tx - i] + tile[ty][tx + i]);
		}

		// Store the output value
		if (validw)
			output[outputIndex] = value;
	}
}

//-------------------------------------------------------------------------------------------------------

Grid_Cuda_Optim::Grid_Cuda_Optim(Grid_type gridTypeIn) : Grid_Cuda(gridTypeIn)
{
	printDebug(MID_DEBUG, "IN Grid_Cuda_Optim::Grid_Cuda_Optim");

	gridMode = GRID_MODE_CUDA_OPTIM ;

	printDebug(MID_DEBUG, "OUT Grid_Cuda_Optim::Grid_Cuda_Optim");
}

//-------------------------------------------------------------------------------------------------------

Grid_Cuda_Optim::Grid_Cuda_Optim(Grid_type gridTypeIn, Dim_type dimIn,
		Myint64 n1InnerIn, Myint64 n2InnerIn, Myint64 n3InnerIn) : Grid_Cuda(gridTypeIn, dimIn,
				n1InnerIn, n2InnerIn, n3InnerIn)
{
	printDebug(MID_DEBUG, "IN Grid_Cuda_Optim::Grid_Cuda_Optim");

	gridMode = GRID_MODE_CUDA_OPTIM ;

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

	printDebug(FULL_DEBUG, "Out Grid_Cuda_Optim::initializeGrid") ;
	return(RTN_CODE_OK) ;
}


//-------------------------------------------------------------------------------------------------------

void Grid_Cuda_Optim::info(void)
{
	printDebug(FULL_DEBUG, "IN Grid_Cuda_Optim::info");

	// parent class info
	Grid_Cuda::info() ;

	// additional info
	printInfo(MASTER, " Block size 1", gpuBlkSize1) ;
	printInfo(MASTER, " Block size 2", gpuBlkSize2) ;
	printInfo(MASTER, " Block size 3", gpuBlkSize3) ;

	printDebug(FULL_DEBUG, "OUT Grid_Cuda_Optim::info");
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_Cuda_Optim::FD_D2_N1(Point_type pointType, const Grid& Wgrid, Myint fdOrder)
{
	printDebug(MID_DEBUG, "IN Grid_Cuda_Optim::FD_D2_N1");

	// TO DO
	Grid_Cuda::FD_D2_N1(pointType, Wgrid, fdOrder) ;
	// check grids are same size
	if (this->sameSize(Wgrid) != true)
	{
		printError("Grid_Cuda_Optim::FD_D2_N1, grids have not same size") ;
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

	if (fdOrder == 8)
	{
		gpuBlkSize3 = 1 ;
		int GridSize1 = n1 / gpuBlkSize1 + 1 ;
		int GridSize2 = n2 / gpuBlkSize2 + 1 ;
		//int GridSize3 = n3 / gpuBlkSize3 + 1 ;
		int GridSize3 = 1 ;

		dim3 BlkSize(gpuBlkSize1, gpuBlkSize2, gpuBlkSize3) ;
		dim3 GridSize(GridSize1, GridSize2, GridSize3) ;

		int dimx = n1Inner ;
		int dimy = n2Inner ;
		int dimz = n3Inner ;

		// Copy the coefficients to the device coefficient buffer
		float* coeff = (float *)malloc((RADIUS + 1) * sizeof(float));

		// Create coefficients
		coeff[0] = FD_D2_O8_A0 * inv2_d1 ;
		coeff[1] = FD_D2_O8_A1 * inv2_d1 ;
		coeff[2] = FD_D2_O8_A2 * inv2_d1 ;
		coeff[3] = FD_D2_O8_A3 * inv2_d1 ;
		coeff[4] = FD_D2_O8_A4 * inv2_d1 ;

		//checkCudaErrors(cudaMemcpyToSymbol(stencil, (void *)coeff, (RADIUS + 1) * sizeof(float)));
		cudaMemcpyToSymbol(stencil, (void *)coeff, (RADIUS + 1) * sizeof(float)) ;

		kernelOpt_FD_N1_O8<<<GridSize, BlkSize>>>(fdOrder, d_w, d_u,inv2_d1,inv2_d2,inv2_d3,
				n1,n2,n3,i1Start,i1End,i2Start,i2End,i3Start,i3End, dimx, dimy, dimz);

		cudaCheckError();
		cudaDeviceSynchronize();
	}
	else
	{
		Grid_Cuda::FD_D2_N1(pointType, Wgrid, fdOrder) ;
	}

	printDebug(MID_DEBUG, "OUT Grid_Cuda_Optim::FD_D2_N1");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_Cuda_Optim::FD_D2_N2(Point_type pointType, const Grid& Wgrid, Myint fdOrder)
{
	printDebug(MID_DEBUG, "IN Grid_Cuda_Optim::FD_D2_N2");

	// check grids are same size
	if (this->sameSize(Wgrid) != true)
	{
		printError("Grid_Cuda_Optim::FD_D2_N2, grids have not same size") ;
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

	if (fdOrder == 8)
	{
		gpuBlkSize3 = 1 ;
		int GridSize1 = n1 / gpuBlkSize1 + 1 ;
		int GridSize2 = n2 / gpuBlkSize2 + 1 ;
		//int GridSize3 = n3 / gpuBlkSize3 + 1 ;
		int GridSize3 = 1 ;

		dim3 BlkSize(gpuBlkSize1, gpuBlkSize2, gpuBlkSize3) ;
		dim3 GridSize(GridSize1, GridSize2, GridSize3) ;

		int dimx = n1Inner ;
		int dimy = n2Inner ;
		int dimz = n3Inner ;

		// Copy the coefficients to the device coefficient buffer
		float* coeff = (float *)malloc((RADIUS + 1) * sizeof(float));

		// Create coefficients
		coeff[0] = FD_D2_O8_A0 * inv2_d2 ;
		coeff[1] = FD_D2_O8_A1 * inv2_d2 ;
		coeff[2] = FD_D2_O8_A2 * inv2_d2 ;
		coeff[3] = FD_D2_O8_A3 * inv2_d2 ;
		coeff[4] = FD_D2_O8_A4 * inv2_d2 ;

		//checkCudaErrors(cudaMemcpyToSymbol(stencil, (void *)coeff, (RADIUS + 1) * sizeof(float)));
		cudaMemcpyToSymbol(stencil, (void *)coeff, (RADIUS + 1) * sizeof(float)) ;

		kernelOpt_FD_N2_O8<<<GridSize, BlkSize>>>(fdOrder, d_w, d_u,inv2_d1,inv2_d2,inv2_d3,
				n1,n2,n3,i1Start,i1End,i2Start,i2End,i3Start,i3End, dimx, dimy, dimz);

		cudaCheckError();
		cudaDeviceSynchronize();
	}
	else
	{
		Grid_Cuda::FD_D2_N2(pointType, Wgrid, fdOrder) ;
	}

	printDebug(MID_DEBUG, "OUT Grid_Cuda_Optim::FD_D2_N2");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_Cuda_Optim::FD_D2_N3(Point_type pointType, const Grid& Wgrid, Myint fdOrder)
{
	printDebug(MID_DEBUG, "IN Grid_Cuda_Optim::FD_D2_N3");

	// check grids are same size
	if (this->sameSize(Wgrid) != true)
	{
		printError("Grid_Cuda_Optim::FD_D2_N3, grids have not same size") ;
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

	if (fdOrder == 8)
	{
		gpuBlkSize3 = 1 ;
		int GridSize1 = n1 / gpuBlkSize1 + 1 ;
		int GridSize2 = n2 / gpuBlkSize2 + 1 ;
		//int GridSize3 = n3 / gpuBlkSize3 + 1 ;
		int GridSize3 = 1 ;

		dim3 BlkSize(gpuBlkSize1, gpuBlkSize2, gpuBlkSize3) ;
		dim3 GridSize(GridSize1, GridSize2, GridSize3) ;

		int dimx = n1Inner ;
		int dimy = n2Inner ;
		int dimz = n3Inner ;

		// Copy the coefficients to the device coefficient buffer
		float* coeff = (float *)malloc((RADIUS + 1) * sizeof(float));

		// Create coefficients
		coeff[0] = FD_D2_O8_A0 * inv2_d3 ;
		coeff[1] = FD_D2_O8_A1 * inv2_d3 ;
		coeff[2] = FD_D2_O8_A2 * inv2_d3 ;
		coeff[3] = FD_D2_O8_A3 * inv2_d3 ;
		coeff[4] = FD_D2_O8_A4 * inv2_d3 ;

		//checkCudaErrors(cudaMemcpyToSymbol(stencil, (void *)coeff, (RADIUS + 1) * sizeof(float)));
		cudaMemcpyToSymbol(stencil, (void *)coeff, (RADIUS + 1) * sizeof(float)) ;

		kernelOpt_FD_N3_O8<<<GridSize, BlkSize>>>(fdOrder, d_w, d_u,inv2_d1,inv2_d2,inv2_d3,
				n1,n2,n3,i1Start,i1End,i2Start,i2End,i3Start,i3End, dimx, dimy, dimz);

		cudaCheckError();
		cudaDeviceSynchronize();
	}
	else
	{
		Grid_Cuda::FD_D2_N3(pointType, Wgrid, fdOrder) ;
	}

	printDebug(MID_DEBUG, "OUT Grid_Cuda_Optim::FD_D2_N3");
	return(RTN_CODE_OK) ;
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

	if ((dim == DIM3) && (fdOrder == 8))
	{
		gpuBlkSize3 = 1 ;
		int GridSize1 = n1 / gpuBlkSize1 + 1 ;
		int GridSize2 = n2 / gpuBlkSize2 + 1 ;
		//int GridSize3 = n3 / gpuBlkSize3 + 1 ;
		int GridSize3 = 1 ;

		dim3 BlkSize(gpuBlkSize1, gpuBlkSize2, gpuBlkSize3) ;
		dim3 GridSize(GridSize1, GridSize2, GridSize3) ;

		int dimx = n1Inner ;
		int dimy = n2Inner ;
		int dimz = n3Inner ;

		// Copy the coefficients to the device coefficient buffer
		float* coeff = (float *)malloc((RADIUS + 1) * sizeof(float));

		// Create coefficients
		coeff[0] = FD_D2_O8_A0 * inv2_d1 ;
		coeff[1] = FD_D2_O8_A1 * inv2_d1 ;
		coeff[2] = FD_D2_O8_A2 * inv2_d1 ;
		coeff[3] = FD_D2_O8_A3 * inv2_d1 ;
		coeff[4] = FD_D2_O8_A4 * inv2_d1 ;

		//checkCudaErrors(cudaMemcpyToSymbol(stencil, (void *)coeff, (RADIUS + 1) * sizeof(float)));
		cudaMemcpyToSymbol(stencil, (void *)coeff, (RADIUS + 1) * sizeof(float)) ;

		kernelOpt_FD_LAPLACIAN_O8<<<GridSize, BlkSize>>>(fdOrder, d_w, d_u,inv2_d1,inv2_d2,inv2_d3,
				n1,n2,n3,i1Start,i1End,i2Start,i2End,i3Start,i3End, dimx, dimy, dimz);

		cudaCheckError();
		cudaDeviceSynchronize();
	}
	else
	{
		Grid_Cuda::FD_LAPLACIAN(pointType, Wgrid, fdOrder) ;
	}

	printDebug(MID_DEBUG, "OUT Grid_Cuda_Optim::FD_LAPLACIAN");
	return(RTN_CODE_OK) ;
}

} // namespace hpcscan
