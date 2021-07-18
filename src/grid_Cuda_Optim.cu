
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

// start code from code samples
#include <cooperative_groups.h>

namespace cg = cooperative_groups;

// TODO get rid of these hardcoded values
#define MaxBlockDimX 64
#define MaxBlockDimY 16
// end code from code samples

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

// With max fd order = 16, there are at most 9 fd coef in the stencil
__constant__ float stencil[9];

//-------------------------------------------------------------------------------------------------------
// compute Laplacian
//   u
// output w

__launch_bounds__(MaxBlockDimY*MaxBlockDimX)

__global__ void kernelOpt_FD_LAPLACIAN_O8(const Myint fdOrder, Myfloat *output, Myfloat *input,
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
	__shared__ float tile[MaxBlockDimY + 2 * RADIUS][MaxBlockDimX + 2 * RADIUS];
	
	const int stride_y = n1 ;
	const int stride_z = stride_y * n2 ;
	
	//const int offsetZ = (int) i3Start - RADIUS ;
	//const int offsetZ = 0 ;
	
	//printf("offsetZ=%d\n", offsetZ) ;

	//int inputIndex  = 0;
	//int outputIndex = 0;
	int inputIndex  = (gtidz) * stride_z ;
    int outputIndex = (gtidz) * stride_z ;

	// Advance inputIndex to start of inner volume
	inputIndex += i2Start * stride_y + i1Start ;

	// Advance inputIndex to target element
	inputIndex += gtidy * stride_y + gtidx;

	Myfloat infront[RADIUS];
	Myfloat behind[RADIUS];
	Myfloat current;
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

	// set max element to process along z (n3)
	//const int maxZ = min(gtidz + dimz, (int) i3End-RADIUS+1) ;
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
		Myfloat value = fdCoef0 * current ;
#pragma unroll 4

		for (int i = 1 ; i <= RADIUS ; i++)
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
	__shared__ float tile[MaxBlockDimY + 2 * RADIUS][MaxBlockDimX + 2 * RADIUS];

	const int stride_y = dimx + 2 * RADIUS;
	const int stride_z = stride_y * (dimy + 2 * RADIUS);
	
	const int offsetZ = (int) i3Start - RADIUS ;

	//int inputIndex  = 0;
	//int outputIndex = 0;
	int inputIndex  = gtidz * stride_z ;
    int outputIndex = gtidz * stride_z ;

	// Advance inputIndex to start of inner volume
	inputIndex += RADIUS * stride_y + RADIUS;

	// Advance inputIndex to target element
	inputIndex += gtidy * stride_y + gtidx;

	float infront[RADIUS];
	float behind[RADIUS];
	float current;
	const Myfloat fdCoef0 = stencil[0] * (inv2_d1 + inv2_d2 + inv2_d3) ;

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
	
	// Copy the coefficients to the device coefficient buffer
		float* coeff = (float *)malloc((RADIUS + 1) * sizeof(float));

		// Create coefficients
		coeff[0] = FD_D2_O8_A0 ;
		coeff[1] = FD_D2_O8_A1 ;
		coeff[2] = FD_D2_O8_A2 ;
		coeff[3] = FD_D2_O8_A3 ;
		coeff[4] = FD_D2_O8_A4 ;

		//checkCudaErrors(cudaMemcpyToSymbol(stencil, (void *)coeff, (RADIUS + 1) * sizeof(float)));
		cudaMemcpyToSymbol(stencil, (void *)coeff, (RADIUS + 1) * sizeof(float)) ;

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
	printInfo(MASTER, " FD block size 1", Config::Instance()->cb1) ;
	printInfo(MASTER, " FD block size 2", Config::Instance()->cb2) ;
	printInfo(MASTER, " FD block size 3", Config::Instance()->cb3) ;	
	
			// TODO overrides gpuBlkSize
		int gpuBlkSize1 = Config::Instance()->cb1 ;
		int gpuBlkSize2 = Config::Instance()->cb2 ;				
		int gpuBlkSize3 = 1 ;
		
		int dimx = n1Inner ;
		int dimy = n2Inner ;
		int dimz = -1 ;	
		
		int GridSize1 = n1Inner / gpuBlkSize1 ;
		if (n1Inner % gpuBlkSize1) GridSize1++ ;
		int GridSize2 = n2Inner / gpuBlkSize2 ;
		if (n2Inner % gpuBlkSize2) GridSize2++ ;
		int GridSize3 = -1 ;
		if (Config::Instance()->cb3 >= n3Inner)
		{
		GridSize3 = 1 ;
		dimz      = n3Inner ;
		}
		else
		{		
		GridSize3 = n3Inner / Config::Instance()->cb3 ;
		if (n3Inner % Config::Instance()->cb3) GridSize3++ ;
		dimz      = Config::Instance()->cb3 ;
		}

		printInfo(ALL, "GridSize1", GridSize1) ;
		printInfo(ALL, "GridSize2", GridSize2) ;
		printInfo(ALL, "GridSize3", GridSize3) ;
		
		printInfo(ALL, "gpuBlkSize1", gpuBlkSize1) ;
		printInfo(ALL, "gpuBlkSize2", gpuBlkSize2) ;
		printInfo(ALL, "gpuBlkSize3", gpuBlkSize3) ;
	
	
	// max number of threads
	struct cudaFuncAttributes funcAttrib ;		
	cudaFuncGetAttributes(&funcAttrib, kernelOpt_FD_LAPLACIAN_O8) ;
	printInfo(MASTER, " Max threads/blk Lapla", funcAttrib.maxThreadsPerBlock) ;
	if ((Config::Instance()->cb1 * Config::Instance()->cb2) > funcAttrib.maxThreadsPerBlock)
	{
	   printError("Grid_Cuda_Optim::info, FD block are too large") ;
	   return ;
	}
	
	cudaFuncGetAttributes(&funcAttrib, kernelOpt_computePressureWithFD_3D_O8) ;
	printInfo(MASTER, " Max threads/blk Propa", funcAttrib.maxThreadsPerBlock) ;
	if ((Config::Instance()->cb1 * Config::Instance()->cb2) > funcAttrib.maxThreadsPerBlock)
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

	if ((dim == DIM3) && (fdOrder == 8))
	{		  	    
		// TODO overrides gpuBlkSize
		int gpuBlkSize1 = Config::Instance()->cb1 ;
		int gpuBlkSize2 = Config::Instance()->cb2 ;				
		int gpuBlkSize3 = 1 ;
		
		int dimx = n1Inner ;
		int dimy = n2Inner ;
		int dimz = -1 ;	
		
		int GridSize1 = n1Inner / gpuBlkSize1 ;
		if (n1Inner % gpuBlkSize1) GridSize1++ ;
		int GridSize2 = n2Inner / gpuBlkSize2 ;
		if (n2Inner % gpuBlkSize2) GridSize2++ ;
		int GridSize3 = -1 ;
		if (Config::Instance()->cb3 >= n3Inner)
		{
		GridSize3 = 1 ;
		dimz      = n3Inner ;
		}
		else
		{		
		GridSize3 = n3Inner / Config::Instance()->cb3 ;
		if (n3Inner % Config::Instance()->cb3) GridSize3++ ;
		dimz      = Config::Instance()->cb3 ;
		}
		
		dim3 BlkSize(gpuBlkSize1, gpuBlkSize2, gpuBlkSize3) ;		
		dim3 GridSize(GridSize1, GridSize2, GridSize3) ;						

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

			// TODO overrides gpuBlkSize
		int gpuBlkSize1 = Config::Instance()->cb1 ;
		int gpuBlkSize2 = Config::Instance()->cb2 ;				
		int gpuBlkSize3 = 1 ;
		
		int dimx = n1Inner ;
		int dimy = n2Inner ;
		int dimz = -1 ;	
		
		int GridSize1 = n1Inner / gpuBlkSize1 ;
		if (n1Inner % gpuBlkSize1) GridSize1++ ;
		int GridSize2 = n2Inner / gpuBlkSize2 ;
		if (n2Inner % gpuBlkSize2) GridSize2++ ;
		int GridSize3 = -1 ;
		if (Config::Instance()->cb3 >= n3Inner)
		{
		GridSize3 = 1 ;
		dimz      = n3Inner ;
		}
		else
		{		
		GridSize3 = n3Inner / Config::Instance()->cb3 ;
		if (n3Inner % Config::Instance()->cb3) GridSize3++ ;
		dimz      = Config::Instance()->cb3 ;
		}		
		
		dim3 BlkSize(gpuBlkSize1, gpuBlkSize2, gpuBlkSize3) ;		
		dim3 GridSize(GridSize1, GridSize2, GridSize3) ;						
	
	if ((dim == DIM3) && (fdOrder == 8))
	{			   
			kernelOpt_computePressureWithFD_3D_O8<<<GridSize, BlkSize>>>(dim, fdOrder, d_grid_3d, prc_d_grid_3d, coef_d_grid_3d,inv2_d1,inv2_d2,inv2_d3,
					n1,n2,n3,i1Start,i1End,i2Start,i2End,i3Start,i3End, dimx, dimy, dimz);
	}

	cudaCheckError();
	cudaDeviceSynchronize();

	printDebug(FULL_DEBUG, "Out Grid_Cuda_Optim::computePressureWithFD") ;
	return(RTN_CODE_OK) ;
}


} // namespace hpcscan
