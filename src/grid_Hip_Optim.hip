#include "hip/hip_runtime.h"

//-------------------------------------------------------------------------------------------------------
// This grid is activated with command line option -testMode HIP_Optim
// Derived class from Grid_Hip
// HIP implementation (target GPU)
//-------------------------------------------------------------------------------------------------------

#include "grid_Hip_Optim.h"

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

namespace hpcscan {

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
// compute derivative along axis 1
// input u
// output w

__global__ void kernelOpt_FD_D2_N1_O8(const Myint fdOrder, Myfloat *w, Myfloat *u,
		const Myfloat inv2_d1, const Myfloat inv2_d2, const Myfloat inv2_d3,
		const Myint n1, const Myint n2, const Myint n3,
		const Myint64 i1Start, const Myint64 i1End, const Myint64 i2Start, const Myint64 i2End, const Myint64 i3Start, const Myint64 i3End)
{
	unsigned int i1 = threadIdx.x + blockIdx.x * blockDim.x ;
	unsigned int i2 = threadIdx.y + blockIdx.y * blockDim.y ;
	unsigned int i3 = threadIdx.z + blockIdx.z * blockDim.z ;

	if (i1 >= i1Start && i1 <= i1End &&
			i2 >= i2Start && i2 <= i2End &&
			i3 >= i3Start && i3 <= i3End   )
	{
		w[i1+i2*n1+i3*n1*n2] =
				FD_D2_O8_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
	}
}

//-------------------------------------------------------------------------------------------------------
// compute derivative along axis 2
// input u
// output w

__global__ void kernelOpt_FD_D2_N2_O8(const Myint fdOrder, Myfloat *w, Myfloat *u,
		const Myfloat inv2_d1, const Myfloat inv2_d2, const Myfloat inv2_d3,
		const Myint n1, const Myint n2, const Myint n3,
		const Myint64 i1Start, const Myint64 i1End, const Myint64 i2Start, const Myint64 i2End, const Myint64 i3Start, const Myint64 i3End)
{
	unsigned int i1 = threadIdx.x + blockIdx.x * blockDim.x ;
	unsigned int i2 = threadIdx.y + blockIdx.y * blockDim.y ;
	unsigned int i3 = threadIdx.z + blockIdx.z * blockDim.z ;

	if (i1 >= i1Start && i1 <= i1End &&
			i2 >= i2Start && i2 <= i2End &&
			i3 >= i3Start && i3 <= i3End   )
	{
		w[i1+i2*n1+i3*n1*n2] =
				FD_D2_O8_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
	}
}

//-------------------------------------------------------------------------------------------------------
// compute Laplacian
// input u
// output w

__global__ void kernelOpt_FD_LAPLACIAN_O8(const Myint fdOrder, Myfloat *w, Myfloat *u,
		const Myfloat inv2_d1, const Myfloat inv2_d2, const Myfloat inv2_d3,
		const Myint n1, const Myint n2, const Myint n3,
		const Myint64 i1Start, const Myint64 i1End, const Myint64 i2Start, const Myint64 i2End, const Myint64 i3Start, const Myint64 i3End)
{
	unsigned int i1 = threadIdx.x + blockIdx.x * blockDim.x ;
	unsigned int i2 = threadIdx.y + blockIdx.y * blockDim.y ;
	unsigned int i3 = threadIdx.z + blockIdx.z * blockDim.z ;

	if (i1 >= i1Start && i1 <= i1End &&
			i2 >= i2Start && i2 <= i2End &&
			i3 >= i3Start && i3 <= i3End   )
	{
		w[i1+i2*n1+i3*n1*n2] =
				FD_D2_O8_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
				+ FD_D2_O8_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
				+ FD_D2_O8_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
	}
}
//-------------------------------------------------------------------------------------------------------
// compute derivative along axis 3
// input u
// output w

__global__ void kernelOpt_FD_D2_N3_O8(const Myint fdOrder, Myfloat *w, Myfloat *u,
		const Myfloat inv2_d1, const Myfloat inv2_d2, const Myfloat inv2_d3,
		const Myint n1, const Myint n2, const Myint n3,
		const Myint64 i1Start, const Myint64 i1End, const Myint64 i2Start, const Myint64 i2End, const Myint64 i3Start, const Myint64 i3End)
{
	unsigned int i1 = threadIdx.x + blockIdx.x * blockDim.x ;
	unsigned int i2 = threadIdx.y + blockIdx.y * blockDim.y ;
	unsigned int i3 = threadIdx.z + blockIdx.z * blockDim.z ;

	if (i1 >= i1Start && i1 <= i1End &&
			i2 >= i2Start && i2 <= i2End &&
			i3 >= i3Start && i3 <= i3End   )
	{
		w[i1+i2*n1+i3*n1*n2] =
				FD_D2_O8_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
	}
}

//-------------------------------------------------------------------------------------------------------

Grid_Hip_Optim::Grid_Hip_Optim(Grid_type gridTypeIn) : Grid_Hip(gridTypeIn)
{
	printDebug(MID_DEBUG, "IN Grid_Hip_Optim::Grid_Hip_Optim");

	gridMode = GRID_MODE_HIP_OPTIM ;

	gpuBlkSize1  = Config::Instance()->gpuBlkSize1 ;
	gpuBlkSize2  = Config::Instance()->gpuBlkSize2 ;
	gpuBlkSize3  = Config::Instance()->gpuBlkSize3 ;

	printDebug(MID_DEBUG, "OUT Grid_Hip_Optim::Grid_Hip_Optim");
}

//-------------------------------------------------------------------------------------------------------

Grid_Hip_Optim::Grid_Hip_Optim(Grid_type gridTypeIn, Dim_type dimIn,
		Myint64 n1InnerIn, Myint64 n2InnerIn, Myint64 n3InnerIn) : Grid_Hip(gridTypeIn, dimIn,
				n1InnerIn, n2InnerIn, n3InnerIn)
{
	printDebug(MID_DEBUG, "IN Grid_Hip_Optim::Grid_Hip_Optim");

	gridMode = GRID_MODE_HIP_OPTIM ;

	gpuBlkSize1  = Config::Instance()->gpuBlkSize1 ;
	gpuBlkSize2  = Config::Instance()->gpuBlkSize2 ;
	gpuBlkSize3  = Config::Instance()->gpuBlkSize3 ;

	printDebug(MID_DEBUG, "OUT Grid_Hip_Optim::Grid_Hip_Optim");
}

//-------------------------------------------------------------------------------------------------------

Grid_Hip_Optim::~Grid_Hip_Optim(void)
{
	printDebug(MID_DEBUG, "IN Grid_Hip_Optim::~Grid_Hip_Optim");

	Grid_Hip::~Grid_Hip() ;

	printDebug(MID_DEBUG, "OUT Grid_Hip_Optim::~Grid_Hip_Optim");
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_Hip_Optim::initializeGrid(void)
{
	printDebug(FULL_DEBUG, "In Grid_Hip_Optim::initializeGrid") ;

	Grid_Hip::initializeGrid() ;

	printDebug(FULL_DEBUG, "Out Grid_Hip_Optim::initializeGrid") ;
	return(RTN_CODE_OK) ;
}


//-------------------------------------------------------------------------------------------------------

void Grid_Hip_Optim::info(void)
{
	printDebug(FULL_DEBUG, "IN Grid_Hip_Optim::info");

	// parent class info
	Grid_Hip::info() ;

	// additional info
	printInfo(MASTER, " Block size 1", gpuBlkSize1) ;
	printInfo(MASTER, " Block size 2", gpuBlkSize2) ;
	printInfo(MASTER, " Block size 3", gpuBlkSize3) ;

	printDebug(FULL_DEBUG, "OUT Grid_Hip_Optim::info");
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_Hip_Optim::FD_D2_N1(Point_type pointType, const Grid& Wgrid, Myint fdOrder)
{
	printDebug(MID_DEBUG, "IN Grid_Hip_Optim::FD_D2_N1");

	// TO DO
	Grid_Hip::FD_D2_N1(pointType, Wgrid, fdOrder) ;
	// check grids are same size
	if (this->sameSize(Wgrid) != true)
	{
		printError("Grid_Hip_Optim::FD_D2_N1, grids have not same size") ;
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

	Myfloat * d_w = ((Grid_Hip_Optim&) Wgrid).d_grid_3d ;
	Myfloat * d_u = this->d_grid_3d ;

	if (fdOrder == 8)
	{
		int GridSize1 = n1 / gpuBlkSize1 + 1 ;
		int GridSize2 = n2 / gpuBlkSize2 + 1 ;
		int GridSize3 = n3 / gpuBlkSize3 + 1 ;
		dim3 BlkSize(gpuBlkSize1, gpuBlkSize2, gpuBlkSize3) ;
		dim3 GridSize(GridSize1, GridSize2, GridSize3) ;
		hipLaunchKernelGGL(kernelOpt_FD_D2_N1_O8, dim3(GridSize), dim3(BlkSize), 0, 0, fdOrder, d_w, d_u,inv2_d1,inv2_d2,inv2_d3,
				n1,n2,n3,i1Start,i1End,i2Start,i2End,i3Start,i3End);

		cudaCheckError();
		hipDeviceSynchronize();
	}
	else
	{
		Grid_Hip::FD_D2_N1(pointType, Wgrid, fdOrder) ;
	}

	printDebug(MID_DEBUG, "OUT Grid_Hip_Optim::FD_D2_N1");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_Hip_Optim::FD_D2_N2(Point_type pointType, const Grid& Wgrid, Myint fdOrder)
{
	printDebug(MID_DEBUG, "IN Grid_Hip_Optim::FD_D2_N2");

	// check grids are same size
	if (this->sameSize(Wgrid) != true)
	{
		printError("Grid_Hip_Optim::FD_D2_N2, grids have not same size") ;
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

	Myfloat * d_w = ((Grid_Hip_Optim&) Wgrid).d_grid_3d ;
	Myfloat * d_u = this->d_grid_3d ;

	if (fdOrder == 8)
	{
		int GridSize1 = n1 / gpuBlkSize1 + 1 ;
		int GridSize2 = n2 / gpuBlkSize2 + 1 ;
		int GridSize3 = n3 / gpuBlkSize3 + 1 ;
		dim3 BlkSize(gpuBlkSize1, gpuBlkSize2, gpuBlkSize3) ;
		dim3 GridSize(GridSize1, GridSize2, GridSize3) ;
		hipLaunchKernelGGL(kernelOpt_FD_D2_N2_O8, dim3(GridSize), dim3(BlkSize), 0, 0, fdOrder, d_w, d_u,inv2_d1,inv2_d2,inv2_d3,
				n1,n2,n3,i1Start,i1End,i2Start,i2End,i3Start,i3End);

		cudaCheckError();
		hipDeviceSynchronize();
	}
	else
	{
		Grid_Hip::FD_D2_N2(pointType, Wgrid, fdOrder) ;
	}

	printDebug(MID_DEBUG, "OUT Grid_Hip_Optim::FD_D2_N2");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_Hip_Optim::FD_D2_N3(Point_type pointType, const Grid& Wgrid, Myint fdOrder)
{
	printDebug(MID_DEBUG, "IN Grid_Hip_Optim::FD_D2_N3");

	// check grids are same size
	if (this->sameSize(Wgrid) != true)
	{
		printError("Grid_Hip_Optim::FD_D2_N3, grids have not same size") ;
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

	Myfloat * d_w = ((Grid_Hip_Optim&) Wgrid).d_grid_3d ;
	Myfloat * d_u = this->d_grid_3d ;

	if (fdOrder == 8)
	{
		int GridSize1 = n1 / gpuBlkSize1 + 1 ;
		int GridSize2 = n2 / gpuBlkSize2 + 1 ;
		int GridSize3 = n3 / gpuBlkSize3 + 1 ;
		dim3 BlkSize(gpuBlkSize1, gpuBlkSize2, gpuBlkSize3) ;
		dim3 GridSize(GridSize1, GridSize2, GridSize3) ;
		hipLaunchKernelGGL(kernelOpt_FD_D2_N3_O8, dim3(GridSize), dim3(BlkSize), 0, 0, fdOrder, d_w, d_u,inv2_d1,inv2_d2,inv2_d3,
				n1,n2,n3,i1Start,i1End,i2Start,i2End,i3Start,i3End);

		cudaCheckError();
		hipDeviceSynchronize();
	}
	else
	{
		Grid_Hip::FD_D2_N3(pointType, Wgrid, fdOrder) ;
	}

	printDebug(MID_DEBUG, "OUT Grid_Hip_Optim::FD_D2_N3");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_Hip_Optim::FD_LAPLACIAN(Point_type pointType, const Grid& Wgrid, Myint fdOrder)
{
	printDebug(MID_DEBUG, "IN Grid_Hip_Optim::FD_LAPLACIAN");

	// check grids are same size
	if (this->sameSize(Wgrid) != true)
	{
		printError("Grid_Hip_Optim::FD_LAPLACIAN, grids have not same size") ;
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

	Myfloat * d_w = ((Grid_Hip_Optim&) Wgrid).d_grid_3d ;
	Myfloat * d_u = this->d_grid_3d ;

	if ((dim == DIM3) && (fdOrder == 8))
	{
		int GridSize1 = n1 / gpuBlkSize1 + 1 ;
		int GridSize2 = n2 / gpuBlkSize2 + 1 ;
		int GridSize3 = n3 / gpuBlkSize3 + 1 ;
		dim3 BlkSize(gpuBlkSize1, gpuBlkSize2, gpuBlkSize3) ;
		dim3 GridSize(GridSize1, GridSize2, GridSize3) ;
		hipLaunchKernelGGL(kernelOpt_FD_LAPLACIAN_O8, dim3(GridSize), dim3(BlkSize), 0, 0, fdOrder, d_w, d_u,inv2_d1,inv2_d2,inv2_d3,
				n1,n2,n3,i1Start,i1End,i2Start,i2End,i3Start,i3End);

		cudaCheckError();
		hipDeviceSynchronize();
	}
	else
	{
		Grid_Hip::FD_LAPLACIAN(pointType, Wgrid, fdOrder) ;
	}

	printDebug(MID_DEBUG, "OUT Grid_Hip_Optim::FD_LAPLACIAN");
	return(RTN_CODE_OK) ;
}

} // namespace hpcscan
