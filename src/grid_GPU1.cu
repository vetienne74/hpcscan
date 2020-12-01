
//-------------------------------------------------------------------------------------------------------
// Derived class from Grid
// Optimized for GPU
// Version 1 ??
//-------------------------------------------------------------------------------------------------------

#include "grid_GPU1.h"

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
__global__ void cuda_fill_const(Myfloat *data, Myfloat val, int n1, int n2, int n3, Myint64 i1Start, Myint64 i1End, Myint64 i2Start, Myint64 i2End, Myint64 i3Start, Myint64 i3End)
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
		int t_i3 = tid / (n1*n2);
		int idx = tid-t_i3*n1*n2;
		int t_i2 = idx/n1;
		int t_i1 = idx%n1;

		if (t_i1 >= i1Start && t_i1 <= i1End &&
			t_i2 >= i2Start && t_i2 <= i2End &&
			t_i3 >= i3Start && t_i3 <= i3End   )
		{
			data[tid] = val;
			//printf("data[%d]=%f\n",tid,val);
		}

		tid += blockDim.x * gridDim.x;
	}
}

//-------------------------------------------------------------------------------------------------------

__global__ void cuda_fill_sine(Myfloat *data, Myfloat64 param1, Myfloat64 param2, Myfloat64 param3, Myfloat64 amp, int n1, int n2, int n3, Myint64 i1Start, Myint64 i1End, Myint64 i2Start, Myint64 i2End, Myint64 i3Start, Myint64 i3End, Myfloat Orig1, Myfloat Orig2, Myfloat Orig3, Myfloat64 d1, Myfloat64 d2, Myfloat64 d3 )
{
	// printf("sine %f %f %f %f %f\n",param1,param2,param3,amp,Orig1);

	int size = n1*n2*n3;
	int tid = threadIdx.x + blockIdx.x*blockDim.x;

	while (tid < size)
	{
		int t_i3 = tid / (n1*n2);
		int idx = tid-t_i3*n1*n2;
		int t_i2 = idx/n1;
		int t_i1 = idx%n1;

		if (t_i1 >= i1Start && t_i1 <= i1End &&
			t_i2 >= i2Start && t_i2 <= i2End &&
			t_i3 >= i3Start && t_i3 <= i3End   )
		{
			Myfloat64 coord1 = Myfloat64(Orig1 + t_i1 * d1);
			Myfloat64 coord2 = Myfloat64(Orig2 + t_i2 * d2);
			Myfloat64 coord3 = Myfloat64(Orig3 + t_i3 * d3);

			Myfloat val = amp * sin(coord1 * param1) * sin(coord2 * param2) * sin(coord3 * param3);

			data[tid] = val;
			//printf("data[%d]=%f\n",tid,val);
		}

		tid += blockDim.x * gridDim.x;
	}
}

//-------------------------------------------------------------------------------------------------------

__global__ void cuda_fill_linear(Myfloat *data, Myfloat64 param1, Myfloat64 param2, Myfloat64 param3, Myfloat64 amp, int n1, int n2, int n3, Myint64 i1Start, Myint64 i1End, Myint64 i2Start, Myint64 i2End, Myint64 i3Start, Myint64 i3End, Myfloat Orig1, Myfloat Orig2, Myfloat Orig3, Myfloat64 d1, Myfloat64 d2, Myfloat64 d3)
{
	// printf("linear %f %f %f %f\n",param1,param2,param3,amp);

	int size = n1*n2*n3;
	int tid = threadIdx.x + blockIdx.x*blockDim.x;

	while (tid < size)
	{
		int t_i3 = tid / (n1*n2);
		int idx = tid-t_i3*n1*n2;
		int t_i2 = idx/n1;
		int t_i1 = idx%n1;

		if (t_i1 >= i1Start && t_i1 <= i1End &&
			t_i2 >= i2Start && t_i2 <= i2End &&
			t_i3 >= i3Start && t_i3 <= i3End   )
		{
			Myfloat64 coord1 = Myfloat64(Orig1 + t_i1 * d1);
			Myfloat64 coord2 = Myfloat64(Orig2 + t_i2 * d2);
			Myfloat64 coord3 = Myfloat64(Orig3 + t_i3 * d3);

			Myfloat val = amp * coord1 * coord2 * coord3;

			data[tid] = val;
			//printf("data[%d]=%f\n",tid,val);
		}

		tid += blockDim.x * gridDim.x;
	}
}

//-------------------------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------------------------

__global__ void cuda_diff(Myfloat *data1, Myfloat *data2, Myfloat *dataOut, int n1, int n2, int n3, Myint64 i1Start, Myint64 i1End, Myint64 i2Start, Myint64 i2End, Myint64 i3Start, Myint64 i3End)
{
	int size = n1*n2*n3;
	int tid = threadIdx.x + blockIdx.x*blockDim.x;

	while (tid < size)
	{
		int t_i3 = tid / (n1*n2);
		int idx = tid-t_i3*n1*n2;
		int t_i2 = idx/n1;
		int t_i1 = idx%n1;

		dataOut[tid]=0;

		if (t_i1 >= i1Start && t_i1 <= i1End &&
			t_i2 >= i2Start && t_i2 <= i2End &&
			t_i3 >= i3Start && t_i3 <= i3End   )
		{
			dataOut[tid] = fabsf(data1[tid]-data2[tid]);
		}

		tid += blockDim.x * gridDim.x;
	}
}

//-------------------------------------------------------------------------------------------------------

__global__ void cuda_fabsf(Myfloat *data, Myfloat *dataOut, int n1, int n2, int n3, Myint64 i1Start, Myint64 i1End, Myint64 i2Start, Myint64 i2End, Myint64 i3Start, Myint64 i3End)
{
	int size = n1*n2*n3;
	int tid = threadIdx.x + blockIdx.x*blockDim.x;

	while (tid < size)
	{
		int t_i3 = tid / (n1*n2);
		int idx = tid-t_i3*n1*n2;
		int t_i2 = idx/n1;
		int t_i1 = idx%n1;

		dataOut[tid]=0;

		if (t_i1 >= i1Start && t_i1 <= i1End &&
			t_i2 >= i2Start && t_i2 <= i2End &&
			t_i3 >= i3Start && t_i3 <= i3End   )
		{
			dataOut[tid] = fabsf(data[tid]);
		}

		tid += blockDim.x * gridDim.x;
	}
}

//-------------------------------------------------------------------------------------------------------

__global__ void cuda_mask(Myfloat *data, Myfloat *dataOut, Myfloat val, int n1, int n2, int n3, Myint64 i1Start, Myint64 i1End, Myint64 i2Start, Myint64 i2End, Myint64 i3Start, Myint64 i3End)
{
	int size = n1*n2*n3;
	int tid = threadIdx.x + blockIdx.x*blockDim.x;

	while (tid < size)
	{
		int t_i3 = tid / (n1*n2);
		int idx = tid-t_i3*n1*n2;
		int t_i2 = idx/n1;
		int t_i1 = idx%n1;

		dataOut[tid]=val;

		if (t_i1 >= i1Start && t_i1 <= i1End &&
			t_i2 >= i2Start && t_i2 <= i2End &&
			t_i3 >= i3Start && t_i3 <= i3End   )
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

__global__ void cuda_updatePressure(Myfloat *prn, Myfloat *prc, Myfloat *coef, Myfloat *lapla, int n1, int n2, int n3, Myint64 i1Start, Myint64 i1End, Myint64 i2Start, Myint64 i2End, Myint64 i3Start, Myint64 i3End)
{
	int size = n1*n2*n3;
	int tid = threadIdx.x + blockIdx.x*blockDim.x;

	while (tid < size)
	{
		int t_i3 = tid / (n1*n2);
		int idx = tid-t_i3*n1*n2;
		int t_i2 = idx/n1;
		int t_i1 = idx%n1;

		// dataOut[tid]=val;

		if (t_i1 >= i1Start && t_i1 <= i1End &&
			t_i2 >= i2Start && t_i2 <= i2End &&
			t_i3 >= i3Start && t_i3 <= i3End   )
		{
			prn[tid]=2.0*prc[tid]-prn[tid]+coef[tid]*lapla[tid];
		}

		tid += blockDim.x * gridDim.x;
	}
}

//-------------------------------------------------------------------------------------------------------




















Grid_GPU1::Grid_GPU1(Grid_type gridTypeIn) : Grid(gridTypeIn)
														{
	printDebug(MID_DEBUG, "IN Grid_GPU1::Grid_GPU1");

	gridMode = "GPU1" ;

	printDebug(MID_DEBUG, "OUT Grid_GPU1::Grid_GPU1");
														}

//-------------------------------------------------------------------------------------------------------

Grid_GPU1::Grid_GPU1(Grid_type gridTypeIn, Dim_type dimIn,
		Myint64 n1InnerIn, Myint64 n2InnerIn, Myint64 n3InnerIn) : Grid(gridTypeIn, dimIn,
				n1InnerIn, n2InnerIn, n3InnerIn)
{
	printDebug(MID_DEBUG, "IN Grid_GPU1::Grid_GPU1");

	gridMode = "GPU1" ;

	printDebug(MID_DEBUG, "OUT Grid_GPU1::Grid_GPU1");
}

//-------------------------------------------------------------------------------------------------------

Grid_GPU1::~Grid_GPU1(void)
{
	printDebug(MID_DEBUG, "IN Grid_GPU1::~Grid_GPU1");

	//delete[] grid_3d ;
	cudaFree(d_grid_3d);
	cudaFree(d_help_3d);
	cudaCheckError();

	printDebug(MID_DEBUG, "OUT Grid_GPU1::~Grid_GPU1");
}

//-------------------------------------------------------------------------------------------------------

void Grid_GPU1::info(void)
{
	printDebug(FULL_DEBUG, "IN Grid_GPU1::info");

	// parent class info
	Grid::info() ;

	// additional info
	// TO DO

	printDebug(FULL_DEBUG, "IN Grid_GPU1::info");
}


//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_GPU1::FD_LAPLACIAN(Point_type pointType, const Grid& Wgrid, Myint fdOrder)
{
	printDebug(MID_DEBUG, "IN Grid_GPU1::FD_LAPLACIAN");

	// TO DO
	Grid::FD_LAPLACIAN(pointType, Wgrid, fdOrder) ;

	printDebug(MID_DEBUG, "OUT Grid_GPU1::FD_LAPLACIAN");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_GPU1::computePressureWithFD(Grid& prcGridIn, Grid& coefGridIn, Myint fdOrder)
{

	printDebug(FULL_DEBUG, "In Grid_GPU1::computePressureWithFD") ;

	// TO DO
	Grid::computePressureWithFD(prcGridIn, coefGridIn, fdOrder) ;

	printDebug(FULL_DEBUG, "Out Grid_GPU1::computePressureWithFD") ;
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

void Grid_GPU1::initializeGrid(void)
{
	printDebug(FULL_DEBUG, "In Grid_GPU1::initializeGrid") ;

	Grid::initializeGrid() ; // this sets up halos etc.
	printf("test n1=%d n2=%d n3=%d\n",n1,n2,n3);

	if (d_grid_3d == NULL)
	{
		cudaMalloc( (void**)&d_grid_3d, n1*n2*n3*sizeof(Myfloat) );
		cudaCheckError();

		cudaMalloc( (void**)&d_help_3d, n1*n2*n3*sizeof(Myfloat) );
		cudaCheckError();
	}
	printDebug(FULL_DEBUG, "Out Grid_GPU1::initializeGrid") ;
}

//-------------------------------------------------------------------------------------------------------
void Grid_GPU1::fill(Point_type pointType, Myfloat val)
{
	printDebug(FULL_DEBUG, "In Grid_GPU1::fill") ;

	// Grid::fill(pointType, val) ; // fill CPU memory (remove me)
	
	//pointType
	Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
	Grid::getGridIndex(pointType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End);
	cuda_fill_const<<<1024,128>>>(d_grid_3d,val,n1,n2,n3,i1Start,i1End,i2Start,i2End,i3Start,i3End);
	cudaDeviceSynchronize();
	cudaCheckError();

	printDebug(FULL_DEBUG, "Out Grid_GPU1::fill") ;
}

//-------------------------------------------------------------------------------------------------------
void Grid_GPU1::fill(Point_type pointType, Func_type t1,  Func_type t2, Func_type t3,
		Myfloat64 param1, Myfloat64 param2, Myfloat64 param3, Myfloat64 amp)
{
	printDebug(FULL_DEBUG, "In Grid_GPU1::fill") ;

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

	if (t1 == FUNC_SINE) cuda_fill_sine<<<1024,128>>>  (d_grid_3d,param1,param2,param3,amp,n1,n2,n3,i1Start,i1End,i2Start,i2End,i3Start,i3End,Orig1,Orig2,Orig3,d1,d2,d2);
	else                 cuda_fill_linear<<<1024,128>>>(d_grid_3d,param1,param2,param3,amp,n1,n2,n3,i1Start,i1End,i2Start,i2End,i3Start,i3End,Orig1,Orig2,Orig3,d1,d2,d2);


	printDebug(FULL_DEBUG, "Out Grid_GPU1::fill") ;
}

//-------------------------------------------------------------------------------------------------------
Myfloat Grid_GPU1::getMin(Point_type pointType)
{
	printDebug(FULL_DEBUG, "In Grid_GPU1::getMin") ;

	// TO DO
	// return(Grid::getMin(pointType)) ;

	//pointType
	Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
	Grid::getGridIndex(pointType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End);

	cuda_mask<<<1024,256>>>(d_grid_3d,d_help_3d,999,n1,n2,n3,i1Start,i1End,i2Start,i2End,i3Start,i3End);

	thrust::device_ptr<Myfloat> d_help_3d_ptr = thrust::device_pointer_cast(d_help_3d);
	thrust::device_ptr<hpcscan::Myfloat> vptr = thrust::min_element(thrust::device, d_help_3d_ptr, d_help_3d_ptr + n1*n2*n3);
	float val = *vptr;
	// printf("val %f\n",val);
	return val;

	printDebug(FULL_DEBUG, "Out Grid_GPU1::getMin") ;
}

//-------------------------------------------------------------------------------------------------------
Myfloat Grid_GPU1::getMax(Point_type pointType)
{
	printDebug(FULL_DEBUG, "In Grid_GPU1::getMax") ;

	// TO DO
	// return(Grid::getMax(pointType)) ;

	//pointType
	Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
	Grid::getGridIndex(pointType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End);

	// cuda_mask<<<1024,256>>>(d_grid_3d,d_help_3d,0,n1,n2,n3,i1Start,i1End,i2Start,i2End,i3Start,i3End);

	thrust::device_ptr<Myfloat> d_help_3d_ptr = thrust::device_pointer_cast(d_grid_3d);
	thrust::device_ptr<hpcscan::Myfloat> vptr = thrust::max_element(thrust::device, d_help_3d_ptr, d_help_3d_ptr + n1*n2*n3);
	float val = *vptr;
	// printf("val %f\n",val);
	return val;

	printDebug(FULL_DEBUG, "Out Grid_GPU1::getMax") ;
}

//-------------------------------------------------------------------------------------------------------
Myfloat Grid_GPU1::L1Err(Point_type pointType, const Grid& gridIn) const
{
	printDebug(FULL_DEBUG, "In Grid_GPU1::L1Err") ;

	// TO DO
	// return(Grid::L1Err(pointType, gridIn)) ;

	//pointType
	Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
	Grid::getGridIndex(pointType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End);

	thrust::device_ptr<Myfloat> d_help_3d_ptr;

	cuda_diff<<<1024,256>>>(d_grid_3d,gridIn.d_grid_3d,d_help_3d,n1,n2,n3,i1Start,i1End,i2Start,i2End,i3Start,i3End);
	d_help_3d_ptr = thrust::device_pointer_cast(d_help_3d);
	double totErr = thrust::reduce(thrust::device, d_help_3d_ptr, d_help_3d_ptr + n1*n2*n3);

	double totArr;
	if (false)
	{
		cuda_fabsf<<<1024,256>>>(gridIn.d_grid_3d,d_help_3d,n1,n2,n3,i1Start,i1End,i2Start,i2End,i3Start,i3End);
		totArr = thrust::reduce(thrust::device, d_help_3d_ptr, d_help_3d_ptr + n1*n2*n3);
	}
	else // assuming grid values are positive
	{
		d_help_3d_ptr = thrust::device_pointer_cast(gridIn.d_grid_3d);
		totArr = thrust::reduce(thrust::device, d_help_3d_ptr, d_help_3d_ptr + n1*n2*n3);
	}

	cudaDeviceSynchronize();

	if (totArr < MAX_ERR_FLOAT) totArr = 1.0 * npoint ;

	return totErr/totArr;

	

	printDebug(FULL_DEBUG, "Out Grid_GPU1::L1Err") ;
}
//-------------------------------------------------------------------------------------------------------
Rtn_code Grid_GPU1::updatePressure(Point_type pointType, const Grid& prcGrid,
		const Grid& coefGrid, const Grid& laplaGrid)
{
	printDebug(FULL_DEBUG, "In Grid_GPU1::updatePressure") ;

	// TO DO
	// return(Grid::updatePressure(pointType, prcGrid, coefGrid, laplaGrid)) ;

	//pointType
	Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
	Grid::getGridIndex(pointType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End);

	cuda_updatePressure<<<1024,256>>>(d_grid_3d, prcGrid.d_grid_3d, coefGrid.d_grid_3d, laplaGrid.d_grid_3d, n1, n2, n3, i1Start, i1End, i2Start, i2End, i3Start, i3End);

	cudaDeviceSynchronize();

	printDebug(FULL_DEBUG, "Out Grid_GPU1::updatePressure") ;
}
//-------------------------------------------------------------------------------------------------------
Rtn_code Grid_GPU1::applyBoundaryCondition(BoundCond_type boundCondType)
{
printDebug(FULL_DEBUG, "In Grid_GPU1::applyBoundaryCondition") ;

// TO DO
printf("BNDNDNND\n");
return(Grid::applyBoundaryCondition(boundCondType)) ;

// cuda_updatePressure<<<1024,256>>>(d_grid_3d, prcGrid.d_grid_3d, coefGrid.d_grid_3d, laplaGrid.d_grid_3d, n1, n2, n3, i1Start, i1End, i2Start, i2End, i3Start, i3End);

cudaDeviceSynchronize();

printDebug(FULL_DEBUG, "Out Grid_GPU1::applyBoundaryCondition") ;
}
} // namespace hpcscan
