
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

__global__ void print_from_gpu(void) 
{
    printf("Hello World! from thread [%d,%d] \
        From device\n", threadIdx.x,blockIdx.x);
}

//-------------------------------------------------------------------------------------------------------

// first (wrong) implementation of filling data(0:n1*n2*n3)=val
// TODO we shouldn't ignore pointType
// we could replace it with a thrust:: one-liner
__global__ void cuda_fill_const(Myfloat *data, Myfloat val, int n1, int n2, int n3, Myint64 i1Start, Myint64 i1End, Myint64 i2Start, Myint64 i2End, Myint64 i3Start, Myint64 i3End)
{
        Myint64 s1 = i1End - i1Start + 1;
        Myint64 s2 = i2End - i2Start + 1;
        Myint64 s3 = i3End - i3Start + 1;

        int size = n1*n2*n3;
        int size2= s1*s2*s3;

        printf("size1:%d\n",size);
        printf("size2:%d\n",size2);
	int tid = threadIdx.x + blockIdx.x*blockDim.x;

	while (tid < size) 
    {
		data[tid] = val;
		//printf("data[%d]=%f\n",tid,val);
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

Rtn_code Grid_GPU1::FD_LAPLACIAN(Point_type pType, const Grid& Wgrid, Myint fdOrder)
{
	printDebug(MID_DEBUG, "IN Grid_GPU1::FD_LAPLACIAN");

	// TO DO
	Grid::FD_LAPLACIAN(pType, Wgrid, fdOrder) ;

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

	// TO DO
	print_from_gpu<<<4,1>>>();
	print_from_gpu<<<4,1>>>();
	cudaDeviceSynchronize();

	cudaCheckError();

	printf("test\n");

	Grid::initializeGrid() ; // this sets up halos etc.
	printf("test n1=%d n2=%d n3=%d\n",n1,n2,n3);

	cudaMalloc( (void**)&d_grid_3d, n1*n2*n3*sizeof(Myfloat) );
	cudaCheckError();

	printDebug(FULL_DEBUG, "Out Grid_GPU1::initializeGrid") ;
}

//-------------------------------------------------------------------------------------------------------
void Grid_GPU1::fill(Point_type pointType, Myfloat val)
{
	printDebug(FULL_DEBUG, "In Grid_GPU1::fill") ;

	Grid::fill(pointType, val) ; // fill CPU memory (remove me)
        //pointtype
        Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
        Grid::getGridIndex(INNER_POINTS, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End);
        printf("test i1=%d i2=%d i3=%d\n",i1End,i2End,i3End);
        cuda_fill_const<<<512,64>>>(d_grid_3d,val,n1,n2,n3,i1Start, i1End, i2Start, i2End, i3Start, i3End);
	cudaCheckError();

	printDebug(FULL_DEBUG, "Out Grid_GPU1::fill") ;
}

//-------------------------------------------------------------------------------------------------------
void Grid_GPU1::fill(Point_type pType, Func_type t1,  Func_type t2, Func_type t3,
		Myfloat64 param1, Myfloat64 param2, Myfloat64 param3, Myfloat64 amp)
{
	printDebug(FULL_DEBUG, "In Grid_GPU1::fill") ;

	// TO DO
	Grid::fill(pType, t1,  t2, t3, param1, param2, param3, amp) ;

	printDebug(FULL_DEBUG, "Out Grid_GPU1::fill") ;
}

//-------------------------------------------------------------------------------------------------------
Myfloat Grid_GPU1::getMin(Point_type pType)
{
	printDebug(FULL_DEBUG, "In Grid_GPU1::getMin") ;

	// TO DO
	return(Grid::getMin(pType)) ;

	printDebug(FULL_DEBUG, "Out Grid_GPU1::getMin") ;
}

//-------------------------------------------------------------------------------------------------------
Myfloat Grid_GPU1::getMax(Point_type pType)
{
	printDebug(FULL_DEBUG, "In Grid_GPU1::getMax") ;

	// TO DO
	return(Grid::getMax(pType)) ;

	printDebug(FULL_DEBUG, "Out Grid_GPU1::getMax") ;
}

//-------------------------------------------------------------------------------------------------------
Myfloat Grid_GPU1::L1Err(Point_type pointType, const Grid& gridIn) const
{
	printDebug(FULL_DEBUG, "In Grid_GPU1::L1Err") ;

	// TO DO
	return(Grid::L1Err(pointType, gridIn)) ;

	printDebug(FULL_DEBUG, "Out Grid_GPU1::L1Err") ;
}
//-------------------------------------------------------------------------------------------------------
Rtn_code Grid_GPU1::updatePressure(Point_type pType, const Grid& prcGrid,
		const Grid& coefGrid, const Grid& laplaGrid)
{
	printDebug(FULL_DEBUG, "In Grid_GPU1::updatePressure") ;

	// TO DO
	return(Grid::updatePressure(pType, prcGrid, coefGrid, laplaGrid)) ;

	printDebug(FULL_DEBUG, "Out Grid_GPU1::updatePressure") ;
}

} // namespace hpcscan
