
//-------------------------------------------------------------------------------------------------------
// This grid is activated with command line option -testMode DPC++
// Derived class from Grid
// DPC++ implementation (target CPU, GPU & FPGA)
//-------------------------------------------------------------------------------------------------------

#include "grid_DPCPP.h"

#include <algorithm> // for min and max
#include <cassert>
#include <cfloat>  // for FLT_MAX ;
#include <cmath>   // for fabs
#include <cstddef> // for NULL
#include <fstream>
#include <stdio.h>

#include "config.h"
#include "constant.h"
#include "fdm.h"
#include "global.h"
#include "output_report.h"

using namespace std;

namespace hpcscan {

//-------------------------------------------------------------------------------------------------------

Grid_DPCPP::Grid_DPCPP(Grid_type gridTypeIn) : Grid(gridTypeIn)
														{
	printDebug(MID_DEBUG, "IN Grid_DPCPP::Grid_DPCPP");

	gridMode = GRID_MODE_DPCPP ;

	d_grid_3d   = NULL ;
	d_help_3d   = NULL ;
	d_help_3d_2 = NULL ;

	printDebug(MID_DEBUG, "OUT Grid_DPCPP::Grid_DPCPP");
														}

//-------------------------------------------------------------------------------------------------------

Grid_DPCPP::Grid_DPCPP(Grid_type gridTypeIn, Dim_type dimIn,
		Myint64 n1InnerIn, Myint64 n2InnerIn, Myint64 n3InnerIn) : Grid(gridTypeIn, dimIn,
				n1InnerIn, n2InnerIn, n3InnerIn)
{
	printDebug(MID_DEBUG, "IN Grid_DPCPP::Grid_DPCPP");

	gridMode = GRID_MODE_DPCPP ;

	d_grid_3d   = NULL;
	d_help_3d   = NULL;
	d_help_3d_2 = NULL ;

	printDebug(MID_DEBUG, "OUT Grid_DPCPP::Grid_DPCPP");
}

//-------------------------------------------------------------------------------------------------------

Grid_DPCPP::~Grid_DPCPP(void)
{
	printDebug(MID_DEBUG, "IN Grid_DPCPP::~Grid_DPCPP");

	//delete[] grid_3d ;
	sycl::free(d_grid_3d, Q);
	//free(d_help_3d);
	//free(d_help_3d_2);

	printDebug(MID_DEBUG, "OUT Grid_DPCPP::~Grid_DPCPP");
}

//-------------------------------------------------------------------------------------------------------

void Grid_DPCPP::initializeGrid(void)
{
	printDebug(FULL_DEBUG, "In Grid_DPCPP::initializeGrid") ;

	Grid::initializeGrid() ;

	// initialize device queue
	//*Q = sycl::queue() ;

	if (d_grid_3d == NULL)
	{
		// allocate the grid on the device
		d_grid_3d = sycl::malloc_shared<Myfloat>(npoint, Q) ;

		// allocate 1d array of the device used to perform reduction operation
		//cudaMalloc( (void**)&d_help_3d, (gpuGridSize) * sizeof(Myfloat) );
		//cudaCheckError();

		// allocate 1d array of the device used to perform reduction operation
		//cudaMalloc( (void**)&d_help_3d_2, (gpuGridSize) * sizeof(Myfloat) );
		//cudaCheckError();
	}
	printDebug(FULL_DEBUG, "Out Grid_DPCPP::initializeGrid") ;
}

//-------------------------------------------------------------------------------------------------------

void Grid_DPCPP::info(void)
{
	printDebug(FULL_DEBUG, "IN Grid_DPCPP::info");

	// parent class info
	Grid::info() ;

	// additional info
	printInfo(MASTER, "") ;
	printInfo(MASTER, " * Device parameters * ") ;

	printInfo(MASTER, " Selected device", Q.get_device().get_info<sycl::info::device::name>() ) ;
	printInfo(MASTER, " Device vendor\t", Q.get_device().get_info<sycl::info::device::vendor>() ) ;

	printDebug(FULL_DEBUG, "OUT Grid_DPCPP::info");
}

//-------------------------------------------------------------------------------------------------------

void Grid_DPCPP::fillArray(Myfloat val)
{
	printDebug(MID_DEBUG, "IN Grid_DPCPP::fillArray");

	Myfloat * const w = grid_3d ;
	Myfloat* w_d = d_grid_3d ;

	Q.parallel_for(npoint,[=](int ii)
	{
		w_d[ii] = val ;
	}).wait() ;

	Q.memcpy(w, w_d, npoint*sizeof(Myfloat)).wait() ;

	printDebug(MID_DEBUG, "OUT Grid_DPCPP::fillArray");
}

} // namespace hpcscan
