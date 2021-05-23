
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

#include <CL/sycl.hpp>
#include "mpi.h"

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

	d_grid_3d = NULL;
	d_help_3d = NULL;

	printDebug(MID_DEBUG, "OUT Grid_DPCPP::Grid_DPCPP");
														}

//-------------------------------------------------------------------------------------------------------

Grid_DPCPP::Grid_DPCPP(Grid_type gridTypeIn, Dim_type dimIn,
		Myint64 n1InnerIn, Myint64 n2InnerIn, Myint64 n3InnerIn) : Grid(gridTypeIn, dimIn,
				n1InnerIn, n2InnerIn, n3InnerIn)
{
	printDebug(MID_DEBUG, "IN Grid_DPCPP::Grid_DPCPP");

	gridMode = GRID_MODE_DPCPP ;

	d_grid_3d = NULL;
	d_help_3d = NULL;

	printDebug(MID_DEBUG, "OUT Grid_DPCPP::Grid_DPCPP");
}

//-------------------------------------------------------------------------------------------------------

Grid_DPCPP::~Grid_DPCPP(void)
{
	printDebug(MID_DEBUG, "IN Grid_DPCPP::~Grid_DPCPP");

	//delete[] grid_3d ;
	//free(d_grid_3d);
	//free(d_help_3d);
	//free(d_help_3d_2);

	printDebug(MID_DEBUG, "OUT Grid_DPCPP::~Grid_DPCPP");
}


//-------------------------------------------------------------------------------------------------------

void Grid_DPCPP::info(void)
{
	printDebug(FULL_DEBUG, "IN Grid_DPCPP::info");

	// parent class info
	Grid::info() ;

	// additional info
	printDebug(FULL_DEBUG, "IN Grid_DPCPP::info");
}

//-------------------------------------------------------------------------------------------------------

void Grid_DPCPP::fillArray(Myfloat val)
{
	printDebug(MID_DEBUG, "IN Grid_DPCPP::fillArray");

	Myfloat * const w = grid_3d ;

//#pragma omp parallel for
//	for (Myint64 ii=0; ii<npoint; ii++)
//	{
//		w[ii] = val ;
//	}

	sycl::queue Q ;

	Myfloat* w_d = sycl::malloc_shared<Myfloat>(npoint, Q) ;

	Q.parallel_for(npoint,[=](int ii)
	{
		w_d[ii] = val ;
	}).wait() ;

	Q.memcpy(w, w_d, npoint*sizeof(Myfloat)).wait() ;

	//for (int i=0; i<npoint; i++)
	//{
	//	w[i] = w_d[i] ;
	//}

	printDebug(MID_DEBUG, "OUT Grid_DPCPP::fillArray");
}

} // namespace hpcscan
