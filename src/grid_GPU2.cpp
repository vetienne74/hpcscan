
//-------------------------------------------------------------------------------------------------------
// Derived class from Grid
// Optimized for GPU
// Version 2 ??
//-------------------------------------------------------------------------------------------------------

#include "grid_GPU2.h"

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

Grid_GPU2::Grid_GPU2(Grid_type gridTypeIn) : Grid(gridTypeIn)
														{
	printDebug(MID_DEBUG, "IN Grid_GPU2::Grid_GPU2");

	gridMode = "GPU2" ;

	printDebug(MID_DEBUG, "OUT Grid_GPU2::Grid_GPU2");
														}

//-------------------------------------------------------------------------------------------------------

Grid_GPU2::Grid_GPU2(Grid_type gridTypeIn, Dim_type dimIn,
		Myint64 n1InnerIn, Myint64 n2InnerIn, Myint64 n3InnerIn) : Grid(gridTypeIn, dimIn,
				n1InnerIn, n2InnerIn, n3InnerIn)
{
	printDebug(MID_DEBUG, "IN Grid_GPU2::Grid_GPU2");

	gridMode = "GPU2" ;

	printDebug(MID_DEBUG, "OUT Grid_GPU2::Grid_GPU2");
}

//-------------------------------------------------------------------------------------------------------

Grid_GPU2::~Grid_GPU2(void)
{
	printDebug(MID_DEBUG, "IN Grid_GPU2::~Grid_GPU2");

	//delete[] grid_3d ;

	printDebug(MID_DEBUG, "OUT Grid_GPU2::~Grid_GPU2");
}

//-------------------------------------------------------------------------------------------------------

void Grid_GPU2::info(void)
{
	printDebug(FULL_DEBUG, "IN Grid_GPU2::info");

	// parent class info
	Grid::info() ;

	// additional info
	// TO DO

	printDebug(FULL_DEBUG, "IN Grid_GPU2::info");
}


//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_GPU2::FD_LAPLACIAN(Point_type pType, const Grid& Wgrid, Myint fdOrder)
{
	printDebug(MID_DEBUG, "IN Grid_GPU2::FD_LAPLACIAN");

	// TO DO
	Grid::FD_LAPLACIAN(pType, Wgrid, fdOrder) ;

	printDebug(MID_DEBUG, "OUT Grid_GPU2::FD_LAPLACIAN");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_GPU2::computePressureWithFD(Grid& prcGridIn, Grid& coefGridIn, Myint fdOrder)
{

	printDebug(FULL_DEBUG, "In Grid_GPU2::computePressureWithFD") ;

	// TO DO
	Grid::computePressureWithFD(prcGridIn, coefGridIn, fdOrder) ;

	printDebug(FULL_DEBUG, "Out Grid_GPU2::computePressureWithFD") ;
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

void Grid_GPU2::initializeGrid(void)
{
	printDebug(FULL_DEBUG, "In Grid_GPU2::initializeGrid") ;

	// TO DO
	Grid::initializeGrid() ;

	printDebug(FULL_DEBUG, "Out Grid_GPU2::initializeGrid") ;
}

//-------------------------------------------------------------------------------------------------------
void Grid_GPU2::fill(Point_type pointType, Myfloat val)
{
	printDebug(FULL_DEBUG, "In Grid_GPU2::fill") ;
 
        Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
        getGridIndex(pointType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End) ;
        Myfloat *grid_3d_2 = grid_3d;         
//#pragma acc declare create(grid_3d_2[0:n1*n2*n2])
//#pragma acc data copy(grid_3d_2[0:n1*n2*n3])
#pragma acc data create(grid_3d_2[0:n1*n2*n3])
#pragma acc parallel loop collapse(2)// present(grid_3d_2)
        for (Myint64 i3 = i3Start; i3<= i3End; i3++)
        {       
                for (Myint64 i2 = i2Start; i2<= i2End; i2++)
               {
                      for (Myint64 i1 = i1Start; i1<= i1End; i1++)
                      {       
                              Myint64 ii = i1 + i2*n1 + i3*n2*n1 ;
                              grid_3d_2[ii] = val ;
                      }
                }
        }
//#pragma acc data copyout(grid_3d_2[0:n1*n2*n3])   
	printDebug(FULL_DEBUG, "Out Grid_GPU2::fill") ;
}

//-------------------------------------------------------------------------------------------------------
void Grid_GPU2::fill(Point_type pType, Func_type t1,  Func_type t2, Func_type t3,
		Myfloat64 param1, Myfloat64 param2, Myfloat64 param3, Myfloat64 amp)
{
	printDebug(FULL_DEBUG, "In Grid_GPU2::fill") ;

	// TO DO
	Grid::fill(pType, t1,  t2, t3, param1, param2, param3, amp) ;

	printDebug(FULL_DEBUG, "Out Grid_GPU2::fill") ;
}

//-------------------------------------------------------------------------------------------------------
Myfloat Grid_GPU2::getMin(Point_type pType)
{
	printDebug(FULL_DEBUG, "In Grid_GPU2::getMin") ;

	// TO DO
	return(Grid::getMin(pType)) ;

	printDebug(FULL_DEBUG, "Out Grid_GPU2::getMin") ;
}

//-------------------------------------------------------------------------------------------------------
Myfloat Grid_GPU2::getMax(Point_type pType)
{
	printDebug(FULL_DEBUG, "In Grid_GPU2::getMax") ;

	// TO DO
	return(Grid::getMax(pType)) ;

	printDebug(FULL_DEBUG, "Out Grid_GPU2::getMax") ;
}

//-------------------------------------------------------------------------------------------------------
Myfloat Grid_GPU2::L1Err(Point_type pointType, const Grid& gridIn) const
{
	printDebug(FULL_DEBUG, "In Grid_GPU2::L1Err") ;

	// TO DO
	return(Grid::L1Err(pointType, gridIn)) ;

	printDebug(FULL_DEBUG, "Out Grid_GPU2::L1Err") ;
}
//-------------------------------------------------------------------------------------------------------
Rtn_code Grid_GPU2::updatePressure(Point_type pType, const Grid& prcGrid,
		const Grid& coefGrid, const Grid& laplaGrid)
{
	printDebug(FULL_DEBUG, "In Grid_GPU2::updatePressure") ;

	// TO DO
	return(Grid::updatePressure(pType, prcGrid, coefGrid, laplaGrid)) ;

	printDebug(FULL_DEBUG, "Out Grid_GPU2::updatePressure") ;
}

} // namespace hpcscan
