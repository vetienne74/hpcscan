
//-------------------------------------------------------------------------------------------------------
// This grid is activated with command line option -testMode CacheBlk
// Derived class from Grid
// Optimized with cache blocking techniques (target CPU)
//-------------------------------------------------------------------------------------------------------

#include "grid_CacheBlk.h"

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

Grid_CacheBlk::Grid_CacheBlk(Grid_type gridTypeIn) : Grid(gridTypeIn)
{
	printDebug(MID_DEBUG, "IN Grid_CacheBlk::Grid_CacheBlk");

	gridMode = "CacheBlk" ;

	printDebug(MID_DEBUG, "OUT Grid_CacheBlk::Grid_CacheBlk");
}

//-------------------------------------------------------------------------------------------------------

Grid_CacheBlk::Grid_CacheBlk(Grid_type gridTypeIn, Dim_type dimIn,
		Myint64 n1InnerIn, Myint64 n2InnerIn, Myint64 n3InnerIn) : Grid(gridTypeIn, dimIn,
				n1InnerIn, n2InnerIn, n3InnerIn)
{
	printDebug(MID_DEBUG, "IN Grid_CacheBlk::Grid_CacheBlk");

	gridMode = "CacheBlk" ;

	printDebug(MID_DEBUG, "OUT Grid_CacheBlk::Grid_CacheBlk");
}

//-------------------------------------------------------------------------------------------------------

Grid_CacheBlk::~Grid_CacheBlk(void)
{
	printDebug(MID_DEBUG, "IN Grid_CacheBlk::~Grid_CacheBlk");

	//delete[] grid_3d ;

	printDebug(MID_DEBUG, "OUT Grid_CacheBlk::~Grid_CacheBlk");
}

//-------------------------------------------------------------------------------------------------------

void Grid_CacheBlk::info(void)
{
	printDebug(FULL_DEBUG, "IN Grid_CacheBlk::info");

	// parent class info
	Grid::info() ;

	// additional info
	const Myint64 cb1 = Config::Instance()->cb1 ;
	const Myint64 cb2 = Config::Instance()->cb2 ;
	const Myint64 cb3 = Config::Instance()->cb3 ;
	printInfo(MASTER, " Cache block size axis1", cb1) ;
	printInfo(MASTER, " Cache block size axis2", cb2) ;
	printInfo(MASTER, " Cache block size axis3", cb3) ;

	printDebug(FULL_DEBUG, "IN Grid_CacheBlk::info");
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_CacheBlk::FD_D2_N1(Point_type pType, const Grid& Wgrid, Myint fdOrder)
{
	printDebug(MID_DEBUG, "IN Grid_CacheBlk::FD_D2_N1");

	// check grids are same size
	if (this->sameSize(Wgrid) != true)
	{
		printError("Grid_CacheBlk::FD_D2_N1, grids have not same size") ;
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

	Myfloat * const w = Wgrid.grid_3d ;
	Myfloat * const u = this->grid_3d ;

	const Myint64 cb1 = Config::Instance()->cb1 ;
	const Myint64 cb2 = Config::Instance()->cb2 ;
	const Myint64 cb3 = Config::Instance()->cb3 ;

	// compute FD along N1
	if (fdOrder == 2)
	{
#pragma omp parallel for collapse(3)
		for (Myint64 b3 = i3Start; b3 <= i3End; b3 += cb3)
		{
			for (Myint64 b2 = i2Start; b2 <= i2End; b2 += cb2)
			{
				for (Myint64 b1 = i1Start; b1 <= i1End; b1 += cb1)
				{
					const Myint64 end3 = min(i3End, b3+cb3-1) ;
					const Myint64 end2 = min(i2End, b2+cb2-1) ;
					const Myint64 end1 = min(i1End, b1+cb1-1) ;

					// loop on grid points inside block
					for (Myint64 i3 = b3; i3 <= end3; i3++)
					{
						for (Myint64 i2 = b2; i2 <= end2; i2++)
						{
#pragma omp simd
							for (Myint64 i1 = b1; i1 <= end1; i1++)
							{

								w[i1+i2*n1+i3*n1*n2] =
										FD_D2_O2_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
							}

						}
					}
				}
			}
		}
	}
	else if (fdOrder == 4)
	{
#pragma omp parallel for collapse(3)
		for (Myint64 b3 = i3Start; b3 <= i3End; b3 += cb3)
		{
			for (Myint64 b2 = i2Start; b2 <= i2End; b2 += cb2)
			{
				for (Myint64 b1 = i1Start; b1 <= i1End; b1 += cb1)
				{
					const Myint64 end3 = min(i3End, b3+cb3-1) ;
					const Myint64 end2 = min(i2End, b2+cb2-1) ;
					const Myint64 end1 = min(i1End, b1+cb1-1) ;

					// loop on grid points inside block
					for (Myint64 i3 = b3; i3 <= end3; i3++)
					{
						for (Myint64 i2 = b2; i2 <= end2; i2++)
						{
#pragma omp simd
							for (Myint64 i1 = b1; i1 <= end1; i1++)
							{

								w[i1+i2*n1+i3*n1*n2] =
										FD_D2_O4_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
							}

						}
					}
				}
			}
		}
	}
	else if (fdOrder == 8)
	{
#pragma omp parallel for collapse(3)
		for (Myint64 b3 = i3Start; b3 <= i3End; b3 += cb3)
		{
			for (Myint64 b2 = i2Start; b2 <= i2End; b2 += cb2)
			{
				for (Myint64 b1 = i1Start; b1 <= i1End; b1 += cb1)
				{
					const Myint64 end3 = min(i3End, b3+cb3-1) ;
					const Myint64 end2 = min(i2End, b2+cb2-1) ;
					const Myint64 end1 = min(i1End, b1+cb1-1) ;

					// loop on grid points inside block
					for (Myint64 i3 = b3; i3 <= end3; i3++)
					{
						for (Myint64 i2 = b2; i2 <= end2; i2++)
						{
#pragma omp simd
							for (Myint64 i1 = b1; i1 <= end1; i1++)
							{

								w[i1+i2*n1+i3*n1*n2] =
										FD_D2_O8_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
							}

						}
					}
				}
			}
		}
	}
	else if (fdOrder == 12)
	{
#pragma omp parallel for collapse(3)
		for (Myint64 b3 = i3Start; b3 <= i3End; b3 += cb3)
		{
			for (Myint64 b2 = i2Start; b2 <= i2End; b2 += cb2)
			{
				for (Myint64 b1 = i1Start; b1 <= i1End; b1 += cb1)
				{
					const Myint64 end3 = min(i3End, b3+cb3-1) ;
					const Myint64 end2 = min(i2End, b2+cb2-1) ;
					const Myint64 end1 = min(i1End, b1+cb1-1) ;

					// loop on grid points inside block
					for (Myint64 i3 = b3; i3 <= end3; i3++)
					{
						for (Myint64 i2 = b2; i2 <= end2; i2++)
						{
#pragma omp simd
							for (Myint64 i1 = b1; i1 <= end1; i1++)
							{

								w[i1+i2*n1+i3*n1*n2] =
										FD_D2_O12_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
							}

						}
					}
				}
			}
		}
	}
	else if (fdOrder == 16)
	{
#pragma omp parallel for collapse(3)
		for (Myint64 b3 = i3Start; b3 <= i3End; b3 += cb3)
		{
			for (Myint64 b2 = i2Start; b2 <= i2End; b2 += cb2)
			{
				for (Myint64 b1 = i1Start; b1 <= i1End; b1 += cb1)
				{
					const Myint64 end3 = min(i3End, b3+cb3-1) ;
					const Myint64 end2 = min(i2End, b2+cb2-1) ;
					const Myint64 end1 = min(i1End, b1+cb1-1) ;

					// loop on grid points inside block
					for (Myint64 i3 = b3; i3 <= end3; i3++)
					{
						for (Myint64 i2 = b2; i2 <= end2; i2++)
						{
#pragma omp simd
							for (Myint64 i1 = b1; i1 <= end1; i1++)
							{

								w[i1+i2*n1+i3*n1*n2] =
										FD_D2_O16_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
							}

						}
					}
				}
			}
		}
	}


	printDebug(MID_DEBUG, "OUT Grid_CacheBlk::FD_D2_N1");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_CacheBlk::FD_D2_N2(Point_type pType, const Grid& Wgrid, Myint fdOrder)
{
	printDebug(MID_DEBUG, "IN Grid_CacheBlk::FD_D2_N2");

	// check grids are same size
	if (this->sameSize(Wgrid) != true)
	{
		printError("Grid_CacheBlk::FD_D2_N2, grids have not same size") ;
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

	Myfloat * const w = Wgrid.grid_3d ;
	Myfloat * const u = this->grid_3d ;

	const Myint64 cb1 = Config::Instance()->cb1 ;
	const Myint64 cb2 = Config::Instance()->cb2 ;
	const Myint64 cb3 = Config::Instance()->cb3 ;

	// compute FD along N2
	if (fdOrder == 2)
	{
#pragma omp parallel for collapse(3)
		for (Myint64 b3 = i3Start; b3 <= i3End; b3 += cb3)
		{
			for (Myint64 b2 = i2Start; b2 <= i2End; b2 += cb2)
			{
				for (Myint64 b1 = i1Start; b1 <= i1End; b1 += cb1)
				{
					const Myint64 end3 = min(i3End, b3+cb3-1) ;
					const Myint64 end2 = min(i2End, b2+cb2-1) ;
					const Myint64 end1 = min(i1End, b1+cb1-1) ;

					// loop on grid points inside block
					for (Myint64 i3 = b3; i3 <= end3; i3++)
					{
						for (Myint64 i2 = b2; i2 <= end2; i2++)
						{
#pragma omp simd
							for (Myint64 i1 = b1; i1 <= end1; i1++)
							{

								w[i1+i2*n1+i3*n1*n2] =
										FD_D2_O2_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
							}

						}
					}
				}
			}
		}
	}
	else if (fdOrder == 4)
	{
#pragma omp parallel for collapse(3)
		for (Myint64 b3 = i3Start; b3 <= i3End; b3 += cb3)
		{
			for (Myint64 b2 = i2Start; b2 <= i2End; b2 += cb2)
			{
				for (Myint64 b1 = i1Start; b1 <= i1End; b1 += cb1)
				{
					const Myint64 end3 = min(i3End, b3+cb3-1) ;
					const Myint64 end2 = min(i2End, b2+cb2-1) ;
					const Myint64 end1 = min(i1End, b1+cb1-1) ;

					// loop on grid points inside block
					for (Myint64 i3 = b3; i3 <= end3; i3++)
					{
						for (Myint64 i2 = b2; i2 <= end2; i2++)
						{
#pragma omp simd
							for (Myint64 i1 = b1; i1 <= end1; i1++)
							{

								w[i1+i2*n1+i3*n1*n2] =
										FD_D2_O4_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
							}

						}
					}
				}
			}
		}
	}
	else if (fdOrder == 8)
	{
#pragma omp parallel for collapse(3)
		for (Myint64 b3 = i3Start; b3 <= i3End; b3 += cb3)
		{
			for (Myint64 b2 = i2Start; b2 <= i2End; b2 += cb2)
			{
				for (Myint64 b1 = i1Start; b1 <= i1End; b1 += cb1)
				{
					const Myint64 end3 = min(i3End, b3+cb3-1) ;
					const Myint64 end2 = min(i2End, b2+cb2-1) ;
					const Myint64 end1 = min(i1End, b1+cb1-1) ;

					// loop on grid points inside block
					for (Myint64 i3 = b3; i3 <= end3; i3++)
					{
						for (Myint64 i2 = b2; i2 <= end2; i2++)
						{
#pragma omp simd
							for (Myint64 i1 = b1; i1 <= end1; i1++)
							{

								w[i1+i2*n1+i3*n1*n2] =
										FD_D2_O8_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
							}

						}
					}
				}
			}
		}
	}
	else if (fdOrder == 12)
	{
#pragma omp parallel for collapse(3)
		for (Myint64 b3 = i3Start; b3 <= i3End; b3 += cb3)
		{
			for (Myint64 b2 = i2Start; b2 <= i2End; b2 += cb2)
			{
				for (Myint64 b1 = i1Start; b1 <= i1End; b1 += cb1)
				{
					const Myint64 end3 = min(i3End, b3+cb3-1) ;
					const Myint64 end2 = min(i2End, b2+cb2-1) ;
					const Myint64 end1 = min(i1End, b1+cb1-1) ;

					// loop on grid points inside block
					for (Myint64 i3 = b3; i3 <= end3; i3++)
					{
						for (Myint64 i2 = b2; i2 <= end2; i2++)
						{
#pragma omp simd
							for (Myint64 i1 = b1; i1 <= end1; i1++)
							{

								w[i1+i2*n1+i3*n1*n2] =
										FD_D2_O12_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
							}

						}
					}
				}
			}
		}
	}
	else if (fdOrder == 16)
	{
#pragma omp parallel for collapse(3)
		for (Myint64 b3 = i3Start; b3 <= i3End; b3 += cb3)
		{
			for (Myint64 b2 = i2Start; b2 <= i2End; b2 += cb2)
			{
				for (Myint64 b1 = i1Start; b1 <= i1End; b1 += cb1)
				{
					const Myint64 end3 = min(i3End, b3+cb3-1) ;
					const Myint64 end2 = min(i2End, b2+cb2-1) ;
					const Myint64 end1 = min(i1End, b1+cb1-1) ;

					// loop on grid points inside block
					for (Myint64 i3 = b3; i3 <= end3; i3++)
					{
						for (Myint64 i2 = b2; i2 <= end2; i2++)
						{
#pragma omp simd
							for (Myint64 i1 = b1; i1 <= end1; i1++)
							{

								w[i1+i2*n1+i3*n1*n2] =
										FD_D2_O16_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
							}

						}
					}
				}
			}
		}
	}

	printDebug(MID_DEBUG, "OUT Grid_CacheBlk::FD_D2_N2");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_CacheBlk::FD_D2_N3(Point_type pType, const Grid& Wgrid, Myint fdOrder)
{
	printDebug(MID_DEBUG, "IN Grid_CacheBlk::FD_D2_N3");

	// check grids are same size
	if (this->sameSize(Wgrid) != true)
	{
		printError("Grid_CacheBlk::FD_D2_N3, grids have not same size") ;
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

	Myfloat * const w = Wgrid.grid_3d ;
	Myfloat * const u = this->grid_3d ;

	const Myint64 cb1 = Config::Instance()->cb1 ;
	const Myint64 cb2 = Config::Instance()->cb2 ;
	const Myint64 cb3 = Config::Instance()->cb3 ;

	// compute FD along N3
	if (fdOrder == 2)
	{
#pragma omp parallel for collapse(3)
		for (Myint64 b3 = i3Start; b3 <= i3End; b3 += cb3)
		{
			for (Myint64 b2 = i2Start; b2 <= i2End; b2 += cb2)
			{
				for (Myint64 b1 = i1Start; b1 <= i1End; b1 += cb1)
				{
					const Myint64 end3 = min(i3End, b3+cb3-1) ;
					const Myint64 end2 = min(i2End, b2+cb2-1) ;
					const Myint64 end1 = min(i1End, b1+cb1-1) ;

					// loop on grid points inside block
					for (Myint64 i3 = b3; i3 <= end3; i3++)
					{
						for (Myint64 i2 = b2; i2 <= end2; i2++)
						{
#pragma omp simd
							for (Myint64 i1 = b1; i1 <= end1; i1++)
							{

								w[i1+i2*n1+i3*n1*n2] =
										FD_D2_O2_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
							}

						}
					}
				}
			}
		}
	}
	else if (fdOrder == 4)
	{
#pragma omp parallel for collapse(3)
		for (Myint64 b3 = i3Start; b3 <= i3End; b3 += cb3)
		{
			for (Myint64 b2 = i2Start; b2 <= i2End; b2 += cb2)
			{
				for (Myint64 b1 = i1Start; b1 <= i1End; b1 += cb1)
				{
					const Myint64 end3 = min(i3End, b3+cb3-1) ;
					const Myint64 end2 = min(i2End, b2+cb2-1) ;
					const Myint64 end1 = min(i1End, b1+cb1-1) ;

					// loop on grid points inside block
					for (Myint64 i3 = b3; i3 <= end3; i3++)
					{
						for (Myint64 i2 = b2; i2 <= end2; i2++)
						{
#pragma omp simd
							for (Myint64 i1 = b1; i1 <= end1; i1++)
							{

								w[i1+i2*n1+i3*n1*n2] =
										FD_D2_O4_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
							}

						}
					}
				}
			}
		}
	}
	else if (fdOrder == 8)
	{
#pragma omp parallel for collapse(3)
		for (Myint64 b3 = i3Start; b3 <= i3End; b3 += cb3)
		{
			for (Myint64 b2 = i2Start; b2 <= i2End; b2 += cb2)
			{
				for (Myint64 b1 = i1Start; b1 <= i1End; b1 += cb1)
				{
					const Myint64 end3 = min(i3End, b3+cb3-1) ;
					const Myint64 end2 = min(i2End, b2+cb2-1) ;
					const Myint64 end1 = min(i1End, b1+cb1-1) ;

					// loop on grid points inside block
					for (Myint64 i3 = b3; i3 <= end3; i3++)
					{
						for (Myint64 i2 = b2; i2 <= end2; i2++)
						{
#pragma omp simd
							for (Myint64 i1 = b1; i1 <= end1; i1++)
							{

								w[i1+i2*n1+i3*n1*n2] =
										FD_D2_O8_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
							}

						}
					}
				}
			}
		}
	}
	else if (fdOrder == 12)
	{
#pragma omp parallel for collapse(3)
		for (Myint64 b3 = i3Start; b3 <= i3End; b3 += cb3)
		{
			for (Myint64 b2 = i2Start; b2 <= i2End; b2 += cb2)
			{
				for (Myint64 b1 = i1Start; b1 <= i1End; b1 += cb1)
				{
					const Myint64 end3 = min(i3End, b3+cb3-1) ;
					const Myint64 end2 = min(i2End, b2+cb2-1) ;
					const Myint64 end1 = min(i1End, b1+cb1-1) ;

					// loop on grid points inside block
					for (Myint64 i3 = b3; i3 <= end3; i3++)
					{
						for (Myint64 i2 = b2; i2 <= end2; i2++)
						{
#pragma omp simd
							for (Myint64 i1 = b1; i1 <= end1; i1++)
							{

								w[i1+i2*n1+i3*n1*n2] =
										FD_D2_O12_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
							}

						}
					}
				}
			}
		}
	}
	else if (fdOrder == 16)
	{
#pragma omp parallel for collapse(3)
		for (Myint64 b3 = i3Start; b3 <= i3End; b3 += cb3)
		{
			for (Myint64 b2 = i2Start; b2 <= i2End; b2 += cb2)
			{
				for (Myint64 b1 = i1Start; b1 <= i1End; b1 += cb1)
				{
					const Myint64 end3 = min(i3End, b3+cb3-1) ;
					const Myint64 end2 = min(i2End, b2+cb2-1) ;
					const Myint64 end1 = min(i1End, b1+cb1-1) ;

					// loop on grid points inside block
					for (Myint64 i3 = b3; i3 <= end3; i3++)
					{
						for (Myint64 i2 = b2; i2 <= end2; i2++)
						{
#pragma omp simd
							for (Myint64 i1 = b1; i1 <= end1; i1++)
							{

								w[i1+i2*n1+i3*n1*n2] =
										FD_D2_O16_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
							}

						}
					}
				}
			}
		}
	}

	printDebug(MID_DEBUG, "OUT Grid_CacheBlk::FD_D2_N3");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_CacheBlk::FD_LAPLACIAN(Point_type pType, const Grid& Wgrid, Myint fdOrder)
{
	printDebug(MID_DEBUG, "IN Grid_CacheBlk::FD_LAPLACIAN");

	// check grids are same size
	if (this->sameSize(Wgrid) != true)
	{
		printError("Grid_CacheBlk::FD_LAPLACIAN, grids have not same size") ;
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

	Myfloat * const w = Wgrid.grid_3d ;
	Myfloat * const u = this->grid_3d ;

	const Myint64 cb1 = Config::Instance()->cb1 ;
	const Myint64 cb2 = Config::Instance()->cb2 ;
	const Myint64 cb3 = Config::Instance()->cb3 ;

	// compute FD Laplacian for 1D
	if (dim == DIM1)
	{
		// same as FD_D2_N1
		this->FD_D2_N1(pType, Wgrid, fdOrder) ;
	}

	// compute FD Laplacian for 2D
	else if (dim == DIM2)
	{
		if (fdOrder == 2)
		{
#pragma omp parallel for collapse(3)
			for (Myint64 b3 = i3Start; b3 <= i3End; b3 += cb3)
			{
				for (Myint64 b2 = i2Start; b2 <= i2End; b2 += cb2)
				{
					for (Myint64 b1 = i1Start; b1 <= i1End; b1 += cb1)
					{
						const Myint64 end3 = min(i3End, b3+cb3-1) ;
						const Myint64 end2 = min(i2End, b2+cb2-1) ;
						const Myint64 end1 = min(i1End, b1+cb1-1) ;

						// loop on grid points inside block
						for (Myint64 i3 = b3; i3 <= end3; i3++)
						{
							for (Myint64 i2 = b2; i2 <= end2; i2++)
							{
#pragma omp simd
								for (Myint64 i1 = b1; i1 <= end1; i1++)
								{

									w[i1+i2*n1+i3*n1*n2] =
											FD_D2_O2_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
											+ FD_D2_O2_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
								}

							}
						}
					}
				}
			}
		}
		else if (fdOrder == 4)
		{
#pragma omp parallel for collapse(3)
			for (Myint64 b3 = i3Start; b3 <= i3End; b3 += cb3)
			{
				for (Myint64 b2 = i2Start; b2 <= i2End; b2 += cb2)
				{
					for (Myint64 b1 = i1Start; b1 <= i1End; b1 += cb1)
					{
						const Myint64 end3 = min(i3End, b3+cb3-1) ;
						const Myint64 end2 = min(i2End, b2+cb2-1) ;
						const Myint64 end1 = min(i1End, b1+cb1-1) ;

						// loop on grid points inside block
						for (Myint64 i3 = b3; i3 <= end3; i3++)
						{
							for (Myint64 i2 = b2; i2 <= end2; i2++)
							{
#pragma omp simd
								for (Myint64 i1 = b1; i1 <= end1; i1++)
								{

									w[i1+i2*n1+i3*n1*n2] =
											FD_D2_O4_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
											+ FD_D2_O4_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
								}

							}
						}
					}
				}
			}
		}
		else if (fdOrder == 8)
		{
#pragma omp parallel for collapse(3)
			for (Myint64 b3 = i3Start; b3 <= i3End; b3 += cb3)
			{
				for (Myint64 b2 = i2Start; b2 <= i2End; b2 += cb2)
				{
					for (Myint64 b1 = i1Start; b1 <= i1End; b1 += cb1)
					{
						const Myint64 end3 = min(i3End, b3+cb3-1) ;
						const Myint64 end2 = min(i2End, b2+cb2-1) ;
						const Myint64 end1 = min(i1End, b1+cb1-1) ;

						// loop on grid points inside block
						for (Myint64 i3 = b3; i3 <= end3; i3++)
						{
							for (Myint64 i2 = b2; i2 <= end2; i2++)
							{
#pragma omp simd
								for (Myint64 i1 = b1; i1 <= end1; i1++)
								{

									w[i1+i2*n1+i3*n1*n2] =
											FD_D2_O8_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
											+ FD_D2_O8_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
								}

							}
						}
					}
				}
			}
		}
		else if (fdOrder == 12)
		{
#pragma omp parallel for collapse(3)
			for (Myint64 b3 = i3Start; b3 <= i3End; b3 += cb3)
			{
				for (Myint64 b2 = i2Start; b2 <= i2End; b2 += cb2)
				{
					for (Myint64 b1 = i1Start; b1 <= i1End; b1 += cb1)
					{
						const Myint64 end3 = min(i3End, b3+cb3-1) ;
						const Myint64 end2 = min(i2End, b2+cb2-1) ;
						const Myint64 end1 = min(i1End, b1+cb1-1) ;

						// loop on grid points inside block
						for (Myint64 i3 = b3; i3 <= end3; i3++)
						{
							for (Myint64 i2 = b2; i2 <= end2; i2++)
							{
#pragma omp simd
								for (Myint64 i1 = b1; i1 <= end1; i1++)
								{

									w[i1+i2*n1+i3*n1*n2] =
											FD_D2_O12_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
											+ FD_D2_O12_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
								}

							}
						}
					}
				}
			}
		}
		else if (fdOrder == 16)
		{
#pragma omp parallel for collapse(3)
			for (Myint64 b3 = i3Start; b3 <= i3End; b3 += cb3)
			{
				for (Myint64 b2 = i2Start; b2 <= i2End; b2 += cb2)
				{
					for (Myint64 b1 = i1Start; b1 <= i1End; b1 += cb1)
					{
						const Myint64 end3 = min(i3End, b3+cb3-1) ;
						const Myint64 end2 = min(i2End, b2+cb2-1) ;
						const Myint64 end1 = min(i1End, b1+cb1-1) ;

						// loop on grid points inside block
						for (Myint64 i3 = b3; i3 <= end3; i3++)
						{
							for (Myint64 i2 = b2; i2 <= end2; i2++)
							{
#pragma omp simd
								for (Myint64 i1 = b1; i1 <= end1; i1++)
								{

									w[i1+i2*n1+i3*n1*n2] =
											FD_D2_O16_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
											+ FD_D2_O16_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
								}

							}
						}
					}
				}
			}
		}
	}

	// compute FD Laplacian for 3D
	else if (dim == DIM3)
	{
		if (fdOrder == 2)
		{
#pragma omp parallel for collapse(3)
			for (Myint64 b3 = i3Start; b3 <= i3End; b3 += cb3)
			{
				for (Myint64 b2 = i2Start; b2 <= i2End; b2 += cb2)
				{
					for (Myint64 b1 = i1Start; b1 <= i1End; b1 += cb1)
					{
						const Myint64 end3 = min(i3End, b3+cb3-1) ;
						const Myint64 end2 = min(i2End, b2+cb2-1) ;
						const Myint64 end1 = min(i1End, b1+cb1-1) ;

						// loop on grid points inside block
						for (Myint64 i3 = b3; i3 <= end3; i3++)
						{
							for (Myint64 i2 = b2; i2 <= end2; i2++)
							{
#pragma omp simd
								for (Myint64 i1 = b1; i1 <= end1; i1++)
								{

									w[i1+i2*n1+i3*n1*n2] =
											FD_D2_O2_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
											+ FD_D2_O2_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
											+ FD_D2_O2_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
								}

							}
						}
					}
				}
			}
		}
		else if (fdOrder == 4)
		{
#pragma omp parallel for collapse(3)
			for (Myint64 b3 = i3Start; b3 <= i3End; b3 += cb3)
			{
				for (Myint64 b2 = i2Start; b2 <= i2End; b2 += cb2)
				{
					for (Myint64 b1 = i1Start; b1 <= i1End; b1 += cb1)
					{
						const Myint64 end3 = min(i3End, b3+cb3-1) ;
						const Myint64 end2 = min(i2End, b2+cb2-1) ;
						const Myint64 end1 = min(i1End, b1+cb1-1) ;

						// loop on grid points inside block
						for (Myint64 i3 = b3; i3 <= end3; i3++)
						{
							for (Myint64 i2 = b2; i2 <= end2; i2++)
							{
#pragma omp simd
								for (Myint64 i1 = b1; i1 <= end1; i1++)
								{

									w[i1+i2*n1+i3*n1*n2] =
											FD_D2_O4_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
											+ FD_D2_O4_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
											+ FD_D2_O4_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
								}

							}
						}
					}
				}
			}
		}
		else if (fdOrder == 8)
		{
#pragma omp parallel for collapse(3)
			for (Myint64 b3 = i3Start; b3 <= i3End; b3 += cb3)
			{
				for (Myint64 b2 = i2Start; b2 <= i2End; b2 += cb2)
				{
					for (Myint64 b1 = i1Start; b1 <= i1End; b1 += cb1)
					{
						const Myint64 end3 = min(i3End, b3+cb3-1) ;
						const Myint64 end2 = min(i2End, b2+cb2-1) ;
						const Myint64 end1 = min(i1End, b1+cb1-1) ;

						// loop on grid points inside block
						for (Myint64 i3 = b3; i3 <= end3; i3++)
						{
							for (Myint64 i2 = b2; i2 <= end2; i2++)
							{
#pragma omp simd
								for (Myint64 i1 = b1; i1 <= end1; i1++)
								{

									w[i1+i2*n1+i3*n1*n2] =
											FD_D2_O8_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
											+ FD_D2_O8_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
											+ FD_D2_O8_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
								}

							}
						}
					}
				}
			}
		}
		else if (fdOrder == 12)
		{
#pragma omp parallel for collapse(3)
			for (Myint64 b3 = i3Start; b3 <= i3End; b3 += cb3)
			{
				for (Myint64 b2 = i2Start; b2 <= i2End; b2 += cb2)
				{
					for (Myint64 b1 = i1Start; b1 <= i1End; b1 += cb1)
					{
						const Myint64 end3 = min(i3End, b3+cb3-1) ;
						const Myint64 end2 = min(i2End, b2+cb2-1) ;
						const Myint64 end1 = min(i1End, b1+cb1-1) ;

						// loop on grid points inside block
						for (Myint64 i3 = b3; i3 <= end3; i3++)
						{
							for (Myint64 i2 = b2; i2 <= end2; i2++)
							{
#pragma omp simd
								for (Myint64 i1 = b1; i1 <= end1; i1++)
								{

									w[i1+i2*n1+i3*n1*n2] =
											FD_D2_O12_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
											+ FD_D2_O12_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
											+ FD_D2_O12_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
								}

							}
						}
					}
				}
			}
		}
		else if (fdOrder == 16)
		{
#pragma omp parallel for collapse(3)
			for (Myint64 b3 = i3Start; b3 <= i3End; b3 += cb3)
			{
				for (Myint64 b2 = i2Start; b2 <= i2End; b2 += cb2)
				{
					for (Myint64 b1 = i1Start; b1 <= i1End; b1 += cb1)
					{
						const Myint64 end3 = min(i3End, b3+cb3-1) ;
						const Myint64 end2 = min(i2End, b2+cb2-1) ;
						const Myint64 end1 = min(i1End, b1+cb1-1) ;

						// loop on grid points inside block
						for (Myint64 i3 = b3; i3 <= end3; i3++)
						{
							for (Myint64 i2 = b2; i2 <= end2; i2++)
							{
#pragma omp simd
								for (Myint64 i1 = b1; i1 <= end1; i1++)
								{

									w[i1+i2*n1+i3*n1*n2] =
											FD_D2_O16_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
											+ FD_D2_O16_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
											+ FD_D2_O16_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
								}

							}
						}
					}
				}
			}
		}
	}

	printDebug(MID_DEBUG, "OUT Grid_CacheBlk::FD_LAPLACIAN");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_CacheBlk::computePressureWithFD(Grid& prcGridIn, Grid& coefGridIn, Myint fdOrder)
{

	printDebug(FULL_DEBUG, "In Grid_CacheBlk::computePressureWithFD") ;

	// check grids are same size
	if (this->sameSize(prcGridIn) != true)
	{
		printError("In Grid_CacheBlk::computePressureWithFD, grids have not same size") ;
		return(RTN_CODE_KO) ;
	}

	Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
	getGridIndex(INNER_POINTS, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End) ;

	Myfloat * const prn  = this->grid_3d ;
	Myfloat * const prc  = prcGridIn.grid_3d ;
	Myfloat * const coef = coefGridIn.grid_3d ;

	const Myfloat inv_d1  = Myfloat(1.0) / d1 ;
	const Myfloat inv_d2  = Myfloat(1.0) / d2 ;
	const Myfloat inv_d3  = Myfloat(1.0) / d3 ;

	const Myfloat inv2_d1 = inv_d1 * inv_d1 ;
	const Myfloat inv2_d2 = inv_d2 * inv_d2 ;
	const Myfloat inv2_d3 = inv_d3 * inv_d3 ;

	const Myint64 cb1 = Config::Instance()->cb1 ;
	const Myint64 cb2 = Config::Instance()->cb2 ;
	const Myint64 cb3 = Config::Instance()->cb3 ;

	// compute FD for 1D
	if (dim == DIM1)
	{
		if (fdOrder == 2)
		{
#pragma omp parallel for collapse(3)
			for (Myint64 b3 = i3Start; b3 <= i3End; b3 += cb3)
			{
				for (Myint64 b2 = i2Start; b2 <= i2End; b2 += cb2)
				{
					for (Myint64 b1 = i1Start; b1 <= i1End; b1 += cb1)
					{
						const Myint64 end3 = min(i3End, b3+cb3-1) ;
						const Myint64 end2 = min(i2End, b2+cb2-1) ;
						const Myint64 end1 = min(i1End, b1+cb1-1) ;

						// loop on grid points inside block
						for (Myint64 i3 = b3; i3 <= end3; i3++)
						{
							for (Myint64 i2 = b2; i2 <= end2; i2++)
							{
#pragma omp simd
								for (Myint64 i1 = b1; i1 <= end1; i1++)
								{
									prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
											coef[i1+i2*n1+i3*n1*n2] *
											FD_D2_O2_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
								}
							}
						}
					}
				}
			}
		}
		else if (fdOrder == 4)
		{
#pragma omp parallel for collapse(3)
			for (Myint64 b3 = i3Start; b3 <= i3End; b3 += cb3)
			{
				for (Myint64 b2 = i2Start; b2 <= i2End; b2 += cb2)
				{
					for (Myint64 b1 = i1Start; b1 <= i1End; b1 += cb1)
					{
						const Myint64 end3 = min(i3End, b3+cb3-1) ;
						const Myint64 end2 = min(i2End, b2+cb2-1) ;
						const Myint64 end1 = min(i1End, b1+cb1-1) ;

						// loop on grid points inside block
						for (Myint64 i3 = b3; i3 <= end3; i3++)
						{
							for (Myint64 i2 = b2; i2 <= end2; i2++)
							{
#pragma omp simd
								for (Myint64 i1 = b1; i1 <= end1; i1++)
								{
									prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
											coef[i1+i2*n1+i3*n1*n2] *
											FD_D2_O4_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
								}
							}
						}
					}
				}
			}
		}
		else if (fdOrder == 8)
		{
#pragma omp parallel for collapse(3)
			for (Myint64 b3 = i3Start; b3 <= i3End; b3 += cb3)
			{
				for (Myint64 b2 = i2Start; b2 <= i2End; b2 += cb2)
				{
					for (Myint64 b1 = i1Start; b1 <= i1End; b1 += cb1)
					{
						const Myint64 end3 = min(i3End, b3+cb3-1) ;
						const Myint64 end2 = min(i2End, b2+cb2-1) ;
						const Myint64 end1 = min(i1End, b1+cb1-1) ;

						// loop on grid points inside block
						for (Myint64 i3 = b3; i3 <= end3; i3++)
						{
							for (Myint64 i2 = b2; i2 <= end2; i2++)
							{
#pragma omp simd
								for (Myint64 i1 = b1; i1 <= end1; i1++)
								{
									prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
											coef[i1+i2*n1+i3*n1*n2] *
											FD_D2_O8_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
								}
							}
						}
					}
				}
			}
		}
		else if (fdOrder == 12)
		{
#pragma omp parallel for collapse(3)
			for (Myint64 b3 = i3Start; b3 <= i3End; b3 += cb3)
			{
				for (Myint64 b2 = i2Start; b2 <= i2End; b2 += cb2)
				{
					for (Myint64 b1 = i1Start; b1 <= i1End; b1 += cb1)
					{
						const Myint64 end3 = min(i3End, b3+cb3-1) ;
						const Myint64 end2 = min(i2End, b2+cb2-1) ;
						const Myint64 end1 = min(i1End, b1+cb1-1) ;

						// loop on grid points inside block
						for (Myint64 i3 = b3; i3 <= end3; i3++)
						{
							for (Myint64 i2 = b2; i2 <= end2; i2++)
							{
#pragma omp simd
								for (Myint64 i1 = b1; i1 <= end1; i1++)
								{
									prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
											coef[i1+i2*n1+i3*n1*n2] *
											FD_D2_O12_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
								}
							}
						}
					}
				}
			}
		}
		else if (fdOrder == 16)
		{
#pragma omp parallel for collapse(3)
			for (Myint64 b3 = i3Start; b3 <= i3End; b3 += cb3)
			{
				for (Myint64 b2 = i2Start; b2 <= i2End; b2 += cb2)
				{
					for (Myint64 b1 = i1Start; b1 <= i1End; b1 += cb1)
					{
						const Myint64 end3 = min(i3End, b3+cb3-1) ;
						const Myint64 end2 = min(i2End, b2+cb2-1) ;
						const Myint64 end1 = min(i1End, b1+cb1-1) ;

						// loop on grid points inside block
						for (Myint64 i3 = b3; i3 <= end3; i3++)
						{
							for (Myint64 i2 = b2; i2 <= end2; i2++)
							{
#pragma omp simd
								for (Myint64 i1 = b1; i1 <= end1; i1++)
								{
									prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
											coef[i1+i2*n1+i3*n1*n2] *
											FD_D2_O16_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
								}
							}
						}
					}
				}
			}
		}

	}

	// compute FD for 2D
	else if (dim == DIM2)
	{
		if (fdOrder == 2)
		{
#pragma omp parallel for collapse(3)
			for (Myint64 b3 = i3Start; b3 <= i3End; b3 += cb3)
			{
				for (Myint64 b2 = i2Start; b2 <= i2End; b2 += cb2)
				{
					for (Myint64 b1 = i1Start; b1 <= i1End; b1 += cb1)
					{
						const Myint64 end3 = min(i3End, b3+cb3-1) ;
						const Myint64 end2 = min(i2End, b2+cb2-1) ;
						const Myint64 end1 = min(i1End, b1+cb1-1) ;

						// loop on grid points inside block
						for (Myint64 i3 = b3; i3 <= end3; i3++)
						{
							for (Myint64 i2 = b2; i2 <= end2; i2++)
							{
#pragma omp simd
								for (Myint64 i1 = b1; i1 <= end1; i1++)
								{
									prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
											coef[i1+i2*n1+i3*n1*n2] *
											(FD_D2_O2_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
													+ FD_D2_O2_N2(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)) ;
								}
							}
						}
					}
				}
			}
		}
		else if (fdOrder == 4)
		{
#pragma omp parallel for collapse(3)
			for (Myint64 b3 = i3Start; b3 <= i3End; b3 += cb3)
			{
				for (Myint64 b2 = i2Start; b2 <= i2End; b2 += cb2)
				{
					for (Myint64 b1 = i1Start; b1 <= i1End; b1 += cb1)
					{
						const Myint64 end3 = min(i3End, b3+cb3-1) ;
						const Myint64 end2 = min(i2End, b2+cb2-1) ;
						const Myint64 end1 = min(i1End, b1+cb1-1) ;

						// loop on grid points inside block
						for (Myint64 i3 = b3; i3 <= end3; i3++)
						{
							for (Myint64 i2 = b2; i2 <= end2; i2++)
							{
#pragma omp simd
								for (Myint64 i1 = b1; i1 <= end1; i1++)
								{
									prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
											coef[i1+i2*n1+i3*n1*n2] *
											(FD_D2_O4_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
													+ FD_D2_O4_N2(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)) ;
								}
							}
						}
					}
				}
			}
		}
		else if (fdOrder == 8)
		{
#pragma omp parallel for collapse(3)
			for (Myint64 b3 = i3Start; b3 <= i3End; b3 += cb3)
			{
				for (Myint64 b2 = i2Start; b2 <= i2End; b2 += cb2)
				{
					for (Myint64 b1 = i1Start; b1 <= i1End; b1 += cb1)
					{
						const Myint64 end3 = min(i3End, b3+cb3-1) ;
						const Myint64 end2 = min(i2End, b2+cb2-1) ;
						const Myint64 end1 = min(i1End, b1+cb1-1) ;

						// loop on grid points inside block
						for (Myint64 i3 = b3; i3 <= end3; i3++)
						{
							for (Myint64 i2 = b2; i2 <= end2; i2++)
							{
#pragma omp simd
								for (Myint64 i1 = b1; i1 <= end1; i1++)
								{
									prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
											coef[i1+i2*n1+i3*n1*n2] *
											(FD_D2_O8_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
													+ FD_D2_O8_N2(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)) ;
								}
							}
						}
					}
				}
			}
		}
		else if (fdOrder == 12)
		{
#pragma omp parallel for collapse(3)
			for (Myint64 b3 = i3Start; b3 <= i3End; b3 += cb3)
			{
				for (Myint64 b2 = i2Start; b2 <= i2End; b2 += cb2)
				{
					for (Myint64 b1 = i1Start; b1 <= i1End; b1 += cb1)
					{
						const Myint64 end3 = min(i3End, b3+cb3-1) ;
						const Myint64 end2 = min(i2End, b2+cb2-1) ;
						const Myint64 end1 = min(i1End, b1+cb1-1) ;

						// loop on grid points inside block
						for (Myint64 i3 = b3; i3 <= end3; i3++)
						{
							for (Myint64 i2 = b2; i2 <= end2; i2++)
							{
#pragma omp simd
								for (Myint64 i1 = b1; i1 <= end1; i1++)
								{
									prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
											coef[i1+i2*n1+i3*n1*n2] *
											(FD_D2_O12_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
													+ FD_D2_O12_N2(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)) ;
								}
							}
						}
					}
				}
			}
		}
		else if (fdOrder == 16)
		{
#pragma omp parallel for collapse(3)
			for (Myint64 b3 = i3Start; b3 <= i3End; b3 += cb3)
			{
				for (Myint64 b2 = i2Start; b2 <= i2End; b2 += cb2)
				{
					for (Myint64 b1 = i1Start; b1 <= i1End; b1 += cb1)
					{
						const Myint64 end3 = min(i3End, b3+cb3-1) ;
						const Myint64 end2 = min(i2End, b2+cb2-1) ;
						const Myint64 end1 = min(i1End, b1+cb1-1) ;

						// loop on grid points inside block
						for (Myint64 i3 = b3; i3 <= end3; i3++)
						{
							for (Myint64 i2 = b2; i2 <= end2; i2++)
							{
#pragma omp simd
								for (Myint64 i1 = b1; i1 <= end1; i1++)
								{
									prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
											coef[i1+i2*n1+i3*n1*n2] *
											(FD_D2_O16_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
													+ FD_D2_O16_N2(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)) ;
								}
							}
						}
					}
				}
			}
		}
	}

	// compute FD for 3D
	else if (dim == DIM3)
	{
		if (fdOrder == 2)
		{
#pragma omp parallel for collapse(3)
			for (Myint64 b3 = i3Start; b3 <= i3End; b3 += cb3)
			{
				for (Myint64 b2 = i2Start; b2 <= i2End; b2 += cb2)
				{
					for (Myint64 b1 = i1Start; b1 <= i1End; b1 += cb1)
					{
						const Myint64 end3 = min(i3End, b3+cb3-1) ;
						const Myint64 end2 = min(i2End, b2+cb2-1) ;
						const Myint64 end1 = min(i1End, b1+cb1-1) ;

						// loop on grid points inside block
						for (Myint64 i3 = b3; i3 <= end3; i3++)
						{
							for (Myint64 i2 = b2; i2 <= end2; i2++)
							{
#pragma omp simd
								for (Myint64 i1 = b1; i1 <= end1; i1++)
								{
									prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
											coef[i1+i2*n1+i3*n1*n2] *
											(FD_D2_O2_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
													+ FD_D2_O2_N2(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
													+ FD_D2_O2_N3(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)) ;
								}
							}
						}
					}
				}
			}
		}
		else if (fdOrder == 4)
		{
#pragma omp parallel for collapse(3)
			for (Myint64 b3 = i3Start; b3 <= i3End; b3 += cb3)
			{
				for (Myint64 b2 = i2Start; b2 <= i2End; b2 += cb2)
				{
					for (Myint64 b1 = i1Start; b1 <= i1End; b1 += cb1)
					{
						const Myint64 end3 = min(i3End, b3+cb3-1) ;
						const Myint64 end2 = min(i2End, b2+cb2-1) ;
						const Myint64 end1 = min(i1End, b1+cb1-1) ;

						// loop on grid points inside block
						for (Myint64 i3 = b3; i3 <= end3; i3++)
						{
							for (Myint64 i2 = b2; i2 <= end2; i2++)
							{
#pragma omp simd
								for (Myint64 i1 = b1; i1 <= end1; i1++)
								{
									prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
											coef[i1+i2*n1+i3*n1*n2] *
											(FD_D2_O4_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
													+ FD_D2_O4_N2(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
													+ FD_D2_O4_N3(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)) ;
								}
							}
						}
					}
				}
			}
		}
		else if (fdOrder == 8)
		{
#pragma omp parallel for collapse(3)
			for (Myint64 b3 = i3Start; b3 <= i3End; b3 += cb3)
			{
				for (Myint64 b2 = i2Start; b2 <= i2End; b2 += cb2)
				{
					for (Myint64 b1 = i1Start; b1 <= i1End; b1 += cb1)
					{
						const Myint64 end3 = min(i3End, b3+cb3-1) ;
						const Myint64 end2 = min(i2End, b2+cb2-1) ;
						const Myint64 end1 = min(i1End, b1+cb1-1) ;

						// loop on grid points inside block
						for (Myint64 i3 = b3; i3 <= end3; i3++)
						{
							for (Myint64 i2 = b2; i2 <= end2; i2++)
							{
#pragma omp simd
								for (Myint64 i1 = b1; i1 <= end1; i1++)
								{
									prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
											coef[i1+i2*n1+i3*n1*n2] *
											(FD_D2_O8_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
													+ FD_D2_O8_N2(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
													+ FD_D2_O8_N3(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)) ;
								}
							}
						}
					}
				}
			}
		}
		else if (fdOrder == 12)
		{
#pragma omp parallel for collapse(3)
			for (Myint64 b3 = i3Start; b3 <= i3End; b3 += cb3)
			{
				for (Myint64 b2 = i2Start; b2 <= i2End; b2 += cb2)
				{
					for (Myint64 b1 = i1Start; b1 <= i1End; b1 += cb1)
					{
						const Myint64 end3 = min(i3End, b3+cb3-1) ;
						const Myint64 end2 = min(i2End, b2+cb2-1) ;
						const Myint64 end1 = min(i1End, b1+cb1-1) ;

						// loop on grid points inside block
						for (Myint64 i3 = b3; i3 <= end3; i3++)
						{
							for (Myint64 i2 = b2; i2 <= end2; i2++)
							{
#pragma omp simd
								for (Myint64 i1 = b1; i1 <= end1; i1++)
								{
									prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
											coef[i1+i2*n1+i3*n1*n2] *
											(FD_D2_O12_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
													+ FD_D2_O12_N2(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
													+ FD_D2_O12_N3(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)) ;
								}
							}
						}
					}
				}
			}
		}
		else if (fdOrder == 16)
		{
#pragma omp parallel for collapse(3)
			for (Myint64 b3 = i3Start; b3 <= i3End; b3 += cb3)
			{
				for (Myint64 b2 = i2Start; b2 <= i2End; b2 += cb2)
				{
					for (Myint64 b1 = i1Start; b1 <= i1End; b1 += cb1)
					{
						const Myint64 end3 = min(i3End, b3+cb3-1) ;
						const Myint64 end2 = min(i2End, b2+cb2-1) ;
						const Myint64 end1 = min(i1End, b1+cb1-1) ;

						// loop on grid points inside block
						for (Myint64 i3 = b3; i3 <= end3; i3++)
						{
							for (Myint64 i2 = b2; i2 <= end2; i2++)
							{
#pragma omp simd
								for (Myint64 i1 = b1; i1 <= end1; i1++)
								{
									prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
											coef[i1+i2*n1+i3*n1*n2] *
											(FD_D2_O16_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
													+ FD_D2_O16_N2(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
													+ FD_D2_O16_N3(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)) ;
								}
							}
						}
					}
				}
			}
		}
	}

	printDebug(FULL_DEBUG, "Out Grid_CacheBlk::computePressureWithFD") ;
	return(RTN_CODE_OK) ;
}



} // namespace hpcscan
