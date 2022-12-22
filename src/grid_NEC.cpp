
//-------------------------------------------------------------------------------------------------------
// This grid is activated with command line option -testMode NEC
// Derived class from Grid
// NEC compiler directives (target NEC SX-Aurora TSUBASA Vector Engine)
//-------------------------------------------------------------------------------------------------------

#include "grid_NEC.h"

#include <cassert>
#include <cfloat>  // for FLT_MAX ;
#include <cmath>   // for fabs
#include <cstddef> // for NULL
#include <fcntl.h> // for C open file operations
#include <fstream>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h> // for C read/write file operations


#include "mpi.h"
#include <sca.h>
#ifdef _FTRACE
#include <ftrace.h>
#endif

#include "config.h"
#include "constant.h"
#include "fdm.h"
#include "global.h"
#include "output_report.h"

#ifdef _DOUBLE_PRECISION_
#define sca_utility_optimize_leading sca_utility_optimize_leading_d
#else
#define sca_utility_optimize_leading sca_utility_optimize_leading_s
#endif

using namespace std;

namespace hpcscan {

//-------------------------------------------------------------------------------------------------------

Grid_NEC::Grid_NEC(Grid_type gridTypeIn) : Grid(gridTypeIn)
{
	printDebug(MID_DEBUG, "IN Grid_NEC::Grid_NEC");

	gridMode = GRID_MODE_NEC ;
	flag_packed_stencil = false ;
	tmp_grid_3d = nullptr ;

	printDebug(MID_DEBUG, "OUT Grid_NEC::Grid_NEC");
}

//-------------------------------------------------------------------------------------------------------

Grid_NEC::Grid_NEC(Grid_type gridTypeIn, Dim_type dimIn,
		Myint64 n1InnerIn, Myint64 n2InnerIn, Myint64 n3InnerIn) : Grid(gridTypeIn, dimIn,
				n1InnerIn, n2InnerIn, n3InnerIn)
{
	printDebug(MID_DEBUG, "IN Grid_NEC::Grid_NEC");

	gridMode = GRID_MODE_NEC ;
	flag_packed_stencil = false ;
	tmp_grid_3d = nullptr ;

	printDebug(MID_DEBUG, "OUT Grid_NEC::Grid_NEC");
}

//-------------------------------------------------------------------------------------------------------

Grid_NEC::~Grid_NEC(void)
{
	printDebug(MID_DEBUG, "IN Grid_NEC::~Grid_NEC");

	printDebug(MID_DEBUG, "OUT Grid_NEC::~Grid_NEC");
}

//-------------------------------------------------------------------------------------------------------

void Grid_NEC::info(void)
{
	printDebug(FULL_DEBUG, "IN Grid_NEC::info");

	// parent class info
	Grid::info() ;

	// additional info specific to NEC
	if (flag_packed_stencil)
	{
		printInfo(MASTER," Packed data\t", "ON") ;
	}
	else
	{
		printInfo(MASTER," Packed data\t", "OFF") ;
	}

	printDebug(FULL_DEBUG, "OUT Grid_NEC::info");
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_NEC::FD_D2_N1(Point_type pType, const Grid& Wgrid, Myint fdOrder)
{
	printDebug(MID_DEBUG, "IN Grid_NEC::FD_D2_N1");

	// check grids are same size
	if (this->sameSize(Wgrid) != true)
	{
		printError("Grid_NEC::FD_D2_N1, grids have not same size") ;
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

	// Workaround for compiler not correctly applying packed_stencil directive.
	const int ln1 = n1;
	const int ln2 = n2;
	const int ln3 = n3;

	Myfloat * const w = Wgrid.grid_3d ;
	Myfloat * const u = this->grid_3d ;

	// compute FD along N1
#ifndef _DOUBLE_PRECISION_
	// if (flag_packed_stencil)
	if (flag_packed_stencil)
	{
		if (fdOrder == 2)
		{
#pragma omp parallel for collapse(2)
			for (Myint64 i3 = i3Start; i3<= i3End; i3++)
			{
				for (Myint64 i2 = i2Start; i2<= i2End; i2++)
				{
#pragma _NEC packed_stencil
					for (Myint64 i1 = i1Start; i1<= i1End; i1++)
					{
						w[i1+i2*ln1+i3*ln1*ln2] =
								FD_D2_O2_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, ln1, ln2, ln3) ;
					}
				}
			}
		}
		else if (fdOrder == 4)
		{
#pragma omp parallel for collapse(2)
			for (Myint64 i3 = i3Start; i3<= i3End; i3++)
			{
				for (Myint64 i2 = i2Start; i2<= i2End; i2++)
				{
#pragma _NEC packed_stencil
					for (Myint64 i1 = i1Start; i1<= i1End; i1++)
					{
						w[i1+i2*ln1+i3*ln1*ln2] =
								FD_D2_O4_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, ln1, ln2, ln3) ;
					}
				}
			}
		}
		else if (fdOrder == 6)
		{
#pragma omp parallel for collapse(2)
			for (Myint64 i3 = i3Start; i3<= i3End; i3++)
			{
				for (Myint64 i2 = i2Start; i2<= i2End; i2++)
				{
#pragma _NEC packed_stencil
					for (Myint64 i1 = i1Start; i1<= i1End; i1++)
					{
						w[i1+i2*ln1+i3*ln1*ln2] =
								FD_D2_O6_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, ln1, ln2, ln3) ;
					}
				}
			}
		}
		else if (fdOrder == 8)
		{
#pragma omp parallel for collapse(2)
			for (Myint64 i3 = i3Start; i3<= i3End; i3++)
			{
				for (Myint64 i2 = i2Start; i2<= i2End; i2++)
				{
#pragma _NEC packed_stencil
					for (Myint64 i1 = i1Start; i1<= i1End; i1++)
					{
						w[i1+i2*ln1+i3*ln1*ln2] =
								FD_D2_O8_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, ln1, ln2, ln3) ;
					}
				}
			}
		}
		else if (fdOrder == 10)
		{
#pragma omp parallel for collapse(2)
			for (Myint64 i3 = i3Start; i3<= i3End; i3++)
			{
				for (Myint64 i2 = i2Start; i2<= i2End; i2++)
				{
#pragma _NEC packed_stencil
					for (Myint64 i1 = i1Start; i1<= i1End; i1++)
					{
						w[i1+i2*ln1+i3*ln1*ln2] =
								FD_D2_O10_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, ln1, ln2, ln3) ;
					}
				}
			}
		}
		else if (fdOrder == 12)
		{
#pragma omp parallel for collapse(2)
			for (Myint64 i3 = i3Start; i3<= i3End; i3++)
			{
				for (Myint64 i2 = i2Start; i2<= i2End; i2++)
				{
#pragma _NEC packed_stencil
					for (Myint64 i1 = i1Start; i1<= i1End; i1++)
					{
						w[i1+i2*ln1+i3*ln1*ln2] =
								FD_D2_O12_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, ln1, ln2, ln3) ;
					}
				}
			}
		}
		else if (fdOrder == 14)
		{
#pragma omp parallel for collapse(2)
			for (Myint64 i3 = i3Start; i3<= i3End; i3++)
			{
				for (Myint64 i2 = i2Start; i2<= i2End; i2++)
				{
#pragma _NEC packed_stencil
					for (Myint64 i1 = i1Start; i1<= i1End; i1++)
					{
						w[i1+i2*ln1+i3*ln1*ln2] =
								FD_D2_O14_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, ln1, ln2, ln3) ;
					}
				}
			}
		}
		else if (fdOrder == 16)
		{
#pragma omp parallel for collapse(2)
			for (Myint64 i3 = i3Start; i3<= i3End; i3++)
			{
				for (Myint64 i2 = i2Start; i2<= i2End; i2++)
				{
#pragma _NEC packed_stencil
					for (Myint64 i1 = i1Start; i1<= i1End; i1++)
					{
						w[i1+i2*ln1+i3*ln1*ln2] =
								FD_D2_O16_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, ln1, ln2, ln3) ;
					}
				}
			}
		}
	}
	// else if (flag_packed_stencil)
	else
#endif
	{
		if (fdOrder == 2)
		{
#pragma omp parallel for collapse(2)
			for (Myint64 i3 = i3Start; i3<= i3End; i3++)
			{
				for (Myint64 i2 = i2Start; i2<= i2End; i2++)
				{
					for (Myint64 i1 = i1Start; i1<= i1End; i1++)
					{
						w[i1+i2*ln1+i3*ln1*ln2] =
								FD_D2_O2_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, ln1, ln2, ln3) ;
					}
				}
			}
		}
		else if (fdOrder == 4)
		{
#pragma omp parallel for collapse(2)
			for (Myint64 i3 = i3Start; i3<= i3End; i3++)
			{
				for (Myint64 i2 = i2Start; i2<= i2End; i2++)
				{
					for (Myint64 i1 = i1Start; i1<= i1End; i1++)
					{
						w[i1+i2*ln1+i3*ln1*ln2] =
								FD_D2_O4_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, ln1, ln2, ln3) ;
					}
				}
			}
		}
		else if (fdOrder == 6)
		{
#pragma omp parallel for collapse(2)
			for (Myint64 i3 = i3Start; i3<= i3End; i3++)
			{
				for (Myint64 i2 = i2Start; i2<= i2End; i2++)
				{
					for (Myint64 i1 = i1Start; i1<= i1End; i1++)
					{
						w[i1+i2*ln1+i3*ln1*ln2] =
								FD_D2_O6_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, ln1, ln2, ln3) ;
					}
				}
			}
		}
		else if (fdOrder == 8)
		{
#pragma omp parallel for collapse(2)
			for (Myint64 i3 = i3Start; i3<= i3End; i3++)
			{
				for (Myint64 i2 = i2Start; i2<= i2End; i2++)
				{
					for (Myint64 i1 = i1Start; i1<= i1End; i1++)
					{
						w[i1+i2*ln1+i3*ln1*ln2] =
								FD_D2_O8_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, ln1, ln2, ln3) ;
					}
				}
			}
		}
		else if (fdOrder == 10)
		{
#pragma omp parallel for collapse(2)
			for (Myint64 i3 = i3Start; i3<= i3End; i3++)
			{
				for (Myint64 i2 = i2Start; i2<= i2End; i2++)
				{
					for (Myint64 i1 = i1Start; i1<= i1End; i1++)
					{
						w[i1+i2*ln1+i3*ln1*ln2] =
								FD_D2_O10_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, ln1, ln2, ln3) ;
					}
				}
			}
		}
		else if (fdOrder == 12)
		{
#pragma omp parallel for collapse(2)
			for (Myint64 i3 = i3Start; i3<= i3End; i3++)
			{
				for (Myint64 i2 = i2Start; i2<= i2End; i2++)
				{
					for (Myint64 i1 = i1Start; i1<= i1End; i1++)
					{
						w[i1+i2*ln1+i3*ln1*ln2] =
								FD_D2_O12_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, ln1, ln2, ln3) ;
					}
				}
			}
		}
		else if (fdOrder == 14)
		{
#pragma omp parallel for collapse(2)
			for (Myint64 i3 = i3Start; i3<= i3End; i3++)
			{
				for (Myint64 i2 = i2Start; i2<= i2End; i2++)
				{
					for (Myint64 i1 = i1Start; i1<= i1End; i1++)
					{
						w[i1+i2*ln1+i3*ln1*ln2] =
								FD_D2_O14_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, ln1, ln2, ln3) ;
					}
				}
			}
		}
		else if (fdOrder == 16)
		{
#pragma omp parallel for collapse(2)
			for (Myint64 i3 = i3Start; i3<= i3End; i3++)
			{
				for (Myint64 i2 = i2Start; i2<= i2End; i2++)
				{
					for (Myint64 i1 = i1Start; i1<= i1End; i1++)
					{
						w[i1+i2*ln1+i3*ln1*ln2] =
								FD_D2_O16_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, ln1, ln2, ln3) ;
					}
				}
			}
		}
	}

	printDebug(MID_DEBUG, "OUT Grid_NEC::FD_D2_N1");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_NEC::FD_D2_N2(Point_type pType, const Grid& Wgrid, Myint fdOrder)
{
	printDebug(MID_DEBUG, "IN Grid_NEC::FD_D2_N2");

	// check grids are same size
	if (this->sameSize(Wgrid) != true)
	{
		printError("Grid_NEC::FD_D2_N2, grids have not same size") ;
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

	// compute FD along N2
#ifndef _DOUBLE_PRECISION_
	// if (flag_packed_stencil)
	if (flag_packed_stencil)
	{
		if (fdOrder == 2)
		{
#pragma omp parallel for
			for (Myint64 i3 = i3Start; i3<= i3End; i3++)
			{
#pragma _NEC outerloop_unroll(16)
				for (Myint64 i2 = i2Start; i2<= i2End; i2++)
				{
#pragma _NEC packed_stencil
					for (Myint64 i1 = i1Start; i1<= i1End; i1++)
					{
						w[i1+i2*n1+i3*n1*n2] =
								FD_D2_O2_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
					}
				}
			}
		}
		else if (fdOrder == 4)
		{
#pragma omp parallel for
			for (Myint64 i3 = i3Start; i3<= i3End; i3++)
			{
#pragma _NEC outerloop_unroll(16)
				for (Myint64 i2 = i2Start; i2<= i2End; i2++)
				{
#pragma _NEC packed_stencil
					for (Myint64 i1 = i1Start; i1<= i1End; i1++)
					{
						w[i1+i2*n1+i3*n1*n2] =
								FD_D2_O4_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
					}
				}
			}
		}
		else if (fdOrder == 6)
		{
#pragma omp parallel for
			for (Myint64 i3 = i3Start; i3<= i3End; i3++)
			{
#pragma _NEC outerloop_unroll(16)
				for (Myint64 i2 = i2Start; i2<= i2End; i2++)
				{
#pragma _NEC packed_stencil
					for (Myint64 i1 = i1Start; i1<= i1End; i1++)
					{
						w[i1+i2*n1+i3*n1*n2] =
								FD_D2_O6_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
					}
				}
			}
		}
		else if (fdOrder == 8)
		{
#pragma omp parallel for
			for (Myint64 i3 = i3Start; i3<= i3End; i3++)
			{
#pragma _NEC outerloop_unroll(16)
				for (Myint64 i2 = i2Start; i2<= i2End; i2++)
				{
#pragma _NEC packed_stencil
					for (Myint64 i1 = i1Start; i1<= i1End; i1++)
					{
						w[i1+i2*n1+i3*n1*n2] =
								FD_D2_O8_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
					}
				}
			}
		}
		else if (fdOrder == 10)
		{
#pragma omp parallel for
			for (Myint64 i3 = i3Start; i3<= i3End; i3++)
			{
#pragma _NEC outerloop_unroll(16)
				for (Myint64 i2 = i2Start; i2<= i2End; i2++)
				{
#pragma _NEC packed_stencil
					for (Myint64 i1 = i1Start; i1<= i1End; i1++)
					{
						w[i1+i2*n1+i3*n1*n2] =
								FD_D2_O10_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
					}
				}
			}
		}
		else if (fdOrder == 12)
		{
#pragma omp parallel for
			for (Myint64 i3 = i3Start; i3<= i3End; i3++)
			{
#pragma _NEC outerloop_unroll(16)
				for (Myint64 i2 = i2Start; i2<= i2End; i2++)
				{
#pragma _NEC packed_stencil
					for (Myint64 i1 = i1Start; i1<= i1End; i1++)
					{
						w[i1+i2*n1+i3*n1*n2] =
								FD_D2_O12_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
					}
				}
			}
		}
		else if (fdOrder == 14)
		{
#pragma omp parallel for
			for (Myint64 i3 = i3Start; i3<= i3End; i3++)
			{
#pragma _NEC outerloop_unroll(16)
				for (Myint64 i2 = i2Start; i2<= i2End; i2++)
				{
#pragma _NEC packed_stencil
					for (Myint64 i1 = i1Start; i1<= i1End; i1++)
					{
						w[i1+i2*n1+i3*n1*n2] =
								FD_D2_O14_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
					}
				}
			}
		}
		else if (fdOrder == 16)
		{
#pragma omp parallel for
			for (Myint64 i3 = i3Start; i3<= i3End; i3++)
			{
#pragma _NEC outerloop_unroll(16)
				for (Myint64 i2 = i2Start; i2<= i2End; i2++)
				{
#pragma _NEC packed_stencil
					for (Myint64 i1 = i1Start; i1<= i1End; i1++)
					{
						w[i1+i2*n1+i3*n1*n2] =
								FD_D2_O16_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
					}
				}
			}
		}
	}
	// else if (flag_packed_stencil)
	else
#endif
	{
		if (fdOrder == 2)
		{
#pragma omp parallel for
			for (Myint64 i3 = i3Start; i3<= i3End; i3++)
			{
#pragma _NEC outerloop_unroll(16)
				for (Myint64 i2 = i2Start; i2<= i2End; i2++)
				{
					for (Myint64 i1 = i1Start; i1<= i1End; i1++)
					{
						w[i1+i2*n1+i3*n1*n2] =
								FD_D2_O2_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
					}
				}
			}
		}
		else if (fdOrder == 4)
		{
#pragma omp parallel for
			for (Myint64 i3 = i3Start; i3<= i3End; i3++)
			{
#pragma _NEC outerloop_unroll(16)
				for (Myint64 i2 = i2Start; i2<= i2End; i2++)
				{
					for (Myint64 i1 = i1Start; i1<= i1End; i1++)
					{
						w[i1+i2*n1+i3*n1*n2] =
								FD_D2_O4_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
					}
				}
			}
		}
		else if (fdOrder == 6)
		{
#pragma omp parallel for
			for (Myint64 i3 = i3Start; i3<= i3End; i3++)
			{
#pragma _NEC outerloop_unroll(16)
				for (Myint64 i2 = i2Start; i2<= i2End; i2++)
				{
					for (Myint64 i1 = i1Start; i1<= i1End; i1++)
					{
						w[i1+i2*n1+i3*n1*n2] =
								FD_D2_O6_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
					}
				}
			}
		}
		else if (fdOrder == 8)
		{
#pragma omp parallel for
			for (Myint64 i3 = i3Start; i3<= i3End; i3++)
			{
#pragma _NEC outerloop_unroll(16)
				for (Myint64 i2 = i2Start; i2<= i2End; i2++)
				{
					for (Myint64 i1 = i1Start; i1<= i1End; i1++)
					{
						w[i1+i2*n1+i3*n1*n2] =
								FD_D2_O8_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
					}
				}
			}
		}
		else if (fdOrder == 10)
		{
#pragma omp parallel for
			for (Myint64 i3 = i3Start; i3<= i3End; i3++)
			{
#pragma _NEC outerloop_unroll(16)
				for (Myint64 i2 = i2Start; i2<= i2End; i2++)
				{
					for (Myint64 i1 = i1Start; i1<= i1End; i1++)
					{
						w[i1+i2*n1+i3*n1*n2] =
								FD_D2_O10_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
					}
				}
			}
		}
		else if (fdOrder == 12)
		{
#pragma omp parallel for
			for (Myint64 i3 = i3Start; i3<= i3End; i3++)
			{
#pragma _NEC outerloop_unroll(16)
				for (Myint64 i2 = i2Start; i2<= i2End; i2++)
				{
					for (Myint64 i1 = i1Start; i1<= i1End; i1++)
					{
						w[i1+i2*n1+i3*n1*n2] =
								FD_D2_O12_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
					}
				}
			}
		}
		else if (fdOrder == 14)
		{
#pragma omp parallel for
			for (Myint64 i3 = i3Start; i3<= i3End; i3++)
			{
#pragma _NEC outerloop_unroll(16)
				for (Myint64 i2 = i2Start; i2<= i2End; i2++)
				{
					for (Myint64 i1 = i1Start; i1<= i1End; i1++)
					{
						w[i1+i2*n1+i3*n1*n2] =
								FD_D2_O14_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
					}
				}
			}
		}
		else if (fdOrder == 16)
		{
#pragma omp parallel for
			for (Myint64 i3 = i3Start; i3<= i3End; i3++)
			{
#pragma _NEC outerloop_unroll(16)
				for (Myint64 i2 = i2Start; i2<= i2End; i2++)
				{
					for (Myint64 i1 = i1Start; i1<= i1End; i1++)
					{
						w[i1+i2*n1+i3*n1*n2] =
								FD_D2_O16_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
					}
				}
			}
		}
	}

	printDebug(MID_DEBUG, "OUT Grid_NEC::FD_D2_N2");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_NEC::FD_D2_N3(Point_type pType, const Grid& Wgrid, Myint fdOrder)
{
	printDebug(MID_DEBUG, "IN Grid_NEC::FD_D2_N3");

	// check grids are same size
	if (this->sameSize(Wgrid) != true)
	{
		printError("Grid_NEC::FD_D2_N3, grids have not same size") ;
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

	// compute FD along N3
#ifndef _DOUBLE_PRECISION_
	// if (flag_packed_stencil)
	if (flag_packed_stencil)
	{
		if (fdOrder == 2)
		{
#pragma omp parallel for
			for (Myint64 i2 = i2Start; i2<= i2End; i2++)
			{
#pragma _NEC outerloop_unroll(16)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
#pragma _NEC packed_stencil
					for (Myint64 i1 = i1Start; i1<= i1End; i1++)
					{
						w[i1+i2*n1+i3*n1*n2] =
								FD_D2_O2_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
					}
				}
			}
		}
		else if (fdOrder == 4)
		{
#pragma omp parallel for
			for (Myint64 i2 = i2Start; i2<= i2End; i2++)
			{
#pragma _NEC outerloop_unroll(16)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
#pragma _NEC packed_stencil
					for (Myint64 i1 = i1Start; i1<= i1End; i1++)
					{
						w[i1+i2*n1+i3*n1*n2] =
								FD_D2_O4_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
					}
				}
			}
		}
		else if (fdOrder == 6)
		{
#pragma omp parallel for
			for (Myint64 i2 = i2Start; i2<= i2End; i2++)
			{
#pragma _NEC outerloop_unroll(16)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
#pragma _NEC packed_stencil
					for (Myint64 i1 = i1Start; i1<= i1End; i1++)
					{
						w[i1+i2*n1+i3*n1*n2] =
								FD_D2_O6_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
					}
				}
			}
		}
		else if (fdOrder == 8)
		{
#pragma omp parallel for
			for (Myint64 i2 = i2Start; i2<= i2End; i2++)
			{
#pragma _NEC outerloop_unroll(16)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
#pragma _NEC packed_stencil
					for (Myint64 i1 = i1Start; i1<= i1End; i1++)
					{
						w[i1+i2*n1+i3*n1*n2] =
								FD_D2_O8_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
					}
				}
			}
		}
		else if (fdOrder == 10)
		{
#pragma omp parallel for
			for (Myint64 i2 = i2Start; i2<= i2End; i2++)
			{
#pragma _NEC outerloop_unroll(16)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
#pragma _NEC packed_stencil
					for (Myint64 i1 = i1Start; i1<= i1End; i1++)
					{
						w[i1+i2*n1+i3*n1*n2] =
								FD_D2_O10_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
					}
				}
			}
		}
		else if (fdOrder == 12)
		{
#pragma omp parallel for
			for (Myint64 i2 = i2Start; i2<= i2End; i2++)
			{
#pragma _NEC outerloop_unroll(16)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
#pragma _NEC packed_stencil
					for (Myint64 i1 = i1Start; i1<= i1End; i1++)
					{
						w[i1+i2*n1+i3*n1*n2] =
								FD_D2_O12_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
					}
				}
			}
		}
		else if (fdOrder == 14)
		{
#pragma omp parallel for
			for (Myint64 i2 = i2Start; i2<= i2End; i2++)
			{
#pragma _NEC outerloop_unroll(16)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
#pragma _NEC packed_stencil
					for (Myint64 i1 = i1Start; i1<= i1End; i1++)
					{
						w[i1+i2*n1+i3*n1*n2] =
								FD_D2_O14_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
					}
				}
			}
		}
		else if (fdOrder == 16)
		{
#pragma omp parallel for
			for (Myint64 i2 = i2Start; i2<= i2End; i2++)
			{
#pragma _NEC outerloop_unroll(16)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
#pragma _NEC packed_stencil
					for (Myint64 i1 = i1Start; i1<= i1End; i1++)
					{
						w[i1+i2*n1+i3*n1*n2] =
								FD_D2_O16_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
					}

				}
			}
		}
	}
	// else if (flag_packed_stencil)
	else
#endif
	{
		if (fdOrder == 2)
		{
#pragma omp parallel for
			for (Myint64 i2 = i2Start; i2<= i2End; i2++)
			{
#pragma _NEC outerloop_unroll(16)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i1 = i1Start; i1<= i1End; i1++)
					{
						w[i1+i2*n1+i3*n1*n2] =
								FD_D2_O2_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
					}
				}
			}
		}
		else if (fdOrder == 4)
		{
#pragma omp parallel for
			for (Myint64 i2 = i2Start; i2<= i2End; i2++)
			{
#pragma _NEC outerloop_unroll(16)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i1 = i1Start; i1<= i1End; i1++)
					{
						w[i1+i2*n1+i3*n1*n2] =
								FD_D2_O4_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
					}
				}
			}
		}
		else if (fdOrder == 6)
		{
#pragma omp parallel for
			for (Myint64 i2 = i2Start; i2<= i2End; i2++)
			{
#pragma _NEC outerloop_unroll(16)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i1 = i1Start; i1<= i1End; i1++)
					{
						w[i1+i2*n1+i3*n1*n2] =
								FD_D2_O6_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
					}
				}
			}
		}
		else if (fdOrder == 8)
		{
#pragma omp parallel for
			for (Myint64 i2 = i2Start; i2<= i2End; i2++)
			{
#pragma _NEC outerloop_unroll(16)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i1 = i1Start; i1<= i1End; i1++)
					{
						w[i1+i2*n1+i3*n1*n2] =
								FD_D2_O8_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
					}
				}
			}
		}
		else if (fdOrder == 10)
		{
#pragma omp parallel for
			for (Myint64 i2 = i2Start; i2<= i2End; i2++)
			{
#pragma _NEC outerloop_unroll(16)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i1 = i1Start; i1<= i1End; i1++)
					{
						w[i1+i2*n1+i3*n1*n2] =
								FD_D2_O10_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
					}
				}
			}
		}
		else if (fdOrder == 12)
		{
#pragma omp parallel for
			for (Myint64 i2 = i2Start; i2<= i2End; i2++)
			{
#pragma _NEC outerloop_unroll(16)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i1 = i1Start; i1<= i1End; i1++)
					{
						w[i1+i2*n1+i3*n1*n2] =
								FD_D2_O12_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
					}
				}
			}
		}
		else if (fdOrder == 14)
		{
#pragma omp parallel for
			for (Myint64 i2 = i2Start; i2<= i2End; i2++)
			{
#pragma _NEC outerloop_unroll(16)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i1 = i1Start; i1<= i1End; i1++)
					{
						w[i1+i2*n1+i3*n1*n2] =
								FD_D2_O14_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
					}
				}
			}
		}
		else if (fdOrder == 16)
		{
#pragma omp parallel for
			for (Myint64 i2 = i2Start; i2<= i2End; i2++)
			{
#pragma _NEC outerloop_unroll(16)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i1 = i1Start; i1<= i1End; i1++)
					{
						w[i1+i2*n1+i3*n1*n2] =
								FD_D2_O16_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
					}

				}
			}
		}
	}

	printDebug(MID_DEBUG, "OUT Grid_NEC::FD_D2_N3");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_NEC::FD_LAPLACIAN(Point_type pType, const Grid& Wgrid, Myint fdOrder)
{
	printDebug(MID_DEBUG, "IN Grid_NEC::FD_LAPLACIAN");

	// check grids are same size
	if (this->sameSize(Wgrid) != true)
	{
		printError("Grid_NEC::FD_LAPLACIAN, grids have not same size") ;
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

#ifndef _DOUBLE_PRECISION_
	// if (flag_packed_stencil)
	if (flag_packed_stencil)
	{
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
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma _NEC packed_stencil
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							w[i1+i2*n1+i3*n1*n2] =
									FD_D2_O2_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O2_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
						}
					}
				}
			}
			else if (fdOrder == 4)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma _NEC packed_stencil
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							w[i1+i2*n1+i3*n1*n2] =
									FD_D2_O4_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O4_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
						}
					}
				}
			}
			else if (fdOrder == 6)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma _NEC packed_stencil
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							w[i1+i2*n1+i3*n1*n2] =
									FD_D2_O4_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O6_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
						}
					}
				}
			}
			else if (fdOrder == 8)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma _NEC packed_stencil
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							w[i1+i2*n1+i3*n1*n2] =
									FD_D2_O8_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O8_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
						}
					}
				}
			}
			else if (fdOrder == 10)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma _NEC packed_stencil
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							w[i1+i2*n1+i3*n1*n2] =
									FD_D2_O10_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O10_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
						}
					}
				}
			}
			else if (fdOrder == 12)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma _NEC packed_stencil
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							w[i1+i2*n1+i3*n1*n2] =
									FD_D2_O12_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O12_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
						}
					}
				}
			}
			else if (fdOrder == 14)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma _NEC packed_stencil
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							w[i1+i2*n1+i3*n1*n2] =
									FD_D2_O14_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O14_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
						}
					}
				}
			}
			else if (fdOrder == 16)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma _NEC packed_stencil
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							w[i1+i2*n1+i3*n1*n2] =
									FD_D2_O16_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O16_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
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
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma _NEC packed_stencil
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							w[i1+i2*n1+i3*n1*n2] =
									FD_D2_O2_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O2_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O2_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
						}
					}
				}
			}
			else if (fdOrder == 4)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma _NEC packed_stencil
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							w[i1+i2*n1+i3*n1*n2] =
									FD_D2_O4_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O4_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O4_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
						}
					}
				}
			}
			else if (fdOrder == 6)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma _NEC packed_stencil
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							w[i1+i2*n1+i3*n1*n2] =
									FD_D2_O6_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O6_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O6_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
						}
					}
				}
			}
			else if (fdOrder == 8)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma _NEC packed_stencil
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							w[i1+i2*n1+i3*n1*n2] =
									FD_D2_O8_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O8_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O8_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
						}
					}
				}
			}
			else if (fdOrder == 10)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma _NEC packed_stencil
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							w[i1+i2*n1+i3*n1*n2] =
									FD_D2_O10_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O10_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O10_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
						}
					}
				}
			}
			else if (fdOrder == 12)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma _NEC packed_stencil
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							w[i1+i2*n1+i3*n1*n2] =
									FD_D2_O12_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O12_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O12_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
						}
					}
				}
			}
			else if (fdOrder == 14)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma _NEC packed_stencil
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							w[i1+i2*n1+i3*n1*n2] =
									FD_D2_O14_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O14_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O14_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
						}
					}
				}
			}
			else if (fdOrder == 16)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma _NEC packed_stencil
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
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
	// else if (flag_packed_stencil)
	else
#endif
	{
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
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							w[i1+i2*n1+i3*n1*n2] =
									FD_D2_O2_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O2_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
						}
					}
				}
			}
			else if (fdOrder == 4)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							w[i1+i2*n1+i3*n1*n2] =
									FD_D2_O4_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O4_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
						}
					}
				}
			}
			else if (fdOrder == 6)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							w[i1+i2*n1+i3*n1*n2] =
									FD_D2_O6_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O6_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
						}
					}
				}
			}
			else if (fdOrder == 8)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							w[i1+i2*n1+i3*n1*n2] =
									FD_D2_O8_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O8_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
						}
					}
				}
			}
			else if (fdOrder == 10)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							w[i1+i2*n1+i3*n1*n2] =
									FD_D2_O10_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O10_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
						}
					}
				}
			}
			else if (fdOrder == 12)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							w[i1+i2*n1+i3*n1*n2] =
									FD_D2_O12_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O12_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
						}
					}
				}
			}
			else if (fdOrder == 14)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							w[i1+i2*n1+i3*n1*n2] =
									FD_D2_O14_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O14_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
						}
					}
				}
			}
			else if (fdOrder == 16)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							w[i1+i2*n1+i3*n1*n2] =
									FD_D2_O16_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O16_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
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
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							w[i1+i2*n1+i3*n1*n2] =
									FD_D2_O2_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O2_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O2_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
						}
					}
				}
			}
			else if (fdOrder == 4)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							w[i1+i2*n1+i3*n1*n2] =
									FD_D2_O4_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O4_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O4_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
						}
					}
				}
			}
			else if (fdOrder == 6)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							w[i1+i2*n1+i3*n1*n2] =
									FD_D2_O6_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O6_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O6_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
						}
					}
				}
			}
			else if (fdOrder == 8)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							w[i1+i2*n1+i3*n1*n2] =
									FD_D2_O8_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O8_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O8_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
						}
					}
				}
			}
			else if (fdOrder == 10)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							w[i1+i2*n1+i3*n1*n2] =
									FD_D2_O10_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O10_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O10_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
						}
					}
				}
			}
			else if (fdOrder == 12)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							w[i1+i2*n1+i3*n1*n2] =
									FD_D2_O12_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O12_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O12_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
						}
					}
				}
			}
			else if (fdOrder == 14)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							w[i1+i2*n1+i3*n1*n2] =
									FD_D2_O14_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O14_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
									+ FD_D2_O14_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
						}
					}
				}
			}
			else if (fdOrder == 16)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
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

	printDebug(MID_DEBUG, "OUT Grid_NEC::FD_LAPLACIAN");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Myfloat Grid_NEC::getMin(Point_type pointType)
{
	printDebug(LIGHT_DEBUG, "IN Grid_NEC::getMin");

	Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
	getGridIndex(pointType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End) ;
	Myfloat val = FLT_MAX ;

	Myint64 nn1 = i1End - i1Start + 1;
	Myint64 nn2 = i2End - i2Start + 1;
	Myint64 nn3 = i3End - i3Start + 1;
	Myint64 n = nn1 * nn2 * nn3;

#pragma omp parallel for
	for (Myint64 i3 = i3Start; i3<= i3End; i3++)
	{
		Myint64 tmp_idx = (i3 - i3Start) * nn1 * nn2;
		for (Myint64 i2 = i2Start; i2<= i2End; i2++)
		{
			Myint64 idx = i1Start + i2 * n1 + i3 * n1 * n2;
			for (Myint64 i1 = i1Start; i1<= i1End; i1++)
			{
				tmp_grid_3d[tmp_idx++] = grid_3d[idx++];
			}
		}
	}

#pragma omp parallel for reduction(min:val)
	for (Myint64 i = 0; i < n; i++)
	{
		if (tmp_grid_3d[i] < val) val = tmp_grid_3d[i] ;
	}

	printDebug(LIGHT_DEBUG, "Min val", val);

	printDebug(LIGHT_DEBUG, "OUT Grid_NEC::getMin");
	return(val) ;
}

//-------------------------------------------------------------------------------------------------------

Myfloat Grid_NEC::getMax(Point_type pointType)
{
	printDebug(LIGHT_DEBUG, "IN Grid_NEC::getMax");

	Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
	getGridIndex(pointType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End) ;

	Myfloat val = -FLT_MAX ;

	Myint64 nn1 = i1End - i1Start + 1;
	Myint64 nn2 = i2End - i2Start + 1;
	Myint64 nn3 = i3End - i3Start + 1;
	Myint64 n = nn1 * nn2 * nn3;

#pragma omp parallel for
	for (Myint64 i3 = i3Start; i3<= i3End; i3++)
	{
		Myint64 tmp_idx = (i3 - i3Start) * nn1 * nn2;
		for (Myint64 i2 = i2Start; i2<= i2End; i2++)
		{
			Myint64 idx = i1Start + i2 * n1 + i3 * n1 * n2;
			for (Myint64 i1 = i1Start; i1<= i1End; i1++)
			{
				tmp_grid_3d[tmp_idx++] = grid_3d[idx++];
			}
		}
	}

#pragma omp parallel for reduction(max:val)
	for (Myint64 i = 0; i < n; i++)
	{
		if (tmp_grid_3d[i] > val) val = tmp_grid_3d[i] ;
	}

	printDebug(LIGHT_DEBUG, "Max val", val);

	printDebug(LIGHT_DEBUG, "OUT Grid_NEC::getMax");
	return(val) ;
}

//-------------------------------------------------------------------------------------------------------

Myfloat Grid_NEC::maxErr(Point_type pointType, const Grid& gridIn) const
{
	printDebug(FULL_DEBUG, "IN Grid_NEC::maxErr");

	// check grids have same size
	if (!(this->sameSize(gridIn)))
	{
		printError("Grid_NEC::maxErr, grids have different size") ;
		return(-1.0) ;
	}

	Myfloat err = -FLT_MAX, err2 = 0.0 ;

	Myfloat* u1 = this->grid_3d ;
	Myfloat* u2 = gridIn.grid_3d ;

	Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
	getGridIndex(pointType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End) ;

#pragma omp parallel for reduction(max:err)
	for (Myint64 i3 = i3Start; i3<= i3End; i3++)
	{
		for (Myint64 i2 = i1Start + i2Start * n1; i2<= i1End + i2End * n1; i2++)
		{
			Myint64 ii = i2 + i3*n2*n1 ;
			Myint64 i1 = (i2 - i1Start) % n1 + i1Start;
			// prevent divide by 0
			if (fabs(u2[ii]) < MAX_ERR_FLOAT)
			{
				err2 = fabs(u1[ii] - u2[ii]) ;
			}
			else
			{
				err2 = fabs(u1[ii] - u2[ii]) / u2[ii] ;
			}

			if (i1 >= i1Start && i1 <= i1End && err2 > err)
			{
				err = err2 ;
			}
		}
	}

	printDebug(FULL_DEBUG, "OUT Grid_NEC::maxErr");
	return(err) ;
}

//-------------------------------------------------------------------------------------------------------

void Grid_NEC::padGridn1(void)
{
	printDebug(MID_DEBUG, "IN Grid_NEC::padGridn1");

	if (Config::Instance()->autoPad == true)
	{
		sca_int_t tmp_n1, tmp_n2, tmp_n3, m1, m2, m3;
		tmp_n1 = i1Halo2End + 1 ;
		tmp_n2 = 10;
		tmp_n3 = 10;
		sca_utility_optimize_leading(tmp_n1, tmp_n2, tmp_n3, 1, &m1, &m2, &m3);
		if (tmp_n1 != m1)
		{
			i1PadStart = i1Halo2End + 1 ;
			i1PadEnd   = i1PadStart + (m1 - tmp_n1) - 1 ;
		}
		else
		{
			i1PadStart = i1Halo2End ;
			i1PadEnd   = i1Halo2End ;
		}
	}
	else if (Config::Instance()->n1AddPad != UNSPECIFIED)
	{
		i1PadStart = i1Halo2End + 1 ;
		i1PadEnd   = i1PadStart + Config::Instance()->n1AddPad - 1 ;
	}
	else if (Config::Instance()->n1MulPad != UNSPECIFIED)
	{
		i1PadStart = i1Halo2End + 1 ;
		Myint nTmp = (i1PadStart + 1) / Config::Instance()->n1MulPad ;
		if ((Config::Instance()->n1MulPad * nTmp) < (i1PadStart))
		{
			Myint nTot = (nTmp + 1) * Config::Instance()->n1MulPad ;
			i1PadEnd = nTot - 1 ;
		}
		else
		{
			Myint nTot = (nTmp) * Config::Instance()->n1MulPad ;
			i1PadEnd = nTot - 1 ;
		}
	}
	else
	{
		i1PadStart = i1Halo2End ;
		i1PadEnd   = i1Halo2End ;
	}

	printDebug(MID_DEBUG, "OUT Grid_NEC::padGridn1");
}

//-------------------------------------------------------------------------------------------------------

void Grid_NEC::padGridn2(void)
{
	printDebug(MID_DEBUG, "IN Grid_NEC::padGridn2");

	if (Config::Instance()->autoPad == true)
	{
		sca_int_t tmp_n2, tmp_n3, m1, m2, m3;
		tmp_n2 = i2Halo2End + 1 ;
		tmp_n3 = 10;
		sca_utility_optimize_leading(n1, tmp_n2, tmp_n3, 1, &m1, &m2, &m3);
		if (tmp_n2 != m2)
		{
			i2PadStart = i2Halo2End + 1 ;
			i2PadEnd   = i2PadStart + (m2 - tmp_n2) - 1 ;
		}
		else
		{
			i2PadStart = i2Halo2End ;
			i2PadEnd   = i2Halo2End ;
		}
	}
	else if (Config::Instance()->n2AddPad != UNSPECIFIED)
	{
		i2PadStart = i2Halo2End + 1 ;
		i2PadEnd   = i2PadStart + Config::Instance()->n2AddPad - 1 ;
	}
	else if (Config::Instance()->n2MulPad != UNSPECIFIED)
	{
		i2PadStart = i2Halo2End + 1 ;
		Myint nTmp = (i2PadStart + 1) / Config::Instance()->n2MulPad ;
		if ((Config::Instance()->n2MulPad * nTmp) < (i2PadStart))
		{
			Myint nTot = (nTmp + 1) * Config::Instance()->n2MulPad ;
			i2PadEnd = nTot - 1 ;
		}
		else
		{
			Myint nTot = (nTmp) * Config::Instance()->n2MulPad ;
			i2PadEnd = nTot - 1 ;
		}
	}
	else
	{
		i2PadStart = i2Halo2End ;
		i2PadEnd   = i2Halo2End ;
	}

	printDebug(MID_DEBUG, "OUT Grid_NEC::padGridn2");
}

//-------------------------------------------------------------------------------------------------------

void Grid_NEC::padGridn3(void)
{
	printDebug(MID_DEBUG, "IN Grid_NEC::padGridn3");

	if (Config::Instance()->autoPad == true)
	{
		sca_int_t tmp_n3, m1, m2, m3;
		tmp_n3 = i3Halo2End + 1 ;
		sca_utility_optimize_leading(n1, n2, tmp_n3, 1, &m1, &m2, &m3);
		if (tmp_n3 != m3)
		{
			i3PadStart = i3Halo2End + 1 ;
			i3PadEnd   = i3PadStart + (m3 - tmp_n3) - 1 ;
		}
		else
		{
			i3PadStart = i3Halo2End ;
			i3PadEnd   = i3Halo2End ;
		}
	}
	else if (Config::Instance()->n3AddPad != UNSPECIFIED)
	{
		i3PadStart = i3Halo2End + 1 ;
		i3PadEnd   = i3PadStart + Config::Instance()->n3AddPad - 1 ;
	}
	else if (Config::Instance()->n3MulPad != UNSPECIFIED)
	{
		i3PadStart = i3Halo2End + 1 ;
		Myint nTmp = (i3PadStart + 1) / Config::Instance()->n3MulPad ;
		if ((Config::Instance()->n3MulPad * nTmp) < (i3PadStart))
		{
			Myint nTot = (nTmp + 1) * Config::Instance()->n3MulPad ;
			i3PadEnd = nTot - 1 ;
		}
		else
		{
			Myint nTot = (nTmp) * Config::Instance()->n3MulPad ;
			i3PadEnd = nTot - 1 ;
		}
	}
	else
	{
		i3PadStart = i3Halo2End ;
		i3PadEnd   = i3Halo2End ;
	}

	printDebug(MID_DEBUG, "OUT Grid_NEC::padGridn3");
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_NEC::initializeGrid(void)
{
	printDebug(MID_DEBUG, "IN Grid_NEC::initializeGrid");

	Grid::initializeGrid() ;

#ifndef _DOUBLE_PRECISION_
	Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
	Myfloat * const u = this->grid_3d ;
	Myint nlayer = Config::Instance()->nlayer ;	
	getGridIndex(INNER_POINTS, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End) ;

	/*
	 * Packed stencil can be used only if:
	 * . first accessed point is 8 bytes aligned.
	 * . loop length is even number.
	 *
	 * Assume all arrays have the same shape.
	 */
	flag_packed_stencil = (((reinterpret_cast<std::uintptr_t>(&(u[i1Start+i2Start*n1+i3Start*n1*n2]))) & 0x7) == 0) && ((i1End-i1Start+1) % 2 == 0) && (nlayer % 2 == 0);
#endif

	tmp_grid_3d = new Myfloat[npoint] ;

	printDebug(MID_DEBUG, "OUT Grid_NEC::initializeGrid");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_NEC::applyBoundaryCondition(BoundCond_type boundCondType)
{
	printDebug(FULL_DEBUG, "IN Grid_NEC::applyBoundaryCondition");

	if (boundCondType == NO_BOUND_COND)
	{
		// nothing to do
	}

	else if (boundCondType == BOUND_COND_ANTI_MIRROR)
	{
		// anti-mirroring value in halos
		Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;

		if (dim >= DIM1)
		{
			// I1HALO1
			if (getNeighbourProc(I1HALO1) == MPI_PROC_NULL)
			{
				getGridIndex(I1HALO1, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End) ;
				Myint64 iInner1 = i1End+1 ;
#pragma omp parallel
				{
#pragma omp for
					for (Myint64 i3 = i3Start; i3<= i3End; i3++)
					{
						for (Myint64 i2 = i2Start; i2<= i2End; i2++)
						{
							// set inner point to 0
							grid_3d[iInner1+i2*n1+i3*n1*n2] = 0.0 ;
						}
					}
#pragma omp for
					for (Myint64 i3 = i3Start; i3<= i3End; i3++)
					{
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
#pragma _NEC ivdep
							for (Myint64 i2 = i2Start; i2<= i2End; i2++)
							{
								// set symetrical point to minus inner point value
								grid_3d[i1+i2*n1+i3*n1*n2] = -grid_3d[(iInner1+iInner1-i1)+i2*n1+i3*n1*n2] ;
							}
						}
					}
				}
			}

			// I1HALO2
			if (getNeighbourProc(I1HALO2) == MPI_PROC_NULL)
			{
				getGridIndex(I1HALO2, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End) ;
				Myint64 iInner1 = i1Start-1 ;
#pragma omp parallel
				{
#pragma omp for
					for (Myint64 i3 = i3Start; i3<= i3End; i3++)
					{
						for (Myint64 i2 = i2Start; i2<= i2End; i2++)
						{
							// set inner point to 0
							grid_3d[iInner1+i2*n1+i3*n1*n2] = 0.0 ;
						}
					}
#pragma omp for
					for (Myint64 i3 = i3Start; i3<= i3End; i3++)
					{
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
#pragma _NEC ivdep
							for (Myint64 i2 = i2Start; i2<= i2End; i2++)
							{
								// set symetrical point to minus inner point value
								grid_3d[i1+i2*n1+i3*n1*n2] = -grid_3d[(iInner1-(i1-iInner1))+i2*n1+i3*n1*n2] ;
							}
						}
					}
				}
			}
		}

		if (dim >= DIM2)
		{
			// I2HALO1
			if (getNeighbourProc(I2HALO1) == MPI_PROC_NULL)
			{
				getGridIndex(I2HALO1, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End) ;
				Myint64 iInner2 = i2End+1 ;
#pragma omp parallel
				{
#pragma omp for
					for (Myint64 i3 = i3Start; i3<= i3End; i3++)
					{
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							// set inner point to 0
							grid_3d[i1+iInner2*n1+i3*n1*n2] = 0.0 ;
						}
					}
#pragma omp for
					for (Myint64 i3 = i3Start; i3<= i3End; i3++)
					{
						for (Myint64 i2 = i2Start; i2<= i2End; i2++)
						{
#pragma _NEC ivdep
							for (Myint64 i1 = i1Start; i1<= i1End; i1++)
							{
								// set symetrical point to minus inner point value
								grid_3d[i1+i2*n1+i3*n1*n2] = -grid_3d[i1+(iInner2+iInner2-i2)*n1+i3*n1*n2] ;
							}
						}
					}
				}
			}

			// I2HALO2
			if (getNeighbourProc(I2HALO2) == MPI_PROC_NULL)
			{
				getGridIndex(I2HALO2, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End) ;
				Myint64 iInner2 = i2Start-1 ;
#pragma omp parallel
				{
#pragma omp for
					for (Myint64 i3 = i3Start; i3<= i3End; i3++)
					{
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							// set inner point to 0
							grid_3d[i1+iInner2*n1+i3*n1*n2] = 0.0 ;
						}
					}
#pragma omp for
					for (Myint64 i3 = i3Start; i3<= i3End; i3++)
					{
						for (Myint64 i2 = i2Start; i2<= i2End; i2++)
						{
#pragma _NEC ivdep
							for (Myint64 i1 = i1Start; i1<= i1End; i1++)
							{
								// set symetrical point to minus inner point value
								grid_3d[i1+i2*n1+i3*n1*n2] = -grid_3d[i1+(iInner2-(i2-iInner2))*n1+i3*n1*n2] ;
							}
						}
					}
				}
			}
		}

		if (dim >= DIM3)
		{
			// I3HALO1
			if (getNeighbourProc(I3HALO1) == MPI_PROC_NULL)
			{
				getGridIndex(I3HALO1, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End) ;
				Myint64 iInner3 = i3End+1 ;
#pragma omp parallel
				{
#pragma omp for
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							// set inner point to 0
							grid_3d[i1+i2*n1+iInner3*n1*n2] = 0.0 ;
						}
					}
#pragma omp for
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
						for (Myint64 i3 = i3Start; i3<= i3End; i3++)
						{
#pragma _NEC ivdep
							for (Myint64 i1 = i1Start; i1<= i1End; i1++)
							{
								// set symetrical point to minus inner point value
								grid_3d[i1+i2*n1+i3*n1*n2] = -grid_3d[i1+i2*n1+(iInner3+iInner3-i3)*n1*n2] ;
							}
						}
					}
				}
			}

			// I3HALO2
			if (getNeighbourProc(I3HALO2) == MPI_PROC_NULL)
			{
				getGridIndex(I3HALO2, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End) ;
				Myint64 iInner3 = i3Start-1 ;
#pragma omp parallel
				{
#pragma omp for
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							// set inner point to 0
							grid_3d[i1+i2*n1+iInner3*n1*n2] = 0.0 ;
						}
					}
#pragma omp for
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
						for (Myint64 i3 = i3Start; i3<= i3End; i3++)
						{
#pragma _NEC ivdep
							for (Myint64 i1 = i1Start; i1<= i1End; i1++)
							{
								// set symetrical point to minus inner point value
								grid_3d[i1+i2*n1+i3*n1*n2] = -grid_3d[i1+i2*n1+(iInner3-(i3-iInner3))*n1*n2] ;
							}
						}
					}
				}
			}
		}
	}
	else
	{
		printError("IN Grid_NEC::applyBoundaryCondition, invalid boundCondType", boundCondType) ;
		return(RTN_CODE_KO) ;
	}

	printDebug(FULL_DEBUG, "OUT Grid_NEC::applyBoundaryCondition");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_NEC::updatePressure(Point_type pType, const Grid& prcGrid,
		const Grid& coefGrid, const Grid& laplaGrid)
{
	printDebug(MID_DEBUG, "IN Grid_NEC::updatePressure");

	Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
	getGridIndex(INNER_POINTS, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End) ;

	Myfloat * prn   = this->grid_3d ;
	Myfloat * prc   = prcGrid.grid_3d ;
	Myfloat * lapla = laplaGrid.grid_3d ;
	Myfloat * coef  = coefGrid.grid_3d ;

#ifndef _DOUBLE_PRECISION_
	if (flag_packed_stencil)
	{
#pragma omp parallel for
		for (Myint64 i3 = i3Start; i3<= i3End; i3++)
		{
			for (Myint64 i2 = i2Start; i2<= i2End; i2++)
			{
#pragma _NEC packed_stencil
				for (Myint64 i1 = i1Start; i1<= i1End; i1++)
				{
					prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
							coef[i1+i2*n1+i3*n1*n2] * lapla[i1+i2*n1+i3*n1*n2] ;
				}
			}
		}
	}
	else
#endif
	{
#pragma omp parallel for collapse(2)
		for (Myint64 i3 = i3Start; i3<= i3End; i3++)
		{
			for (Myint64 i2 = i2Start; i2<= i2End; i2++)
			{
				for (Myint64 i1 = i1Start; i1<= i1End; i1++)
				{
					prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
							coef[i1+i2*n1+i3*n1*n2] * lapla[i1+i2*n1+i3*n1*n2] ;
				}
			}
		}
	}

	printDebug(MID_DEBUG, "OUT Grid_NEC::updatePressure");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_NEC::computePressureWithFD(Grid& prcGridIn, Grid& coefGridIn, Myint fdOrder)
{

	printDebug(FULL_DEBUG, "In Grid_NEC::computePressureWithFD") ;

	// check grids are same size
	if (this->sameSize(prcGridIn) != true)
	{
		printError("In Grid_NEC::computePressureWithFD, grids have not same size") ;
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

#ifndef _DOUBLE_PRECISION_
	// if (flag_packed_stencil)
	if (flag_packed_stencil)
	{
		// compute FD for 1D
		if (dim == DIM1)
		{
			if (fdOrder == 2)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma _NEC packed_stencil
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
									coef[i1+i2*n1+i3*n1*n2] *
									FD_D2_O2_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
						}

					}
				}
			}
			else if (fdOrder == 4)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma _NEC packed_stencil
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
									coef[i1+i2*n1+i3*n1*n2] *
									FD_D2_O4_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
						}

					}
				}
			}
			else if (fdOrder == 6)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma _NEC packed_stencil
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
									coef[i1+i2*n1+i3*n1*n2] *
									FD_D2_O6_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
						}

					}
				}
			}
			else if (fdOrder == 8)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma _NEC packed_stencil
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
									coef[i1+i2*n1+i3*n1*n2] *
									FD_D2_O8_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
						}

					}
				}
			}
			else if (fdOrder == 10)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma _NEC packed_stencil
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
									coef[i1+i2*n1+i3*n1*n2] *
									FD_D2_O10_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
						}

					}
				}
			}
			else if (fdOrder == 12)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma _NEC packed_stencil
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
									coef[i1+i2*n1+i3*n1*n2] *
									FD_D2_O12_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
						}

					}
				}
			}
			else if (fdOrder == 14)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma _NEC packed_stencil
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
									coef[i1+i2*n1+i3*n1*n2] *
									FD_D2_O14_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
						}

					}
				}
			}
			else if (fdOrder == 16)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma _NEC packed_stencil
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
									coef[i1+i2*n1+i3*n1*n2] *
									FD_D2_O16_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
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
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma _NEC packed_stencil
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
									coef[i1+i2*n1+i3*n1*n2] *
									(FD_D2_O2_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
											+ FD_D2_O2_N2(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)) ;
						}

					}
				}
			}
			else if (fdOrder == 4)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma _NEC packed_stencil
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
									coef[i1+i2*n1+i3*n1*n2] *
									(FD_D2_O4_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
											+ FD_D2_O4_N2(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)) ;
						}

					}
				}
			}
			else if (fdOrder == 6)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma _NEC packed_stencil
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
									coef[i1+i2*n1+i3*n1*n2] *
									(FD_D2_O6_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
											+ FD_D2_O6_N2(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)) ;
						}

					}
				}
			}
			else if (fdOrder == 8)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma _NEC packed_stencil
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
									coef[i1+i2*n1+i3*n1*n2] *
									(FD_D2_O8_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
											+ FD_D2_O8_N2(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)) ;
						}

					}
				}
			}
			else if (fdOrder == 10)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma _NEC packed_stencil
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
									coef[i1+i2*n1+i3*n1*n2] *
									(FD_D2_O10_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
											+ FD_D2_O10_N2(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)) ;
						}

					}
				}
			}
			else if (fdOrder == 12)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma _NEC packed_stencil
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
									coef[i1+i2*n1+i3*n1*n2] *
									(FD_D2_O12_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
											+ FD_D2_O12_N2(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)) ;
						}

					}
				}
			}
			else if (fdOrder == 14)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma _NEC packed_stencil
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
									coef[i1+i2*n1+i3*n1*n2] *
									(FD_D2_O14_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
											+ FD_D2_O14_N2(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)) ;
						}

					}
				}
			}
			else if (fdOrder == 16)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma _NEC packed_stencil
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
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

		// compute FD for 3D
		else if (dim == DIM3)
		{
			if (fdOrder == 2)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma _NEC packed_stencil
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
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
			else if (fdOrder == 4)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma _NEC packed_stencil
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
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
			else if (fdOrder == 6)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma _NEC packed_stencil
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
									coef[i1+i2*n1+i3*n1*n2] *
									(FD_D2_O6_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
											+ FD_D2_O6_N2(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
											+ FD_D2_O6_N3(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)) ;
						}

					}
				}
			}
			else if (fdOrder == 8)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma _NEC packed_stencil
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
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
			else if (fdOrder == 10)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma _NEC packed_stencil
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
									coef[i1+i2*n1+i3*n1*n2] *
									(FD_D2_O10_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
											+ FD_D2_O10_N2(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
											+ FD_D2_O10_N3(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)) ;
						}

					}
				}
			}
			else if (fdOrder == 12)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma _NEC packed_stencil
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
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
			else if (fdOrder == 14)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma _NEC packed_stencil
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
									coef[i1+i2*n1+i3*n1*n2] *
									(FD_D2_O14_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
											+ FD_D2_O14_N2(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
											+ FD_D2_O14_N3(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)) ;
						}

					}
				}
			}
			else if (fdOrder == 16)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma _NEC packed_stencil
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
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
	// else if (flag_packed_stencil)
	else
#endif
	{

		// compute FD for 1D
		if (dim == DIM1)
		{
			if (fdOrder == 2)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma omp simd
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
									coef[i1+i2*n1+i3*n1*n2] *
									FD_D2_O2_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
						}

					}
				}
			}
			else if (fdOrder == 4)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma omp simd
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
									coef[i1+i2*n1+i3*n1*n2] *
									FD_D2_O4_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
						}

					}
				}
			}
			else if (fdOrder == 6)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma omp simd
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
									coef[i1+i2*n1+i3*n1*n2] *
									FD_D2_O6_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
						}

					}
				}
			}
			else if (fdOrder == 8)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma omp simd
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
									coef[i1+i2*n1+i3*n1*n2] *
									FD_D2_O8_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
						}

					}
				}
			}
			else if (fdOrder == 10)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma omp simd
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
									coef[i1+i2*n1+i3*n1*n2] *
									FD_D2_O10_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
						}

					}
				}
			}
			else if (fdOrder == 12)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma omp simd
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
									coef[i1+i2*n1+i3*n1*n2] *
									FD_D2_O12_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
						}

					}
				}
			}
			else if (fdOrder == 14)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma omp simd
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
									coef[i1+i2*n1+i3*n1*n2] *
									FD_D2_O14_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
						}

					}
				}
			}
			else if (fdOrder == 16)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma omp simd
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
									coef[i1+i2*n1+i3*n1*n2] *
									FD_D2_O16_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
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
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma omp simd
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
									coef[i1+i2*n1+i3*n1*n2] *
									(FD_D2_O2_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
											+ FD_D2_O2_N2(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)) ;
						}

					}
				}
			}
			else if (fdOrder == 4)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma omp simd
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
									coef[i1+i2*n1+i3*n1*n2] *
									(FD_D2_O4_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
											+ FD_D2_O4_N2(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)) ;
						}

					}
				}
			}
			else if (fdOrder == 6)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma omp simd
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
									coef[i1+i2*n1+i3*n1*n2] *
									(FD_D2_O6_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
											+ FD_D2_O6_N2(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)) ;
						}

					}
				}
			}
			else if (fdOrder == 8)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma omp simd
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
									coef[i1+i2*n1+i3*n1*n2] *
									(FD_D2_O8_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
											+ FD_D2_O8_N2(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)) ;
						}

					}
				}
			}
			else if (fdOrder == 10)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma omp simd
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
									coef[i1+i2*n1+i3*n1*n2] *
									(FD_D2_O10_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
											+ FD_D2_O10_N2(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)) ;
						}

					}
				}
			}
			else if (fdOrder == 12)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma omp simd
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
									coef[i1+i2*n1+i3*n1*n2] *
									(FD_D2_O12_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
											+ FD_D2_O12_N2(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)) ;
						}

					}
				}
			}
			else if (fdOrder == 14)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma omp simd
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
									coef[i1+i2*n1+i3*n1*n2] *
									(FD_D2_O14_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
											+ FD_D2_O14_N2(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)) ;
						}

					}
				}
			}
			else if (fdOrder == 16)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma omp simd
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
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

		// compute FD for 3D
		else if (dim == DIM3)
		{
			if (fdOrder == 2)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma omp simd
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
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
			else if (fdOrder == 4)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma omp simd
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
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
			else if (fdOrder == 6)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma omp simd
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
									coef[i1+i2*n1+i3*n1*n2] *
									(FD_D2_O6_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
											+ FD_D2_O6_N2(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
											+ FD_D2_O6_N3(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)) ;
						}

					}
				}
			}
			else if (fdOrder == 8)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma omp simd
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
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
			else if (fdOrder == 10)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma omp simd
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
									coef[i1+i2*n1+i3*n1*n2] *
									(FD_D2_O10_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
											+ FD_D2_O10_N2(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
											+ FD_D2_O10_N3(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)) ;
						}

					}
				}
			}
			else if (fdOrder == 12)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma omp simd
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
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
			else if (fdOrder == 14)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma omp simd
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
									coef[i1+i2*n1+i3*n1*n2] *
									(FD_D2_O14_N1(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
											+ FD_D2_O14_N2(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
											+ FD_D2_O14_N3(prc, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)) ;
						}

					}
				}
			}
			else if (fdOrder == 16)
			{
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
#pragma omp simd
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
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

	printDebug(FULL_DEBUG, "Out Grid_NEC::computePressureWithFD") ;
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_NEC::write(string fileName)
{
	printDebug(LIGHT_DEBUG, "IN Grid_NEC::write");

	// each proc write its own file

	if (Config::Instance()->writeGrid)
	{
		Myint nb_wrote;
		double t0 = MPI_Wtime() ;

		// write grid in binary format
		string fileName_bin = fileName + ".proc" + to_string(myMpiRank) + ".grid.bin" ;
		printInfo(MASTER, "* NEC::write on disk\t", fileName_bin.c_str());
		printDebug(LIGHT_DEBUG, "* write on disk\t", fileName_bin.c_str());
		
		int out_file;
		out_file = open(fileName_bin.c_str(), O_WRONLY | O_APPEND | O_CREAT, 0644);
		if(out_file == -1){
			printError("IN Grid_NEC::write, error while opening file") ;
			return(RTN_CODE_KO) ;
		}

		nb_wrote=::write(out_file, grid_3d, npoint * sizeof(Myfloat));
		if(nb_wrote < 0){
			printError("IN Grid_NEC::write, error while writing") ;
			return(RTN_CODE_KO) ;

		}

		close(out_file);

		// write grid info
		string fileName_info = fileName + ".proc" + to_string(myMpiRank) + ".grid.info" ;
		printDebug(LIGHT_DEBUG, "* write on disk\t", fileName_info.c_str());
		ofstream out_file2(fileName_info) ;
		out_file2 << this->n1 << "\n" ;
		out_file2 << this->n2 << "\n" ;
		out_file2 << this->n3 << "\n" ;
		out_file2.close() ;

		double t1 = MPI_Wtime() ;
		grid_write_time = t1-t0;
		printInfo(MASTER, "* time to write (s)", t1-t0) ;
	}
	else
	{
		grid_write_time = 0.0;
	}
	printDebug(LIGHT_DEBUG, "OUT Grid_NEC::write");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_NEC::write(Point_type pointType, string fileName)
{
	printDebug(LIGHT_DEBUG, "IN Grid_NEC::write");

	// each proc write its own file

	if (Config::Instance()->writeGrid)
	{		
		// particular case ALL_POINT
		if (pointType == ALL_POINTS)
		{
			this->write(fileName);
		}
		else
		{
			double t0 = MPI_Wtime() ;

			// write grid in binary format
			string fileName_bin = fileName + ".proc" + to_string(myMpiRank) + ".grid.bin";
			printInfo(MASTER, "* write on disk NEC_C:write\t", fileName_bin.c_str());
			printDebug(LIGHT_DEBUG, "* write on disk\t", fileName_bin.c_str());
			
			int out_file;
			out_file = open(fileName_bin.c_str(), O_WRONLY | O_APPEND | O_CREAT, 0644);
			if(out_file == -1){
				printError("IN Grid_NEC::write, error while opening file for writing") ;
				return(RTN_CODE_KO) ;
			}

			// get index range of desired points
			Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End;
			getGridIndex(pointType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End);

			//#pragma omp parallel for collapse(2)
			Myint64 n1Write = i1End - i1Start + 1;
			Myint64 n2Write = i2End - i2Start + 1;
			Myint64 n3Write = i3End - i3Start + 1;
			ssize_t nb_wrote = 0;
			for (Myint64 i3 = i3Start; i3 <= i3End; i3++)
			{
				for (Myint64 i2 = i2Start; i2 <= i2End; i2++)
				{
					// read file and fill grid along n1
					Myint64 i1 = i1Start;
					Myint64 idxGrid = i1 + i2 * n1 + i3 * n2 * n1;
					Myint64 idxFile = 1;

					nb_wrote=::write(out_file, &grid_3d[idxGrid], n1Write * sizeof(Myfloat));
					if(nb_wrote < 0){
						printError("IN Grid_NEC::write, Error while writing") ;
						return(RTN_CODE_KO) ;
					}
				}
			}

			close(out_file);

			// write grid info
			string fileName_info = fileName + ".proc" + to_string(myMpiRank) + ".grid.info";
			printDebug(LIGHT_DEBUG, "* write on disk\t", fileName_info.c_str());
			ofstream out_file2(fileName_info);
			out_file2 << n1Write << "\n";
			out_file2 << n2Write << "\n";
			out_file2 << n3Write << "\n";
			out_file2.close();

			double t1 = MPI_Wtime();
			grid_write_time = t1 - t0;
			printInfo(MASTER, "* time to write (s)", grid_write_time);
		}
	}
	else
	{
		grid_write_time = 0.0;
	}

	printDebug(LIGHT_DEBUG, "OUT Grid_NEC::write");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_NEC::read_all_at_once_then_loop(Point_type pointType, string fileName)
{
	//read all inside file overthrust
	
	printDebug(LIGHT_DEBUG, "IN Grid_NEC::read_all_at_once");
	printInfo(MASTER, "* read on disk\t", fileName.c_str());
	printDebug(LIGHT_DEBUG, "* read on disk\t", fileName.c_str());
	
	// open file
	int out_file;
	out_file= open(fileName.c_str(),O_RDONLY);
	if(out_file < 0){
		printError("IN Grid_NEC::read_all_at_once_then_loop, Error while opening file") ;
		return(RTN_CODE_KO) ;
	}

	// get index range of desired points
	Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End;
	getGridIndex(pointType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End);

	Myint64 nRead = i1End - i1Start + 1 ;

	// read all VPModel in one array
	int read_nb=0;
	int file_end_byte=0;

	file_end_byte=lseek(out_file,0,SEEK_END);
	Myfloat * all_modelVP = new Myfloat[file_end_byte/sizeof(Myfloat)];
	read_nb = pread(out_file, all_modelVP, file_end_byte ,0);
	if(read_nb == -1){
		printError("IN Grid_NEC::read_all_at_once_then_loop, Reading file failed") ;
		return(RTN_CODE_KO) ;
	}
	
	//#pragma omp parallel for collapse(2)
	for (Myint64 i3 = i3Start; i3 <= i3End; i3++)
	{
		for (Myint64 i2 = i2Start; i2 <= i2End; i2++)
		{
			// file is already read, copy data into grid_3d
			Myint64 i2file = i2-i2Start;	
			Myint64 i3file = i3-i3Start;			
			Myint64 i1 = i1Start ;
			Myint64 idxGrid = i1 + i2 * n1 + i3 * n2 * n1;

			Myint64 idxModelVP = (i1OffsetGlobInner + 
			Config::Instance()->n1 * (i2OffsetGlobInner + i2file) + 
			Config::Instance()->n1 * Config::Instance()->n2 * (i3OffsetGlobInner + i3file));  
			
			for(i1 = 0; i1 < nRead; i1++)
				grid_3d[idxGrid+i1] = all_modelVP[idxModelVP+i1];
		}
	}
	delete[] all_modelVP;
	close(out_file);
	printDebug(LIGHT_DEBUG, "OUT Grid_NEC::read_all_at_once");

}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_NEC::read_line_by_line(Point_type pointType, string fileName)
{
	
	printDebug(LIGHT_DEBUG, "IN Grid_NEC::read_line_by_line");
	// open file in binary format
	printInfo(MASTER, "* read on disk\t", fileName.c_str());
	printDebug(LIGHT_DEBUG, "* read on disk\t", fileName.c_str());
	
	int out_file;
	out_file= open(fileName.c_str(),O_RDONLY);
	if(out_file < 0){
		printError("IN Grid_NEC::read_line_by_line, error while opening file") ;
		return(RTN_CODE_KO) ;
	}

	// get index range of desired points
	Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End;
	getGridIndex(pointType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End);

	//#pragma omp parallel for collapse(2)
	Myint64 nRead = i1End - i1Start + 1 ;

	#ifdef _FTRACE
	ftrace_region_begin("lecture complete");
	#endif
	for (Myint64 i3 = i3Start; i3 <= i3End; i3++)
	{
		for (Myint64 i2 = i2Start; i2 <= i2End; i2++)
		{
			// read file and fill grid along n1
			Myint64 i2file = i2-i2Start;
			Myint64 i3file = i3-i3Start;
			Myint64 i1 = i1Start ;
			Myint64 idxGrid = i1 + i2 * n1 + i3 * n2 * n1;
			Myint64 myOffset = sizeof(Myfloat)*
			(i1OffsetGlobInner + 
			Config::Instance()->n1 * (i2OffsetGlobInner + i2file) + 
			Config::Instance()->n1 * Config::Instance()->n2 * (i3OffsetGlobInner + i3file));
			
			#ifdef _FTRACE
			ftrace_region_begin("read function time");
			#endif
					
			int read_nb=0;
			read_nb = pread(out_file, &grid_3d[idxGrid], nRead * sizeof(Myfloat) ,myOffset);
			if(read_nb == -1){
				printError("IN Grid_NEC::read_line_by_line, Reading file failed") ;
				return(RTN_CODE_KO) ;
			}

			#ifdef _FTRACE
			ftrace_region_end("read function time");
			#endif
		}
	}	
	// close file
	close(out_file);

	#ifdef _FTRACE
	ftrace_region_end("lecture complete");
	#endif
	printDebug(LIGHT_DEBUG, "OUT Grid_NEC::read_line_by_line");
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_NEC::read_plane_by_plane(Point_type pointType, string fileName)
{
	// read plane by plane IN 3D ONLY!!
	printDebug(LIGHT_DEBUG, "IN Grid_NEC::read_plane_by_plane");

	double t0 = MPI_Wtime();
	i_o_pread_time = 0 ;

	if(dim == DIM3){

		printInfo(MASTER, "* read on disk\t", fileName.c_str());
		printDebug(LIGHT_DEBUG, "* read on disk\t", fileName.c_str());
		
		double t0_open_time = MPI_Wtime();
		int out_file;
		out_file = open(fileName.c_str(),O_RDONLY);
		if(out_file == -1){
			printError("IN Grid_NEC::read_plane_by_plane, error while opening file") ;
			return(RTN_CODE_KO) ;
		}

		double t1_open_time = MPI_Wtime();
		i_o_pread_time += t1_open_time - t0_open_time;

		// get index range of desired points
		Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End;
		getGridIndex(pointType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End); //innerpoint

		//#pragma omp parallel for collapse(2)
		Myint64 n1Read = i1End - i1Start + 1 ; //=n1InnerLocal
		Myint64 nRead = n1Inner * n2Inner ; //points in a plane

		Myfloat * plane1_modelVP = new Myfloat[nRead];
		for(int i=0; i<nRead; i++)
			plane1_modelVP[i]=0;

		Myint64 read_nb=0;

		for (Myint64 i3 = i3Start; i3 <= i3End; i3++)
		{
			Myint64 myOffset = sizeof(Myfloat) * (i3-i3Start + i3OffsetGlobInner) * nRead ;

			double t0_pread_time = MPI_Wtime();

			read_nb = pread(out_file, plane1_modelVP, nRead * sizeof(Myfloat), myOffset);
			if(read_nb < 0){
				printError("IN Grid_NEC::read_plane_by_plane, error while reading file") ;
				return(RTN_CODE_KO) ;
			}

			double t1_pread_time = MPI_Wtime();
			i_o_pread_time += t1_pread_time - t0_pread_time ;

			for (Myint64 i2 = i2Start; i2 <= i2End; i2++)
			{
				// file is already read, copy data into grid_3d
				Myint i2file = i2-i2Start;
				Myint i1 = i1Start ;
				Myint i1File = i1 - i1Start;
				Myint64 idxGrid = i1 + i2 * n1 + i3 * n2 * n1;

				Myint64 idxModelVP = (i1OffsetGlobInner + i1File +
				n1Inner * (i2OffsetGlobInner + i2file));
				
				for(Myint64 i1plane = 0; i1plane < n1Read; i1plane++)
					grid_3d[idxGrid + i1plane] = plane1_modelVP[idxModelVP + i1plane];
			}
		}
		// close file
		delete[] plane1_modelVP;
		double t0_close_time = MPI_Wtime() ;
		close(out_file);
		double t1_close_time = MPI_Wtime() ;
		i_o_pread_time += t1_close_time - t0_close_time;
	}
	else{
		printInfo(MASTER, "Cannot read plane by plane because grid dim != 3 !!\n Reading line by line instead");
		double t0_read_line_time = MPI_Wtime();
		read_line_by_line(pointType, fileName);
		double t1_read_line_time = MPI_Wtime();
		i_o_pread_time += t1_read_line_time - t0_read_line_time;

	}
	printDebug(LIGHT_DEBUG, "OUT Grid_NEC::read_plane_by_plane");
	return(RTN_CODE_OK);

}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_NEC::read(string fileName)
{
	printDebug(LIGHT_DEBUG, "IN Grid_NEC::read");
	double t0 = MPI_Wtime() ;
	int nb_read=0;

	// read grid in binary format
	printInfo(MASTER, "* read on disk\t", fileName.c_str());
	printDebug(LIGHT_DEBUG, "* read on disk\t", fileName.c_str());
	
	int out_file;
	out_file = open(fileName.c_str(), O_RDONLY);
	if(out_file < 0){
		printError("IN Grid_NEC::read, error while Error while opening file") ;
		return(RTN_CODE_KO) ;
	}
	nb_read = ::read(out_file,grid_3d,npoint * sizeof(Myfloat));
	if(nb_read <0){
		printError("IN Grid_NEC::read, error while Error while reading") ;
		return(RTN_CODE_KO) ;
	}

	close(out_file);

	double t1 = MPI_Wtime();
	double readTime = t1 - t0 ;
	printInfo(MASTER, "* read time (s)\t", readTime);
	printInfo(MASTER, "* read bwth (GB/s)", npoint * sizeof(Myfloat) / readTime / 1e9);

	printDebug(LIGHT_DEBUG, "OUT Grid_NEC::read");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_NEC::read(Point_type pointType, string fileName)
{
	#ifdef _FTRACE
	ftrace_region_begin("all read");
	#endif
	printDebug(LIGHT_DEBUG, "IN Grid_NEC::read");
	double t0 = MPI_Wtime() ;

	// particular case ALL_POINT
	if (pointType == ALL_POINTS)
	{
		this->read(fileName);
	}
	else
	{
		read_plane_by_plane(pointType,fileName);
	}	

	double t1 = MPI_Wtime();
	double readTime = t1 - t0 ;
	printInfo(MASTER, "* read time (s)\t", readTime);
	printInfo(MASTER, "* read bwth (GB/s)", npoint * sizeof(Myfloat) / readTime / 1e9);

	printDebug(LIGHT_DEBUG, "OUT Grid_NEC::read");
	#ifdef _FTRACE
	ftrace_region_end("all read");
	#endif
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_NEC::exchangeHalos(MPI_comm_mode_type commMode)
{
	Myint64 li1s, li1e, li2s, li2e, li3s, li3e, li1n, li2n, li3n, n;
	Myfloat *bufSend, *bufRecv;
	MPI_Request requests[12];
	MPI_Status statuses[12];
	Myint64 req_li1s[12], req_li1e[12], req_li2s[12], req_li2e[12], req_li3s[12], req_li3e[12], req_li1n[12], req_li2n[12];
	Myfloat *req_bufRecv[12];
	Myint nreq = 0;
	Myint64 ofs = 0;

	printDebug(FULL_DEBUG, "IN Grid_NEC::exchangeHalos");

	if (commMode == MPI_COMM_MODE_SENDRECV)
	{
		printDebug(MID_DEBUG, "MPI_COMM_MODE_SENDRECV") ;

		// exchange the 6 halos
#if 0
		if (exchangeHalo(MPI_COMM_MODE_SENDRECV, I1HALO1) != RTN_CODE_OK) return (RTN_CODE_KO) ;
		if (exchangeHalo(MPI_COMM_MODE_SENDRECV, I1HALO2) != RTN_CODE_OK) return (RTN_CODE_KO) ;
		if (exchangeHalo(MPI_COMM_MODE_SENDRECV, I2HALO1) != RTN_CODE_OK) return (RTN_CODE_KO) ;
		if (exchangeHalo(MPI_COMM_MODE_SENDRECV, I2HALO2) != RTN_CODE_OK) return (RTN_CODE_KO) ;
		if (exchangeHalo(MPI_COMM_MODE_SENDRECV, I3HALO1) != RTN_CODE_OK) return (RTN_CODE_KO) ;
		if (exchangeHalo(MPI_COMM_MODE_SENDRECV, I3HALO2) != RTN_CODE_OK) return (RTN_CODE_KO) ;
#else
	// I1HALO1 - send/recv to i1ProcIdStart
	if (i1ProcIdStart != MPI_PROC_NULL) {
		li1s = i1InnerStart;
		li1e = i1InnerStart + haloWidth - 1;
		li2s = i2InnerStart;
		li2e = i2InnerEnd;
		li3s = i3InnerStart;
		li3e = i3InnerEnd;
		li1n = li1e - li1s + 1;
		li2n = li2e - li2s + 1;
		li3n = li3e - li3s + 1;
		n = li1n * li2n * li3n;
		bufSend = &(tmp_grid_3d[ofs]); ofs += n;
		bufRecv = &(tmp_grid_3d[ofs]); ofs += n;
		for (Myint64 i3 = li3s; i3<= li3e; i3++) {
#pragma _NEC interchange
			for (Myint64 i2 = li2s; i2<= li2e; i2++) {
				for (Myint64 i1 = li1s; i1<= li1e; i1++) {
					bufSend[(i1-li1s) + (i2-li2s) * li1n + (i3-li3s) * li1n * li2n] = grid_3d[i1 + i2*n1 + i3*n1*n2];
				}
			}
		}
		MPI_Irecv(bufRecv, n, MPI_MYFLOAT, i1ProcIdStart, 0, MPI_COMM_WORLD, &requests[nreq]);
		MPI_Isend(bufSend, n, MPI_MYFLOAT, i1ProcIdStart, 0, MPI_COMM_WORLD, &requests[nreq+1]);
		req_li1s[nreq] = i1Halo1Start;
		req_li1e[nreq] = i1Halo1End;
		req_li2s[nreq] = li2s;
		req_li2e[nreq] = li2e;
		req_li3s[nreq] = li3s;
		req_li3e[nreq] = li3e;
		req_li1n[nreq] = li1n;
		req_li2n[nreq] = li2n;
		req_bufRecv[nreq] = bufRecv;
		req_bufRecv[nreq+1] = NULL;
		nreq += 2;
	}

	// I1HALO2 - send/recv to i1ProcIdEnd
	if (i1ProcIdEnd != MPI_PROC_NULL) {
		li1s = i1InnerEnd - haloWidth + 1;
		li1e = i1InnerEnd;
		li2s = i2InnerStart;
		li2e = i2InnerEnd;
		li3s = i3InnerStart;
		li3e = i3InnerEnd;
		li1n = li1e - li1s + 1;
		li2n = li2e - li2s + 1;
		li3n = li3e - li3s + 1;
		n = li1n * li2n * li3n;
		bufSend = &(tmp_grid_3d[ofs]); ofs += n;
		bufRecv = &(tmp_grid_3d[ofs]); ofs += n;
		for (Myint64 i3 = li3s; i3<= li3e; i3++) {
#pragma _NEC interchange
			for (Myint64 i2 = li2s; i2<= li2e; i2++) {
				for (Myint64 i1 = li1s; i1<= li1e; i1++) {
					bufSend[(i1-li1s) + (i2-li2s) * li1n + (i3-li3s) * li1n * li2n] = grid_3d[i1 + i2*n1 + i3*n1*n2];
				}
			}
		}
		MPI_Irecv(bufRecv, n, MPI_MYFLOAT, i1ProcIdEnd, 0, MPI_COMM_WORLD, &requests[nreq]);
		MPI_Isend(bufSend, n, MPI_MYFLOAT, i1ProcIdEnd, 0, MPI_COMM_WORLD, &requests[nreq+1]);
		req_li1s[nreq] = i1Halo2Start;
		req_li1e[nreq] = i1Halo2End;
		req_li2s[nreq] = li2s;
		req_li2e[nreq] = li2e;
		req_li3s[nreq] = li3s;
		req_li3e[nreq] = li3e;
		req_li1n[nreq] = li1n;
		req_li2n[nreq] = li2n;
		req_bufRecv[nreq] = bufRecv;
		req_bufRecv[nreq+1] = NULL;
		nreq += 2;
	}

	// I2HALO1 - send/recv to i2ProcIdStart
	if (i2ProcIdStart != MPI_PROC_NULL) {
		li1s = i1InnerStart;
		li1e = i1InnerEnd;
		li2s = i2InnerStart;
		li2e = i2InnerStart + haloWidth - 1;
		li3s = i3InnerStart;
		li3e = i3InnerEnd;
		li1n = li1e - li1s + 1;
		li2n = li2e - li2s + 1;
		li3n = li3e - li3s + 1;
		n = li1n * li2n * li3n;
		bufSend = &(tmp_grid_3d[ofs]); ofs += n;
		bufRecv = &(tmp_grid_3d[ofs]); ofs += n;
		for (Myint64 i3 = li3s; i3<= li3e; i3++) {
			for (Myint64 i2 = li2s; i2<= li2e; i2++) {
				for (Myint64 i1 = li1s; i1<= li1e; i1++) {
					bufSend[(i1-li1s) + (i2-li2s) * li1n + (i3-li3s) * li1n * li2n] = grid_3d[i1 + i2*n1 + i3*n1*n2];
				}
			}
		}
		MPI_Irecv(bufRecv, n, MPI_MYFLOAT, i2ProcIdStart, 0, MPI_COMM_WORLD, &requests[nreq]);
		MPI_Isend(bufSend, n, MPI_MYFLOAT, i2ProcIdStart, 0, MPI_COMM_WORLD, &requests[nreq+1]);
		req_li1s[nreq] = li1s;
		req_li1e[nreq] = li1e;
		req_li2s[nreq] = i2Halo1Start;
		req_li2e[nreq] = i2Halo1End;
		req_li3s[nreq] = li3s;
		req_li3e[nreq] = li3e;
		req_li1n[nreq] = li1n;
		req_li2n[nreq] = li2n;
		req_bufRecv[nreq] = bufRecv;
		req_bufRecv[nreq+1] = NULL;
		nreq += 2;
	}

	// I2HALO2 - send/recv to i2ProcIdEnd
	if (i2ProcIdEnd != MPI_PROC_NULL) {
		li1s = i1InnerStart;
		li1e = i1InnerEnd;
		li2s = i2InnerEnd - haloWidth + 1;
		li2e = i2InnerEnd;
		li3s = i3InnerStart;
		li3e = i3InnerEnd;
		li1n = li1e - li1s + 1;
		li2n = li2e - li2s + 1;
		li3n = li3e - li3s + 1;
		n = li1n * li2n * li3n;
		bufSend = &(tmp_grid_3d[ofs]); ofs += n;
		bufRecv = &(tmp_grid_3d[ofs]); ofs += n;
		for (Myint64 i3 = li3s; i3<= li3e; i3++) {
			for (Myint64 i2 = li2s; i2<= li2e; i2++) {
				for (Myint64 i1 = li1s; i1<= li1e; i1++) {
					bufSend[(i1-li1s) + (i2-li2s) * li1n + (i3-li3s) * li1n * li2n] = grid_3d[i1 + i2*n1 + i3*n1*n2];
				}
			}
		}
		MPI_Irecv(bufRecv, n, MPI_MYFLOAT, i2ProcIdEnd, 0, MPI_COMM_WORLD, &requests[nreq]);
		MPI_Isend(bufSend, n, MPI_MYFLOAT, i2ProcIdEnd, 0, MPI_COMM_WORLD, &requests[nreq+1]);
		req_li1s[nreq] = li1s;
		req_li1e[nreq] = li1e;
		req_li2s[nreq] = i2Halo2Start;
		req_li2e[nreq] = i2Halo2End;
		req_li3s[nreq] = li3s;
		req_li3e[nreq] = li3e;
		req_li1n[nreq] = li1n;
		req_li2n[nreq] = li2n;
		req_bufRecv[nreq] = bufRecv;
		req_bufRecv[nreq+1] = NULL;
		nreq += 2;
	}

	// I3HALO1 - send/recv to i3ProcIdStart
	if (i3ProcIdStart != MPI_PROC_NULL) {
		li1s = i1InnerStart;
		li1e = i1InnerEnd;
		li2s = i2InnerStart;
		li2e = i2InnerEnd;
		li3s = i3InnerStart;
		li3e = i3InnerStart + haloWidth - 1;
		li1n = li1e - li1s + 1;
		li1n = li1e - li1s + 1;
		li2n = li2e - li2s + 1;
		li3n = li3e - li3s + 1;
		n = li1n * li2n * li3n;
		bufSend = &(tmp_grid_3d[ofs]); ofs += n;
		bufRecv = &(tmp_grid_3d[ofs]); ofs += n;
		for (Myint64 i3 = li3s; i3<= li3e; i3++) {
			for (Myint64 i2 = li2s; i2<= li2e; i2++) {
				for (Myint64 i1 = li1s; i1<= li1e; i1++) {
					bufSend[(i1-li1s) + (i2-li2s) * li1n + (i3-li3s) * li1n * li2n] = grid_3d[i1 + i2*n1 + i3*n1*n2];
				}
			}
		}
		MPI_Irecv(bufRecv, n, MPI_MYFLOAT, i3ProcIdStart, 0, MPI_COMM_WORLD, &requests[nreq]);
		MPI_Isend(bufSend, n, MPI_MYFLOAT, i3ProcIdStart, 0, MPI_COMM_WORLD, &requests[nreq+1]);
		req_li1s[nreq] = li1s;
		req_li1e[nreq] = li1e;
		req_li2s[nreq] = li2s;
		req_li2e[nreq] = li2e;
		req_li3s[nreq] = i3Halo1Start;
		req_li3e[nreq] = i3Halo1End;
		req_li1n[nreq] = li1n;
		req_li2n[nreq] = li2n;
		req_bufRecv[nreq] = bufRecv;
		req_bufRecv[nreq+1] = NULL;
		nreq += 2;
	}

	// I3HALO2 - send/recv to i3ProcIdEnd
	if (i3ProcIdEnd != MPI_PROC_NULL) {
		li1s = i1InnerStart;
		li1e = i1InnerEnd;
		li2s = i2InnerStart;
		li2e = i2InnerEnd;
		li3s = i3InnerEnd - haloWidth + 1;
		li3e = i3InnerEnd;
		li1n = li1e - li1s + 1;
		li2n = li2e - li2s + 1;
		li3n = li3e - li3s + 1;
		n = li1n * li2n * li3n;
		bufSend = &(tmp_grid_3d[ofs]); ofs += n;
		bufRecv = &(tmp_grid_3d[ofs]); ofs += n;
		for (Myint64 i3 = li3s; i3<= li3e; i3++) {
			for (Myint64 i2 = li2s; i2<= li2e; i2++) {
				for (Myint64 i1 = li1s; i1<= li1e; i1++) {
					bufSend[(i1-li1s) + (i2-li2s) * li1n + (i3-li3s) * li1n * li2n] = grid_3d[i1 + i2*n1 + i3*n1*n2];
				}
			}
		}
		MPI_Irecv(bufRecv, n, MPI_MYFLOAT, i3ProcIdEnd, 0, MPI_COMM_WORLD, &requests[nreq]);
		MPI_Isend(bufSend, n, MPI_MYFLOAT, i3ProcIdEnd, 0, MPI_COMM_WORLD, &requests[nreq+1]);
		req_li1s[nreq] = li1s;
		req_li1e[nreq] = li1e;
		req_li2s[nreq] = li2s;
		req_li2e[nreq] = li2e;
		req_li3s[nreq] = i3Halo2Start;
		req_li3e[nreq] = i3Halo2End;
		req_li1n[nreq] = li1n;
		req_li2n[nreq] = li2n;
		req_bufRecv[nreq] = bufRecv;
		req_bufRecv[nreq+1] = NULL;
		nreq += 2;
	}

	n = nreq;
	while (n > 0) {
		int indx;
		MPI_Waitany(nreq, requests, &indx, statuses);
		if (req_bufRecv[indx]) {
			li1s = req_li1s[indx];
			li1e = req_li1e[indx];
			li2s = req_li2s[indx];
			li2e = req_li2e[indx];
			li3s = req_li3s[indx];
			li3e = req_li3e[indx];
			li1n = req_li1n[indx];
			li2n = req_li2n[indx];
			bufRecv = req_bufRecv[indx];
			if (li1n > 32) {
				for (Myint64 i3 = li3s; i3<= li3e; i3++) {
					for (Myint64 i2 = li2s; i2<= li2e; i2++) {
						for (Myint64 i1 = li1s; i1<= li1e; i1++) {
							grid_3d[i1 + i2*n1 + i3*n1*n2] = bufRecv[(i1-li1s) + (i2-li2s) * li1n + (i3-li3s) * li1n * li2n];
						}
					}
				}
			} else {
				for (Myint64 i3 = li3s; i3<= li3e; i3++) {
#pragma _NEC interchange
					for (Myint64 i2 = li2s; i2<= li2e; i2++) {
						for (Myint64 i1 = li1s; i1<= li1e; i1++) {
							grid_3d[i1 + i2*n1 + i3*n1*n2] = bufRecv[(i1-li1s) + (i2-li2s) * li1n + (i3-li3s) * li1n * li2n];
						}
					}
				}
			}
		}
		n--;
	}
#endif
	}
	else
	{
		printError("IN Grid_NEC::exchangeHalos, invalid commMode", commMode) ;
		return(RTN_CODE_KO) ;
	}

	printDebug(FULL_DEBUG, "OUT Grid_NEC::exchangeHalos");
	return(RTN_CODE_OK) ;
}

} // namespace hpcscan
