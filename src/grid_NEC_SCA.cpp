
//-------------------------------------------------------------------------------------------------------
// Derived class from Grid
// Optimized for with NEC SCA
//-------------------------------------------------------------------------------------------------------

#include "grid_NEC_SCA.h"

#include <cassert>
#include <cfloat>  // for FLT_MAX ;
#include <cmath>   // for fabs
#include <cstddef> // for NULL
#include <fstream>
#include <stdio.h>

#include "mpi.h"
#ifdef __NEC__
#include <sca.h>
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

#ifdef _DOUBLE_PRECISION_
#define sca_stencil_create sca_stencil_create_d
#define sca_stencil_set_factor sca_stencil_set_factor_d
#define sca_stencil_set_input_array sca_stencil_set_input_array_d
#define sca_stencil_set_output_array sca_stencil_set_output_array_d
#else
#define sca_stencil_create sca_stencil_create_s
#define sca_stencil_set_factor sca_stencil_set_factor_s
#define sca_stencil_set_input_array sca_stencil_set_input_array_s
#define sca_stencil_set_output_array sca_stencil_set_output_array_s
#endif

using namespace std;

namespace hpcscan {

//-------------------------------------------------------------------------------------------------------

Grid_NEC_SCA::Grid_NEC_SCA(Grid_type gridTypeIn) : Grid(gridTypeIn)
{
	printDebug(MID_DEBUG, "IN Grid_NEC_SCA::Grid_NEC_SCA");

	gridMode = "NEC_SCA" ;

	flag_code_FD_D2_N1  = false ;
	flag_code_FD_D2_N2  = false ;
	flag_code_FD_D2_N3  = false ;
	flag_code_FD_LAPLACIAN = false ;

	printDebug(MID_DEBUG, "OUT Grid_NEC_SCA::Grid_NEC_SCA");
}

//-------------------------------------------------------------------------------------------------------

Grid_NEC_SCA::Grid_NEC_SCA(Grid_type gridTypeIn, Dim_type dimIn,
		Myint64 n1InnerIn, Myint64 n2InnerIn, Myint64 n3InnerIn) : Grid(gridTypeIn, dimIn,
				n1InnerIn, n2InnerIn, n3InnerIn)
{
	printDebug(MID_DEBUG, "IN Grid_NEC_SCA::Grid_NEC_SCA");

	gridMode = "CacheBlk" ;

	flag_code_FD_D2_N1  = false ;
	flag_code_FD_D2_N2  = false ;
	flag_code_FD_D2_N3  = false ;
	flag_code_FD_LAPLACIAN = false ;

	printDebug(MID_DEBUG, "OUT Grid_NEC_SCA::Grid_NEC_SCA");
}

//-------------------------------------------------------------------------------------------------------

Grid_NEC_SCA::~Grid_NEC_SCA(void)
{
	printDebug(MID_DEBUG, "IN Grid_NEC_SCA::~Grid_NEC_SCA");

	//====================================================
	// Destroy Stencil Code
	//====================================================
#ifdef __NEC__
	if (flag_code_FD_D2_N1) sca_code_destroy(code_FD_D2_N1);
	if (flag_code_FD_D2_N2) sca_code_destroy(code_FD_D2_N2);
	if (flag_code_FD_D2_N3) sca_code_destroy(code_FD_D2_N3);
	if (flag_code_FD_LAPLACIAN) sca_code_destroy(code_FD_LAPLACIAN);
#endif

	printDebug(MID_DEBUG, "OUT Grid_NEC_SCA::~Grid_NEC_SCA");
}

//-------------------------------------------------------------------------------------------------------

void Grid_NEC_SCA::info(void)
{
	printDebug(FULL_DEBUG, "IN Grid_NEC_SCA::info");

	// parent class info
	Grid::info() ;

	// additional info
	printInfo(MASTER, " TO BE COMPLETED") ;

	printDebug(FULL_DEBUG, "IN Grid_NEC_SCA::info");
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_NEC_SCA::FD_D2_N1(Point_type pType, const Grid& Wgrid, Myint fdOrder)
{
	printDebug(MID_DEBUG, "IN Grid_NEC_SCA::FD_D2_N1");

	// check grids are same size
	if (this->sameSize(Wgrid) != true)
	{
		printError("Grid_NEC_SCA::FD_D2_N1, grids have not same size") ;
		return(RTN_CODE_KO) ;
	}

	// initialize code_FD_D2_N1 at first call
	if (!flag_code_FD_D2_N1) initialize_code_FD_D2_N1(pType, Wgrid, fdOrder) ;

	// compute FD along N1
#ifdef __NEC__
	sca_code_execute(code_FD_D2_N1);
#endif

	printDebug(MID_DEBUG, "OUT Grid_NEC_SCA::FD_D2_N1");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_NEC_SCA::FD_D2_N2(Point_type pType, const Grid& Wgrid, Myint fdOrder)
{
	printDebug(MID_DEBUG, "IN Grid_NEC_SCA::FD_D2_N2");

	// check grids are same size
	if (this->sameSize(Wgrid) != true)
	{
		printError("Grid_NEC_SCA::FD_D2_N2, grids have not same size") ;
		return(RTN_CODE_KO) ;
	}

	// initialize code_FD_D2_N2 at first call
	if (!flag_code_FD_D2_N2) initialize_code_FD_D2_N2(pType, Wgrid, fdOrder) ;

	// compute FD along N2
#ifdef __NEC__
	sca_code_execute(code_FD_D2_N2);
#endif

	printDebug(MID_DEBUG, "OUT Grid_NEC_SCA::FD_D2_N2");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_NEC_SCA::FD_D2_N3(Point_type pType, const Grid& Wgrid, Myint fdOrder)
{
	printDebug(MID_DEBUG, "IN Grid_NEC_SCA::FD_D2_N3");

	// check grids are same size
	if (this->sameSize(Wgrid) != true)
	{
		printError("Grid_NEC_SCA::FD_D2_N3, grids have not same size") ;
		return(RTN_CODE_KO) ;
	}

	// initialize code_FD_D2_N3 at first call
	if (!flag_code_FD_D2_N3) initialize_code_FD_D2_N3(pType, Wgrid, fdOrder) ;

	// compute FD along N3
#ifdef __NEC__
	sca_code_execute(code_FD_D2_N3);
#endif

	printDebug(MID_DEBUG, "OUT Grid_NEC_SCA::FD_D2_N3");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_NEC_SCA::FD_LAPLACIAN(Point_type pType, const Grid& Wgrid, Myint fdOrder)
{
	printDebug(MID_DEBUG, "IN Grid_NEC_SCA::FD_LAPLACIAN");

	// check grids are same size
	if (this->sameSize(Wgrid) != true)
	{
		printError("Grid_NEC_SCA::FD_LAPLACIAN, grids have not same size") ;
		return(RTN_CODE_KO) ;
	}

	// initialize code_FD_LAPLACIAN at first call
	if (!flag_code_FD_LAPLACIAN) initialize_code_FD_LAPLACIAN(pType, Wgrid, fdOrder) ;

	// compute FD along N3
#ifdef __NEC__
	sca_code_execute(code_FD_LAPLACIAN);
#endif

	printDebug(MID_DEBUG, "OUT Grid_NEC_SCA::FD_LAPLACIAN");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Myfloat Grid_NEC_SCA::getMin(Point_type pointType)
{
	printDebug(LIGHT_DEBUG, "IN Grid_NEC_SCA::getMin");

	Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
	getGridIndex(pointType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End) ;
	Myfloat val = FLT_MAX ;

#pragma omp parallel for reduction(min:val)
	for (Myint64 i3 = i3Start; i3<= i3End; i3++)
	{
		for (Myint64 i2 = i1Start + i2Start * n1; i2<= i1End + i2End * n1; i2++)
		{
			Myint64 ii = i2 + i3*n2*n1 ;
			Myint64 i1 = (i2 - i1Start) % n1 + i1Start;
			if (i1 >= i1Start && i1 <= i1End && grid_3d[ii] < val) val = grid_3d[ii] ;
		}
	}

	printDebug(LIGHT_DEBUG, "Min val", val);

	printDebug(LIGHT_DEBUG, "OUT Grid_NEC_SCA::getMin");
	return(val) ;
}

//-------------------------------------------------------------------------------------------------------

Myfloat Grid_NEC_SCA::getMax(Point_type pointType)
{
	printDebug(LIGHT_DEBUG, "IN Grid_NEC_SCA::getMax");

	Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
	getGridIndex(pointType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End) ;

	Myfloat val = -FLT_MAX ;

#pragma omp parallel for reduction(max:val)
	for (Myint64 i3 = i3Start; i3<= i3End; i3++)
	{
		for (Myint64 i2 = i1Start + i2Start * n1; i2<= i1End + i2End * n1; i2++)
		{
			Myint64 ii = i2 + i3*n2*n1 ;
			Myint64 i1 = (i2 - i1Start) % n1 + i1Start;
			if (i1 >= i1Start && i1 <= i1End && grid_3d[ii] > val) val = grid_3d[ii] ;
		}
	}

	printDebug(LIGHT_DEBUG, "Max val", val);

	printDebug(LIGHT_DEBUG, "OUT Grid_NEC_SCA::getMax");
	return(val) ;
}

//-------------------------------------------------------------------------------------------------------

Myfloat Grid_NEC_SCA::maxErr(Point_type pointType, const Grid& gridIn) const
{
	printDebug(FULL_DEBUG, "IN Grid_NEC_SCA::maxErr");

	// check grids have same size
	if (!(this->sameSize(gridIn)))
	{
		printError("Grid_NEC_SCA::maxErr, grids have different size") ;
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

	printDebug(FULL_DEBUG, "OUT Grid_NEC_SCA::maxErr");
	return(err) ;
}

//-------------------------------------------------------------------------------------------------------

void Grid_NEC_SCA::padGridn1(void)
{
	printDebug(MID_DEBUG, "IN Grid_NEC_SCA::padGridn1");

	if (Config::Instance()->autoPad == true)
	{
#ifndef __NEC__
		printWarning("Grid_NEC_SCA::padGridn1, autoPad is not available on your platform") ;
		i1PadStart = i1Halo2End ;
		i1PadEnd   = i1Halo2End ;
#else
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
#endif
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

	printDebug(MID_DEBUG, "OUT Grid_NEC_SCA::padGridn1");
}

//-------------------------------------------------------------------------------------------------------

void Grid_NEC_SCA::padGridn2(void)
{
	printDebug(MID_DEBUG, "IN Grid_NEC_SCA::padGridn2");

	if (Config::Instance()->autoPad == true)
	{
#ifndef __NEC__
		printWarning("Grid_NEC_SCA::padGridn2, autoPad is not available on your platform") ;
		i2PadStart = i2Halo2End ;
		i2PadEnd   = i2Halo2End ;
#else
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
#endif
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

	printDebug(MID_DEBUG, "OUT Grid_NEC_SCA::padGridn2");
}

//-------------------------------------------------------------------------------------------------------

void Grid_NEC_SCA::padGridn3(void)
{
	printDebug(MID_DEBUG, "IN Grid_NEC_SCA::padGridn3");

	if (Config::Instance()->autoPad == true)
	{
#ifndef __NEC__
		printWarning("Grid_NEC_SCA::padGridn3, autoPad is not available on your platform") ;
		i3PadStart = i3Halo2End ;
		i3PadEnd   = i3Halo2End ;
#else
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
#endif
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

	printDebug(MID_DEBUG, "OUT Grid_NEC_SCA::padGridn3");
}

//-------------------------------------------------------------------------------------------------------

void Grid_NEC_SCA::initializeGrid(void)
{
	printDebug(MID_DEBUG, "IN Grid_NEC_SCA::initializeGrid");

	Grid::initializeGrid() ;

	printDebug(MID_DEBUG, "OUT Grid_NEC_SCA::initializeGrid");
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_NEC_SCA::initialize_code_FD_D2_N1(Point_type pType, const Grid& Wgrid, Myint fdOrder)
{
	printDebug(MID_DEBUG, "IN Grid_NEC_SCA::initialize_code_FD_D2_N1");

#ifdef __NEC__
	if (!flag_code_FD_D2_N1)
	{
		sca_stencil_t sten ;
		sca_stencil_create(&sten);
		Myint nPtPerStencil = getPtPerStencilFD_D2(fdOrder) ;
		sca_stencil_append_elements(sten, nPtPerStencil);

		const Myfloat inv_d1  = Myfloat(1.0) / d1 ;
		const Myfloat inv_d2  = Myfloat(1.0) / d2 ;
		const Myfloat inv_d3  = Myfloat(1.0) / d3 ;

		const Myfloat inv2_d1 = inv_d1 * inv_d1 ;
		const Myfloat inv2_d2 = inv_d2 * inv_d2 ;
		const Myfloat inv2_d3 = inv_d3 * inv_d3 ;

		// sca_stencil_set_location for all points in stencil
		if (fdOrder == 2)
		{
			sca_stencil_set_location(sten, 0, -1, 0, 0, 0); sca_stencil_set_factor(sten, 0, FD_D2_O2_A1 * inv2_d1); // i1-1, i2, i3
			sca_stencil_set_location(sten, 1,  0, 0, 0, 0); sca_stencil_set_factor(sten, 1, FD_D2_O2_A0 * inv2_d1); // i1  , i2, i3
			sca_stencil_set_location(sten, 2,  1, 0, 0, 0); sca_stencil_set_factor(sten, 2, FD_D2_O2_A1 * inv2_d1); // i1+1, i2, i3
		}
		else if (fdOrder == 4)
		{
			sca_stencil_set_location(sten, 0, -2, 0, 0, 0); sca_stencil_set_factor(sten, 0, FD_D2_O4_A2 * inv2_d1); // i1-2, i2, i3
			sca_stencil_set_location(sten, 1, -1, 0, 0, 0); sca_stencil_set_factor(sten, 1, FD_D2_O4_A1 * inv2_d1); // i1-1, i2, i3
			sca_stencil_set_location(sten, 2,  0, 0, 0, 0); sca_stencil_set_factor(sten, 2, FD_D2_O4_A0 * inv2_d1); // i1  , i2, i3
			sca_stencil_set_location(sten, 3,  1, 0, 0, 0); sca_stencil_set_factor(sten, 3, FD_D2_O4_A1 * inv2_d1); // i1+1, i2, i3
			sca_stencil_set_location(sten, 4,  2, 0, 0, 0); sca_stencil_set_factor(sten, 4, FD_D2_O4_A2 * inv2_d1); // i1+2, i2, i3
		}
		else if (fdOrder == 8)
		{
			sca_stencil_set_location(sten, 0, -4, 0, 0, 0); sca_stencil_set_factor(sten, 0, FD_D2_O8_A4 * inv2_d1); // i1-4, i2, i3
			sca_stencil_set_location(sten, 1, -3, 0, 0, 0); sca_stencil_set_factor(sten, 1, FD_D2_O8_A3 * inv2_d1); // i1-3, i2, i3
			sca_stencil_set_location(sten, 2, -2, 0, 0, 0); sca_stencil_set_factor(sten, 2, FD_D2_O8_A2 * inv2_d1); // i1-2, i2, i3
			sca_stencil_set_location(sten, 3, -1, 0, 0, 0); sca_stencil_set_factor(sten, 3, FD_D2_O8_A1 * inv2_d1); // i1-1, i2, i3
			sca_stencil_set_location(sten, 4,  0, 0, 0, 0); sca_stencil_set_factor(sten, 4, FD_D2_O8_A0 * inv2_d1); // i1  , i2, i3
			sca_stencil_set_location(sten, 5,  1, 0, 0, 0); sca_stencil_set_factor(sten, 5, FD_D2_O8_A1 * inv2_d1); // i1+1, i2, i3
			sca_stencil_set_location(sten, 6,  2, 0, 0, 0); sca_stencil_set_factor(sten, 6, FD_D2_O8_A2 * inv2_d1); // i1+2, i2, i3
			sca_stencil_set_location(sten, 7,  3, 0, 0, 0); sca_stencil_set_factor(sten, 7, FD_D2_O8_A3 * inv2_d1); // i1+3, i2, i3
			sca_stencil_set_location(sten, 8,  4, 0, 0, 0); sca_stencil_set_factor(sten, 8, FD_D2_O8_A4 * inv2_d1); // i1+4, i2, i3
		}
		else if (fdOrder == 12)
		{
			sca_stencil_set_location(sten,  0, -6, 0, 0, 0); sca_stencil_set_factor(sten,  0, FD_D2_O12_A6 * inv2_d1); // i1-6, i2, i3
			sca_stencil_set_location(sten,  1, -5, 0, 0, 0); sca_stencil_set_factor(sten,  1, FD_D2_O12_A5 * inv2_d1); // i1-5, i2, i3
			sca_stencil_set_location(sten,  2, -4, 0, 0, 0); sca_stencil_set_factor(sten,  2, FD_D2_O12_A4 * inv2_d1); // i1-4, i2, i3
			sca_stencil_set_location(sten,  3, -3, 0, 0, 0); sca_stencil_set_factor(sten,  3, FD_D2_O12_A3 * inv2_d1); // i1-3, i2, i3
			sca_stencil_set_location(sten,  4, -2, 0, 0, 0); sca_stencil_set_factor(sten,  4, FD_D2_O12_A2 * inv2_d1); // i1-2, i2, i3
			sca_stencil_set_location(sten,  5, -1, 0, 0, 0); sca_stencil_set_factor(sten,  5, FD_D2_O12_A1 * inv2_d1); // i1-1, i2, i3
			sca_stencil_set_location(sten,  6,  0, 0, 0, 0); sca_stencil_set_factor(sten,  6, FD_D2_O12_A0 * inv2_d1); // i1  , i2, i3
			sca_stencil_set_location(sten,  7,  1, 0, 0, 0); sca_stencil_set_factor(sten,  7, FD_D2_O12_A1 * inv2_d1); // i1+1, i2, i3
			sca_stencil_set_location(sten,  8,  2, 0, 0, 0); sca_stencil_set_factor(sten,  8, FD_D2_O12_A2 * inv2_d1); // i1+2, i2, i3
			sca_stencil_set_location(sten,  9,  3, 0, 0, 0); sca_stencil_set_factor(sten,  9, FD_D2_O12_A3 * inv2_d1); // i1+3, i2, i3
			sca_stencil_set_location(sten, 10,  4, 0, 0, 0); sca_stencil_set_factor(sten, 10, FD_D2_O12_A4 * inv2_d1); // i1+4, i2, i3
			sca_stencil_set_location(sten, 11,  5, 0, 0, 0); sca_stencil_set_factor(sten, 11, FD_D2_O12_A5 * inv2_d1); // i1+5, i2, i3
			sca_stencil_set_location(sten, 12,  6, 0, 0, 0); sca_stencil_set_factor(sten, 12, FD_D2_O12_A6 * inv2_d1); // i1+6, i2, i3
		}
		else if (fdOrder == 16)
		{
			sca_stencil_set_location(sten,  0, -8, 0, 0, 0); sca_stencil_set_factor(sten,  0, FD_D2_O16_A8 * inv2_d1); // i1-8, i2, i3
			sca_stencil_set_location(sten,  1, -7, 0, 0, 0); sca_stencil_set_factor(sten,  1, FD_D2_O16_A7 * inv2_d1); // i1-7, i2, i3
			sca_stencil_set_location(sten,  2, -6, 0, 0, 0); sca_stencil_set_factor(sten,  2, FD_D2_O16_A6 * inv2_d1); // i1-6, i2, i3
			sca_stencil_set_location(sten,  3, -5, 0, 0, 0); sca_stencil_set_factor(sten,  3, FD_D2_O16_A5 * inv2_d1); // i1-5, i2, i3
			sca_stencil_set_location(sten,  4, -4, 0, 0, 0); sca_stencil_set_factor(sten,  4, FD_D2_O16_A4 * inv2_d1); // i1-4, i2, i3
			sca_stencil_set_location(sten,  5, -3, 0, 0, 0); sca_stencil_set_factor(sten,  5, FD_D2_O16_A3 * inv2_d1); // i1-3, i2, i3
			sca_stencil_set_location(sten,  6, -2, 0, 0, 0); sca_stencil_set_factor(sten,  6, FD_D2_O16_A2 * inv2_d1); // i1-2, i2, i3
			sca_stencil_set_location(sten,  7, -1, 0, 0, 0); sca_stencil_set_factor(sten,  7, FD_D2_O16_A1 * inv2_d1); // i1-1, i2, i3
			sca_stencil_set_location(sten,  8,  0, 0, 0, 0); sca_stencil_set_factor(sten,  8, FD_D2_O16_A0 * inv2_d1); // i1  , i2, i3
			sca_stencil_set_location(sten,  9,  1, 0, 0, 0); sca_stencil_set_factor(sten,  9, FD_D2_O16_A1 * inv2_d1); // i1+1, i2, i3
			sca_stencil_set_location(sten, 10,  2, 0, 0, 0); sca_stencil_set_factor(sten, 10, FD_D2_O16_A2 * inv2_d1); // i1+2, i2, i3
			sca_stencil_set_location(sten, 11,  3, 0, 0, 0); sca_stencil_set_factor(sten, 11, FD_D2_O16_A3 * inv2_d1); // i1+3, i2, i3
			sca_stencil_set_location(sten, 12,  4, 0, 0, 0); sca_stencil_set_factor(sten, 12, FD_D2_O16_A4 * inv2_d1); // i1+4, i2, i3
			sca_stencil_set_location(sten, 13,  5, 0, 0, 0); sca_stencil_set_factor(sten, 13, FD_D2_O16_A5 * inv2_d1); // i1+5, i2, i3
			sca_stencil_set_location(sten, 14,  6, 0, 0, 0); sca_stencil_set_factor(sten, 14, FD_D2_O16_A6 * inv2_d1); // i1+6, i2, i3
			sca_stencil_set_location(sten, 15,  7, 0, 0, 0); sca_stencil_set_factor(sten, 15, FD_D2_O16_A7 * inv2_d1); // i1+7, i2, i3
			sca_stencil_set_location(sten, 16,  8, 0, 0, 0); sca_stencil_set_factor(sten, 16, FD_D2_O16_A8 * inv2_d1); // i1+8, i2, i3
		}

		Myfloat * const w = Wgrid.grid_3d ;
		Myfloat * const u = this->grid_3d ;

		// sca_code_create
		Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
		getGridIndex(pType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End) ;
		for (Myint i = 0; i < nPtPerStencil; i++) sca_stencil_set_input_array(sten, i, 1, n1, n2, n3, &u[i1Start + i2Start*n1 + i3Start*n2*n1]);
		sca_stencil_set_output_array(sten, 1, n1, n2, n3, &w[i1Start + i2Start*n1 + i3Start*n2*n1]);
		sca_code_create(&code_FD_D2_N1, sten, i1End-i1Start+1, i2End-i2Start+1, i3End-i3Start+1, 1);
		sca_stencil_destroy(sten);
	}

	flag_code_FD_D2_N1 = true ;
	return(RTN_CODE_OK) ;

#else
	printError("Grid_NEC_SCA::initialize_code_FD_D2_N1, not supported on your platform") ;
	return(RTN_CODE_KO) ;
#endif

	printDebug(MID_DEBUG, "OUT Grid_NEC_SCA::initialize_code_FD_D2_N1");
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_NEC_SCA::initialize_code_FD_D2_N2(Point_type pType, const Grid& Wgrid, Myint fdOrder)
{
	printDebug(MID_DEBUG, "IN Grid_NEC_SCA::initialize_code_FD_D2_N2");

#ifdef __NEC__
	if (!flag_code_FD_D2_N2)
	{
		sca_stencil_t sten ;
		sca_stencil_create(&sten);
		Myint nPtPerStencil = getPtPerStencilFD_D2(fdOrder) ;
		sca_stencil_append_elements(sten, nPtPerStencil);

		const Myfloat inv_d1  = Myfloat(1.0) / d1 ;
		const Myfloat inv_d2  = Myfloat(1.0) / d2 ;
		const Myfloat inv_d3  = Myfloat(1.0) / d3 ;

		const Myfloat inv2_d1 = inv_d1 * inv_d1 ;
		const Myfloat inv2_d2 = inv_d2 * inv_d2 ;
		const Myfloat inv2_d3 = inv_d3 * inv_d3 ;

		// sca_stencil_set_location for all points in stencil
		if (fdOrder == 2)
		{
			sca_stencil_set_location(sten, 0, 0, -1, 0, 0); sca_stencil_set_factor(sten, 0, FD_D2_O2_A1 * inv2_d2); // i1, i2-1, i3
			sca_stencil_set_location(sten, 1, 0,  0, 0, 0); sca_stencil_set_factor(sten, 1, FD_D2_O2_A0 * inv2_d2); // i1, i2  , i3
			sca_stencil_set_location(sten, 2, 0,  1, 0, 0); sca_stencil_set_factor(sten, 2, FD_D2_O2_A1 * inv2_d2); // i1, i2+1, i3
		}
		else if (fdOrder == 4)
		{
			sca_stencil_set_location(sten, 0, 0, -2, 0, 0); sca_stencil_set_factor(sten, 0, FD_D2_O4_A2 * inv2_d2); // i1, i2-2, i3
			sca_stencil_set_location(sten, 1, 0, -1, 0, 0); sca_stencil_set_factor(sten, 1, FD_D2_O4_A1 * inv2_d2); // i1, i2-1, i3
			sca_stencil_set_location(sten, 2, 0,  0, 0, 0); sca_stencil_set_factor(sten, 2, FD_D2_O4_A0 * inv2_d2); // i1, i2  , i3
			sca_stencil_set_location(sten, 3, 0,  1, 0, 0); sca_stencil_set_factor(sten, 3, FD_D2_O4_A1 * inv2_d2); // i1, i2+1, i3
			sca_stencil_set_location(sten, 4, 0,  2, 0, 0); sca_stencil_set_factor(sten, 4, FD_D2_O4_A2 * inv2_d2); // i1, i2+2, i3
		}
		else if (fdOrder == 8)
		{
			sca_stencil_set_location(sten, 0, 0, -4, 0, 0); sca_stencil_set_factor(sten, 0, FD_D2_O8_A4 * inv2_d2); // i1, i2-4, i3
			sca_stencil_set_location(sten, 1, 0, -3, 0, 0); sca_stencil_set_factor(sten, 1, FD_D2_O8_A3 * inv2_d2); // i1, i2-3, i3
			sca_stencil_set_location(sten, 2, 0, -2, 0, 0); sca_stencil_set_factor(sten, 2, FD_D2_O8_A2 * inv2_d2); // i1, i2-2, i3
			sca_stencil_set_location(sten, 3, 0, -1, 0, 0); sca_stencil_set_factor(sten, 3, FD_D2_O8_A1 * inv2_d2); // i1, i2-1, i3
			sca_stencil_set_location(sten, 4, 0,  0, 0, 0); sca_stencil_set_factor(sten, 4, FD_D2_O8_A0 * inv2_d2); // i1, i2  , i3
			sca_stencil_set_location(sten, 5, 0,  1, 0, 0); sca_stencil_set_factor(sten, 5, FD_D2_O8_A1 * inv2_d2); // i1, i2+1, i3
			sca_stencil_set_location(sten, 6, 0,  2, 0, 0); sca_stencil_set_factor(sten, 6, FD_D2_O8_A2 * inv2_d2); // i1, i2+2, i3
			sca_stencil_set_location(sten, 7, 0,  3, 0, 0); sca_stencil_set_factor(sten, 7, FD_D2_O8_A3 * inv2_d2); // i1, i2+3, i3
			sca_stencil_set_location(sten, 8, 0,  4, 0, 0); sca_stencil_set_factor(sten, 8, FD_D2_O8_A4 * inv2_d2); // i1, i2+4, i3
		}
		else if (fdOrder == 12)
		{
			sca_stencil_set_location(sten,  0, 0, -6, 0, 0); sca_stencil_set_factor(sten,  0, FD_D2_O12_A6 * inv2_d2); // i1, i2-6, i3
			sca_stencil_set_location(sten,  1, 0, -5, 0, 0); sca_stencil_set_factor(sten,  1, FD_D2_O12_A5 * inv2_d2); // i1, i2-5, i3
			sca_stencil_set_location(sten,  2, 0, -4, 0, 0); sca_stencil_set_factor(sten,  2, FD_D2_O12_A4 * inv2_d2); // i1, i2-4, i3
			sca_stencil_set_location(sten,  3, 0, -3, 0, 0); sca_stencil_set_factor(sten,  3, FD_D2_O12_A3 * inv2_d2); // i1, i2-3, i3
			sca_stencil_set_location(sten,  4, 0, -2, 0, 0); sca_stencil_set_factor(sten,  4, FD_D2_O12_A2 * inv2_d2); // i1, i2-2, i3
			sca_stencil_set_location(sten,  5, 0, -1, 0, 0); sca_stencil_set_factor(sten,  5, FD_D2_O12_A1 * inv2_d2); // i1, i2-1, i3
			sca_stencil_set_location(sten,  6, 0,  0, 0, 0); sca_stencil_set_factor(sten,  6, FD_D2_O12_A0 * inv2_d2); // i1, i2  , i3
			sca_stencil_set_location(sten,  7, 0,  1, 0, 0); sca_stencil_set_factor(sten,  7, FD_D2_O12_A1 * inv2_d2); // i1, i2+1, i3
			sca_stencil_set_location(sten,  8, 0,  2, 0, 0); sca_stencil_set_factor(sten,  8, FD_D2_O12_A2 * inv2_d2); // i1, i2+2, i3
			sca_stencil_set_location(sten,  9, 0,  3, 0, 0); sca_stencil_set_factor(sten,  9, FD_D2_O12_A3 * inv2_d2); // i1, i2+3, i3
			sca_stencil_set_location(sten, 10, 0,  4, 0, 0); sca_stencil_set_factor(sten, 10, FD_D2_O12_A4 * inv2_d2); // i1, i2+4, i3
			sca_stencil_set_location(sten, 11, 0,  5, 0, 0); sca_stencil_set_factor(sten, 11, FD_D2_O12_A5 * inv2_d2); // i1, i2+5, i3
			sca_stencil_set_location(sten, 12, 0,  6, 0, 0); sca_stencil_set_factor(sten, 12, FD_D2_O12_A6 * inv2_d2); // i1, i2+6, i3
		}
		else if (fdOrder == 16)
		{
			sca_stencil_set_location(sten,  0, 0, -8, 0, 0); sca_stencil_set_factor(sten,  0, FD_D2_O16_A8 * inv2_d2); // i1, i2-8, i3
			sca_stencil_set_location(sten,  1, 0, -7, 0, 0); sca_stencil_set_factor(sten,  1, FD_D2_O16_A7 * inv2_d2); // i1, i2-7, i3
			sca_stencil_set_location(sten,  2, 0, -6, 0, 0); sca_stencil_set_factor(sten,  2, FD_D2_O16_A6 * inv2_d2); // i1, i2-6, i3
			sca_stencil_set_location(sten,  3, 0, -5, 0, 0); sca_stencil_set_factor(sten,  3, FD_D2_O16_A5 * inv2_d2); // i1, i2-5, i3
			sca_stencil_set_location(sten,  4, 0, -4, 0, 0); sca_stencil_set_factor(sten,  4, FD_D2_O16_A4 * inv2_d2); // i1, i2-4, i3
			sca_stencil_set_location(sten,  5, 0, -3, 0, 0); sca_stencil_set_factor(sten,  5, FD_D2_O16_A3 * inv2_d2); // i1, i2-3, i3
			sca_stencil_set_location(sten,  6, 0, -2, 0, 0); sca_stencil_set_factor(sten,  6, FD_D2_O16_A2 * inv2_d2); // i1, i2-2, i3
			sca_stencil_set_location(sten,  7, 0, -1, 0, 0); sca_stencil_set_factor(sten,  7, FD_D2_O16_A1 * inv2_d2); // i1, i2-1, i3
			sca_stencil_set_location(sten,  8, 0,  0, 0, 0); sca_stencil_set_factor(sten,  8, FD_D2_O16_A0 * inv2_d2); // i1, i2  , i3
			sca_stencil_set_location(sten,  9, 0,  1, 0, 0); sca_stencil_set_factor(sten,  9, FD_D2_O16_A1 * inv2_d2); // i1, i2+1, i3
			sca_stencil_set_location(sten, 10, 0,  2, 0, 0); sca_stencil_set_factor(sten, 10, FD_D2_O16_A2 * inv2_d2); // i1, i2+2, i3
			sca_stencil_set_location(sten, 11, 0,  3, 0, 0); sca_stencil_set_factor(sten, 11, FD_D2_O16_A3 * inv2_d2); // i1, i2+3, i3
			sca_stencil_set_location(sten, 12, 0,  4, 0, 0); sca_stencil_set_factor(sten, 12, FD_D2_O16_A4 * inv2_d2); // i1, i2+4, i3
			sca_stencil_set_location(sten, 13, 0,  5, 0, 0); sca_stencil_set_factor(sten, 13, FD_D2_O16_A5 * inv2_d2); // i1, i2+5, i3
			sca_stencil_set_location(sten, 14, 0,  6, 0, 0); sca_stencil_set_factor(sten, 14, FD_D2_O16_A6 * inv2_d2); // i1, i2+6, i3
			sca_stencil_set_location(sten, 15, 0,  7, 0, 0); sca_stencil_set_factor(sten, 15, FD_D2_O16_A7 * inv2_d2); // i1, i2+7, i3
			sca_stencil_set_location(sten, 16, 0,  8, 0, 0); sca_stencil_set_factor(sten, 16, FD_D2_O16_A8 * inv2_d2); // i1, i2+8, i3
		}

		Myfloat * const w = Wgrid.grid_3d ;
		Myfloat * const u = this->grid_3d ;

		// sca_code_create
		Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
		getGridIndex(pType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End) ;
		for (Myint i = 0; i < nPtPerStencil; i++) sca_stencil_set_input_array(sten, i, 1, n1, n2, n3, &u[i1Start + i2Start*n1 + i3Start*n2*n1]);
		sca_stencil_set_output_array(sten, 1, n1, n2, n3, &w[i1Start + i2Start*n1 + i3Start*n2*n1]);
		sca_code_create(&code_FD_D2_N2, sten, i1End-i1Start+1, i2End-i2Start+1, i3End-i3Start+1, 1);
		sca_stencil_destroy(sten);
	}

	flag_code_FD_D2_N2 = true ;
	return(RTN_CODE_OK) ;

#else
	printError("Grid_NEC_SCA::initialize_code_FD_D2_N2, not supported on your platform") ;
	return(RTN_CODE_KO) ;
#endif

	printDebug(MID_DEBUG, "OUT Grid_NEC_SCA::initialize_code_FD_D2_N2");
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_NEC_SCA::initialize_code_FD_D2_N3(Point_type pType, const Grid& Wgrid, Myint fdOrder)
{
	printDebug(MID_DEBUG, "IN Grid_NEC_SCA::initialize_code_FD_D2_N3");

#ifdef __NEC__
	if (!flag_code_FD_D2_N3)
	{
		sca_stencil_t sten ;
		sca_stencil_create(&sten);
		Myint nPtPerStencil = getPtPerStencilFD_D2(fdOrder) ;
		sca_stencil_append_elements(sten, nPtPerStencil);

		const Myfloat inv_d1  = Myfloat(1.0) / d1 ;
		const Myfloat inv_d2  = Myfloat(1.0) / d2 ;
		const Myfloat inv_d3  = Myfloat(1.0) / d3 ;

		const Myfloat inv2_d1 = inv_d1 * inv_d1 ;
		const Myfloat inv2_d2 = inv_d2 * inv_d2 ;
		const Myfloat inv2_d3 = inv_d3 * inv_d3 ;

		// sca_stencil_set_location for all points in stencil
		if (fdOrder == 2)
		{
			sca_stencil_set_location(sten, 0, 0, 0, -1, 0); sca_stencil_set_factor(sten, 0, FD_D2_O2_A1 * inv2_d3); // i1, i2, i3-1
			sca_stencil_set_location(sten, 1, 0, 0,  0, 0); sca_stencil_set_factor(sten, 1, FD_D2_O2_A0 * inv2_d3); // i1, i2, i3
			sca_stencil_set_location(sten, 2, 0, 0,  1, 0); sca_stencil_set_factor(sten, 2, FD_D2_O2_A1 * inv2_d3); // i1, i2, i3+1
		}
		else if (fdOrder == 4)
		{
			sca_stencil_set_location(sten, 0, 0, 0, -2, 0); sca_stencil_set_factor(sten, 0, FD_D2_O4_A2 * inv2_d3); // i1, i2, i3-2
			sca_stencil_set_location(sten, 1, 0, 0, -1, 0); sca_stencil_set_factor(sten, 1, FD_D2_O4_A1 * inv2_d3); // i1, i2, i3-1
			sca_stencil_set_location(sten, 2, 0, 0,  0, 0); sca_stencil_set_factor(sten, 2, FD_D2_O4_A0 * inv2_d3); // i1, i2, i3
			sca_stencil_set_location(sten, 3, 0, 0,  1, 0); sca_stencil_set_factor(sten, 3, FD_D2_O4_A1 * inv2_d3); // i1, i2, i3+1
			sca_stencil_set_location(sten, 4, 0, 0,  2, 0); sca_stencil_set_factor(sten, 4, FD_D2_O4_A2 * inv2_d3); // i1, i2, i3+2
		}
		else if (fdOrder == 8)
		{
			sca_stencil_set_location(sten, 0, 0, 0, -4, 0); sca_stencil_set_factor(sten, 0, FD_D2_O8_A4 * inv2_d3); // i1, i2, i3-4
			sca_stencil_set_location(sten, 1, 0, 0, -3, 0); sca_stencil_set_factor(sten, 1, FD_D2_O8_A3 * inv2_d3); // i1, i2, i3-3
			sca_stencil_set_location(sten, 2, 0, 0, -2, 0); sca_stencil_set_factor(sten, 2, FD_D2_O8_A2 * inv2_d3); // i1, i2, i3-2
			sca_stencil_set_location(sten, 3, 0, 0, -1, 0); sca_stencil_set_factor(sten, 3, FD_D2_O8_A1 * inv2_d3); // i1, i2, i3-1
			sca_stencil_set_location(sten, 4, 0, 0,  0, 0); sca_stencil_set_factor(sten, 4, FD_D2_O8_A0 * inv2_d3); // i1, i2, i3
			sca_stencil_set_location(sten, 5, 0, 0,  1, 0); sca_stencil_set_factor(sten, 5, FD_D2_O8_A1 * inv2_d3); // i1, i2, i3+1
			sca_stencil_set_location(sten, 6, 0, 0,  2, 0); sca_stencil_set_factor(sten, 6, FD_D2_O8_A2 * inv2_d3); // i1, i2, i3+2
			sca_stencil_set_location(sten, 7, 0, 0,  3, 0); sca_stencil_set_factor(sten, 7, FD_D2_O8_A3 * inv2_d3); // i1, i2, i3+3
			sca_stencil_set_location(sten, 8, 0, 0,  4, 0); sca_stencil_set_factor(sten, 8, FD_D2_O8_A4 * inv2_d3); // i1, i2, i3+4
		}
		else if (fdOrder == 12)
		{
			sca_stencil_set_location(sten,  0, 0, 0, -6, 0); sca_stencil_set_factor(sten,  0, FD_D2_O12_A6 * inv2_d3); // i1, i2, i3-6
			sca_stencil_set_location(sten,  1, 0, 0, -5, 0); sca_stencil_set_factor(sten,  1, FD_D2_O12_A5 * inv2_d3); // i1, i2, i3-5
			sca_stencil_set_location(sten,  2, 0, 0, -4, 0); sca_stencil_set_factor(sten,  2, FD_D2_O12_A4 * inv2_d3); // i1, i2, i3-4
			sca_stencil_set_location(sten,  3, 0, 0, -3, 0); sca_stencil_set_factor(sten,  3, FD_D2_O12_A3 * inv2_d3); // i1, i2, i3-3
			sca_stencil_set_location(sten,  4, 0, 0, -2, 0); sca_stencil_set_factor(sten,  4, FD_D2_O12_A2 * inv2_d3); // i1, i2, i3-2
			sca_stencil_set_location(sten,  5, 0, 0, -1, 0); sca_stencil_set_factor(sten,  5, FD_D2_O12_A1 * inv2_d3); // i1, i2, i3-1
			sca_stencil_set_location(sten,  6, 0, 0,  0, 0); sca_stencil_set_factor(sten,  6, FD_D2_O12_A0 * inv2_d3); // i1, i2, i3
			sca_stencil_set_location(sten,  7, 0, 0,  1, 0); sca_stencil_set_factor(sten,  7, FD_D2_O12_A1 * inv2_d3); // i1, i2, i3+1
			sca_stencil_set_location(sten,  8, 0, 0,  2, 0); sca_stencil_set_factor(sten,  8, FD_D2_O12_A2 * inv2_d3); // i1, i2, i3+2
			sca_stencil_set_location(sten,  9, 0, 0,  3, 0); sca_stencil_set_factor(sten,  9, FD_D2_O12_A3 * inv2_d3); // i1, i2, i3+3
			sca_stencil_set_location(sten, 10, 0, 0,  4, 0); sca_stencil_set_factor(sten, 10, FD_D2_O12_A4 * inv2_d3); // i1, i2, i3+4
			sca_stencil_set_location(sten, 11, 0, 0,  5, 0); sca_stencil_set_factor(sten, 11, FD_D2_O12_A5 * inv2_d3); // i1, i2, i3+5
			sca_stencil_set_location(sten, 12, 0, 0,  6, 0); sca_stencil_set_factor(sten, 12, FD_D2_O12_A6 * inv2_d3); // i1, i2, i3+6
		}
		else if (fdOrder == 16)
		{
			sca_stencil_set_location(sten,  0, 0, 0, -8, 0); sca_stencil_set_factor(sten,  0, FD_D2_O16_A8 * inv2_d3); // i1, i2, i3-8
			sca_stencil_set_location(sten,  1, 0, 0, -7, 0); sca_stencil_set_factor(sten,  1, FD_D2_O16_A7 * inv2_d3); // i1, i2, i3-7
			sca_stencil_set_location(sten,  2, 0, 0, -6, 0); sca_stencil_set_factor(sten,  2, FD_D2_O16_A6 * inv2_d3); // i1, i2, i3-6
			sca_stencil_set_location(sten,  3, 0, 0, -5, 0); sca_stencil_set_factor(sten,  3, FD_D2_O16_A5 * inv2_d3); // i1, i2, i3-5
			sca_stencil_set_location(sten,  4, 0, 0, -4, 0); sca_stencil_set_factor(sten,  4, FD_D2_O16_A4 * inv2_d3); // i1, i2, i3-4
			sca_stencil_set_location(sten,  5, 0, 0, -3, 0); sca_stencil_set_factor(sten,  5, FD_D2_O16_A3 * inv2_d3); // i1, i2, i3-3
			sca_stencil_set_location(sten,  6, 0, 0, -2, 0); sca_stencil_set_factor(sten,  6, FD_D2_O16_A2 * inv2_d3); // i1, i2, i3-2
			sca_stencil_set_location(sten,  7, 0, 0, -1, 0); sca_stencil_set_factor(sten,  7, FD_D2_O16_A1 * inv2_d3); // i1, i2, i3-1
			sca_stencil_set_location(sten,  8, 0, 0,  0, 0); sca_stencil_set_factor(sten,  8, FD_D2_O16_A0 * inv2_d3); // i1, i2, i3
			sca_stencil_set_location(sten,  9, 0, 0,  1, 0); sca_stencil_set_factor(sten,  9, FD_D2_O16_A1 * inv2_d3); // i1, i2, i3+1
			sca_stencil_set_location(sten, 10, 0, 0,  2, 0); sca_stencil_set_factor(sten, 10, FD_D2_O16_A2 * inv2_d3); // i1, i2, i3+2
			sca_stencil_set_location(sten, 11, 0, 0,  3, 0); sca_stencil_set_factor(sten, 11, FD_D2_O16_A3 * inv2_d3); // i1, i2, i3+3
			sca_stencil_set_location(sten, 12, 0, 0,  4, 0); sca_stencil_set_factor(sten, 12, FD_D2_O16_A4 * inv2_d3); // i1, i2, i3+4
			sca_stencil_set_location(sten, 13, 0, 0,  5, 0); sca_stencil_set_factor(sten, 13, FD_D2_O16_A5 * inv2_d3); // i1, i2, i3+5
			sca_stencil_set_location(sten, 14, 0, 0,  6, 0); sca_stencil_set_factor(sten, 14, FD_D2_O16_A6 * inv2_d3); // i1, i2, i3+6
			sca_stencil_set_location(sten, 15, 0, 0,  7, 0); sca_stencil_set_factor(sten, 15, FD_D2_O16_A7 * inv2_d3); // i1, i2, i3+7
			sca_stencil_set_location(sten, 16, 0, 0,  8, 0); sca_stencil_set_factor(sten, 16, FD_D2_O16_A8 * inv2_d3); // i1, i2, i3+8
		}

		Myfloat * const w = Wgrid.grid_3d ;
		Myfloat * const u = this->grid_3d ;

		// sca_code_create
		Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
		getGridIndex(pType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End) ;
		for (Myint i = 0; i < nPtPerStencil; i++) sca_stencil_set_input_array(sten, i, 1, n1, n2, n3, &u[i1Start + i2Start*n1 + i3Start*n2*n1]);
		sca_stencil_set_output_array(sten, 1, n1, n2, n3, &w[i1Start + i2Start*n1 + i3Start*n2*n1]);
		sca_code_create(&code_FD_D2_N3, sten, i1End-i1Start+1, i2End-i2Start+1, i3End-i3Start+1, 1);
		sca_stencil_destroy(sten);
	}

	flag_code_FD_D2_N3 = true ;
	return(RTN_CODE_OK) ;

#else
	printError("Grid_NEC_SCA::initialize_code_FD_D2_N3, not supported on your platform") ;
	return(RTN_CODE_KO) ;
#endif

	printDebug(MID_DEBUG, "OUT Grid_NEC_SCA::initialize_code_FD_D2_N3");
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_NEC_SCA::initialize_code_FD_LAPLACIAN(Point_type pType, const Grid& Wgrid, Myint fdOrder)
{
	printDebug(MID_DEBUG, "IN Grid_NEC_SCA::initialize_code_FD_LAPLACIAN");

#ifdef __NEC__
	if (!flag_code_FD_LAPLACIAN)
	{
		sca_stencil_t sten ;
		sca_stencil_create(&sten);
		Myint nPtPerStencil = getPtPerStencilFD_LAPLACIAN(fdOrder) ;
		sca_stencil_append_elements(sten, nPtPerStencil);

		const Myfloat inv_d1  = Myfloat(1.0) / d1 ;
		const Myfloat inv_d2  = Myfloat(1.0) / d2 ;
		const Myfloat inv_d3  = Myfloat(1.0) / d3 ;

		const Myfloat inv2_d1 = inv_d1 * inv_d1 ;
		const Myfloat inv2_d2 = inv_d2 * inv_d2 ;
		const Myfloat inv2_d3 = inv_d3 * inv_d3 ;

		// sca_stencil_set_location for all points in stencil
		if (dim == DIM1)
		{
			// same as in initialize_code_FD_D2_N1
			if (fdOrder == 2)
			{
				sca_stencil_set_location(sten, 0, -1, 0, 0, 0); sca_stencil_set_factor(sten, 0, FD_D2_O2_A1 * inv2_d1); // i1-1, i2, i3
				sca_stencil_set_location(sten, 1,  0, 0, 0, 0); sca_stencil_set_factor(sten, 1, FD_D2_O2_A0 * inv2_d1); // i1  , i2, i3
				sca_stencil_set_location(sten, 2,  1, 0, 0, 0); sca_stencil_set_factor(sten, 2, FD_D2_O2_A1 * inv2_d1); // i1+1, i2, i3
			}
			else if (fdOrder == 4)
			{
				sca_stencil_set_location(sten, 0, -2, 0, 0, 0); sca_stencil_set_factor(sten, 0, FD_D2_O4_A2 * inv2_d1); // i1-2, i2, i3
				sca_stencil_set_location(sten, 1, -1, 0, 0, 0); sca_stencil_set_factor(sten, 1, FD_D2_O4_A1 * inv2_d1); // i1-1, i2, i3
				sca_stencil_set_location(sten, 2,  0, 0, 0, 0); sca_stencil_set_factor(sten, 2, FD_D2_O4_A0 * inv2_d1); // i1  , i2, i3
				sca_stencil_set_location(sten, 3,  1, 0, 0, 0); sca_stencil_set_factor(sten, 3, FD_D2_O4_A1 * inv2_d1); // i1+1, i2, i3
				sca_stencil_set_location(sten, 4,  2, 0, 0, 0); sca_stencil_set_factor(sten, 4, FD_D2_O4_A2 * inv2_d1); // i1+2, i2, i3
			}
			else if (fdOrder == 8)
			{
				sca_stencil_set_location(sten, 0, -4, 0, 0, 0); sca_stencil_set_factor(sten, 0, FD_D2_O8_A4 * inv2_d1); // i1-4, i2, i3
				sca_stencil_set_location(sten, 1, -3, 0, 0, 0); sca_stencil_set_factor(sten, 1, FD_D2_O8_A3 * inv2_d1); // i1-3, i2, i3
				sca_stencil_set_location(sten, 2, -2, 0, 0, 0); sca_stencil_set_factor(sten, 2, FD_D2_O8_A2 * inv2_d1); // i1-2, i2, i3
				sca_stencil_set_location(sten, 3, -1, 0, 0, 0); sca_stencil_set_factor(sten, 3, FD_D2_O8_A1 * inv2_d1); // i1-1, i2, i3
				sca_stencil_set_location(sten, 4,  0, 0, 0, 0); sca_stencil_set_factor(sten, 4, FD_D2_O8_A0 * inv2_d1); // i1  , i2, i3
				sca_stencil_set_location(sten, 5,  1, 0, 0, 0); sca_stencil_set_factor(sten, 5, FD_D2_O8_A1 * inv2_d1); // i1+1, i2, i3
				sca_stencil_set_location(sten, 6,  2, 0, 0, 0); sca_stencil_set_factor(sten, 6, FD_D2_O8_A2 * inv2_d1); // i1+2, i2, i3
				sca_stencil_set_location(sten, 7,  3, 0, 0, 0); sca_stencil_set_factor(sten, 7, FD_D2_O8_A3 * inv2_d1); // i1+3, i2, i3
				sca_stencil_set_location(sten, 8,  4, 0, 0, 0); sca_stencil_set_factor(sten, 8, FD_D2_O8_A4 * inv2_d1); // i1+4, i2, i3
			}
			else if (fdOrder == 12)
			{
				sca_stencil_set_location(sten,  0, -6, 0, 0, 0); sca_stencil_set_factor(sten,  0, FD_D2_O12_A6 * inv2_d1); // i1-6, i2, i3
				sca_stencil_set_location(sten,  1, -5, 0, 0, 0); sca_stencil_set_factor(sten,  1, FD_D2_O12_A5 * inv2_d1); // i1-5, i2, i3
				sca_stencil_set_location(sten,  2, -4, 0, 0, 0); sca_stencil_set_factor(sten,  2, FD_D2_O12_A4 * inv2_d1); // i1-4, i2, i3
				sca_stencil_set_location(sten,  3, -3, 0, 0, 0); sca_stencil_set_factor(sten,  3, FD_D2_O12_A3 * inv2_d1); // i1-3, i2, i3
				sca_stencil_set_location(sten,  4, -2, 0, 0, 0); sca_stencil_set_factor(sten,  4, FD_D2_O12_A2 * inv2_d1); // i1-2, i2, i3
				sca_stencil_set_location(sten,  5, -1, 0, 0, 0); sca_stencil_set_factor(sten,  5, FD_D2_O12_A1 * inv2_d1); // i1-1, i2, i3
				sca_stencil_set_location(sten,  6,  0, 0, 0, 0); sca_stencil_set_factor(sten,  6, FD_D2_O12_A0 * inv2_d1); // i1  , i2, i3
				sca_stencil_set_location(sten,  7,  1, 0, 0, 0); sca_stencil_set_factor(sten,  7, FD_D2_O12_A1 * inv2_d1); // i1+1, i2, i3
				sca_stencil_set_location(sten,  8,  2, 0, 0, 0); sca_stencil_set_factor(sten,  8, FD_D2_O12_A2 * inv2_d1); // i1+2, i2, i3
				sca_stencil_set_location(sten,  9,  3, 0, 0, 0); sca_stencil_set_factor(sten,  9, FD_D2_O12_A3 * inv2_d1); // i1+3, i2, i3
				sca_stencil_set_location(sten, 10,  4, 0, 0, 0); sca_stencil_set_factor(sten, 10, FD_D2_O12_A4 * inv2_d1); // i1+4, i2, i3
				sca_stencil_set_location(sten, 11,  5, 0, 0, 0); sca_stencil_set_factor(sten, 11, FD_D2_O12_A5 * inv2_d1); // i1+5, i2, i3
				sca_stencil_set_location(sten, 12,  6, 0, 0, 0); sca_stencil_set_factor(sten, 12, FD_D2_O12_A6 * inv2_d1); // i1+6, i2, i3
			}
			else if (fdOrder == 16)
			{
				sca_stencil_set_location(sten,  0, -8, 0, 0, 0); sca_stencil_set_factor(sten,  0, FD_D2_O16_A8 * inv2_d1); // i1-8, i2, i3
				sca_stencil_set_location(sten,  1, -7, 0, 0, 0); sca_stencil_set_factor(sten,  1, FD_D2_O16_A7 * inv2_d1); // i1-7, i2, i3
				sca_stencil_set_location(sten,  2, -6, 0, 0, 0); sca_stencil_set_factor(sten,  2, FD_D2_O16_A6 * inv2_d1); // i1-6, i2, i3
				sca_stencil_set_location(sten,  3, -5, 0, 0, 0); sca_stencil_set_factor(sten,  3, FD_D2_O16_A5 * inv2_d1); // i1-5, i2, i3
				sca_stencil_set_location(sten,  4, -4, 0, 0, 0); sca_stencil_set_factor(sten,  4, FD_D2_O16_A4 * inv2_d1); // i1-4, i2, i3
				sca_stencil_set_location(sten,  5, -3, 0, 0, 0); sca_stencil_set_factor(sten,  5, FD_D2_O16_A3 * inv2_d1); // i1-3, i2, i3
				sca_stencil_set_location(sten,  6, -2, 0, 0, 0); sca_stencil_set_factor(sten,  6, FD_D2_O16_A2 * inv2_d1); // i1-2, i2, i3
				sca_stencil_set_location(sten,  7, -1, 0, 0, 0); sca_stencil_set_factor(sten,  7, FD_D2_O16_A1 * inv2_d1); // i1-1, i2, i3
				sca_stencil_set_location(sten,  8,  0, 0, 0, 0); sca_stencil_set_factor(sten,  8, FD_D2_O16_A0 * inv2_d1); // i1  , i2, i3
				sca_stencil_set_location(sten,  9,  1, 0, 0, 0); sca_stencil_set_factor(sten,  9, FD_D2_O16_A1 * inv2_d1); // i1+1, i2, i3
				sca_stencil_set_location(sten, 10,  2, 0, 0, 0); sca_stencil_set_factor(sten, 10, FD_D2_O16_A2 * inv2_d1); // i1+2, i2, i3
				sca_stencil_set_location(sten, 11,  3, 0, 0, 0); sca_stencil_set_factor(sten, 11, FD_D2_O16_A3 * inv2_d1); // i1+3, i2, i3
				sca_stencil_set_location(sten, 12,  4, 0, 0, 0); sca_stencil_set_factor(sten, 12, FD_D2_O16_A4 * inv2_d1); // i1+4, i2, i3
				sca_stencil_set_location(sten, 13,  5, 0, 0, 0); sca_stencil_set_factor(sten, 13, FD_D2_O16_A5 * inv2_d1); // i1+5, i2, i3
				sca_stencil_set_location(sten, 14,  6, 0, 0, 0); sca_stencil_set_factor(sten, 14, FD_D2_O16_A6 * inv2_d1); // i1+6, i2, i3
				sca_stencil_set_location(sten, 15,  7, 0, 0, 0); sca_stencil_set_factor(sten, 15, FD_D2_O16_A7 * inv2_d1); // i1+7, i2, i3
				sca_stencil_set_location(sten, 16,  8, 0, 0, 0); sca_stencil_set_factor(sten, 16, FD_D2_O16_A8 * inv2_d1); // i1+8, i2, i3
			}
		}

		else if (dim == DIM2)
		{
			if (fdOrder == 2)
			{
				sca_stencil_set_location(sten, 0, -1,  0, 0, 0); sca_stencil_set_factor(sten, 0, FD_D2_O2_A1 * inv2_d1); // i1-1, i2,   i3
				sca_stencil_set_location(sten, 1,  1,  0, 0, 0); sca_stencil_set_factor(sten, 1, FD_D2_O2_A1 * inv2_d1); // i1+1, i2,   i3
				sca_stencil_set_location(sten, 2,  0, -1, 0, 0); sca_stencil_set_factor(sten, 2, FD_D2_O2_A1 * inv2_d2); // i1,   i2-1, i3
				sca_stencil_set_location(sten, 3,  0,  1, 0, 0); sca_stencil_set_factor(sten, 3, FD_D2_O2_A1 * inv2_d2); // i1,   i2+1, i3
				sca_stencil_set_location(sten, 4,  0,  0, 0, 0); sca_stencil_set_factor(sten, 4, FD_D2_O2_A0 * inv2_d1 + FD_D2_O2_A0 * inv2_d2); // i1, i2, i3
			}
			else if (fdOrder == 4)
			{
				sca_stencil_set_location(sten, 0, -2,  0, 0, 0); sca_stencil_set_factor(sten, 0, FD_D2_O4_A2 * inv2_d1); // i1-2, i2,   i3
				sca_stencil_set_location(sten, 1, -1,  0, 0, 0); sca_stencil_set_factor(sten, 1, FD_D2_O4_A1 * inv2_d1); // i1-1, i2,   i3
				sca_stencil_set_location(sten, 2,  1,  0, 0, 0); sca_stencil_set_factor(sten, 2, FD_D2_O4_A1 * inv2_d1); // i1+1, i2,   i3
				sca_stencil_set_location(sten, 3,  2,  0, 0, 0); sca_stencil_set_factor(sten, 3, FD_D2_O4_A2 * inv2_d1); // i1+2, i2,   i3
				sca_stencil_set_location(sten, 4,  0, -2, 0, 0); sca_stencil_set_factor(sten, 4, FD_D2_O4_A2 * inv2_d2); // i1,   i2-2, i3
				sca_stencil_set_location(sten, 5,  0, -1, 0, 0); sca_stencil_set_factor(sten, 5, FD_D2_O4_A1 * inv2_d2); // i1,   i2-1, i3
				sca_stencil_set_location(sten, 6,  0,  1, 0, 0); sca_stencil_set_factor(sten, 6, FD_D2_O4_A1 * inv2_d2); // i1,   i2+1, i3
				sca_stencil_set_location(sten, 7,  0,  2, 0, 0); sca_stencil_set_factor(sten, 7, FD_D2_O4_A2 * inv2_d2); // i1,   i2+2, i3
				sca_stencil_set_location(sten, 8,  0,  0, 0, 0); sca_stencil_set_factor(sten, 8, FD_D2_O4_A0 * inv2_d1 + FD_D2_O4_A0 * inv2_d2); // i1, i2, i3
			}
			else if (fdOrder == 8)
			{
				sca_stencil_set_location(sten,  0, -4,  0, 0, 0); sca_stencil_set_factor(sten,  0, FD_D2_O8_A4 * inv2_d1); // i1-4, i2,   i3
				sca_stencil_set_location(sten,  1, -3,  0, 0, 0); sca_stencil_set_factor(sten,  1, FD_D2_O8_A3 * inv2_d1); // i1-3, i2,   i3
				sca_stencil_set_location(sten,  2, -2,  0, 0, 0); sca_stencil_set_factor(sten,  2, FD_D2_O8_A2 * inv2_d1); // i1-2, i2,   i3
				sca_stencil_set_location(sten,  3, -1,  0, 0, 0); sca_stencil_set_factor(sten,  3, FD_D2_O8_A1 * inv2_d1); // i1-1, i2,   i3
				sca_stencil_set_location(sten,  4,  1,  0, 0, 0); sca_stencil_set_factor(sten,  4, FD_D2_O8_A1 * inv2_d1); // i1+1, i2,   i3
				sca_stencil_set_location(sten,  5,  2,  0, 0, 0); sca_stencil_set_factor(sten,  5, FD_D2_O8_A2 * inv2_d1); // i1+2, i2,   i3
				sca_stencil_set_location(sten,  6,  3,  0, 0, 0); sca_stencil_set_factor(sten,  6, FD_D2_O8_A3 * inv2_d1); // i1+3, i2,   i3
				sca_stencil_set_location(sten,  7,  4,  0, 0, 0); sca_stencil_set_factor(sten,  7, FD_D2_O8_A4 * inv2_d1); // i1+4, i2,   i3
				sca_stencil_set_location(sten,  8,  0, -4, 0, 0); sca_stencil_set_factor(sten,  8, FD_D2_O8_A4 * inv2_d2); // i1,   i2-4, i3
				sca_stencil_set_location(sten,  9,  0, -3, 0, 0); sca_stencil_set_factor(sten,  9, FD_D2_O8_A3 * inv2_d2); // i1,   i2-3, i3
				sca_stencil_set_location(sten, 10,  0, -2, 0, 0); sca_stencil_set_factor(sten, 10, FD_D2_O8_A2 * inv2_d2); // i1,   i2-2, i3
				sca_stencil_set_location(sten, 11,  0, -1, 0, 0); sca_stencil_set_factor(sten, 11, FD_D2_O8_A1 * inv2_d2); // i1,   i2-1, i3
				sca_stencil_set_location(sten, 12,  0,  1, 0, 0); sca_stencil_set_factor(sten, 12, FD_D2_O8_A1 * inv2_d2); // i1,   i2+1, i3
				sca_stencil_set_location(sten, 13,  0,  2, 0, 0); sca_stencil_set_factor(sten, 13, FD_D2_O8_A2 * inv2_d2); // i1,   i2+2, i3
				sca_stencil_set_location(sten, 14,  0,  3, 0, 0); sca_stencil_set_factor(sten, 14, FD_D2_O8_A3 * inv2_d2); // i1,   i2+3, i3
				sca_stencil_set_location(sten, 15,  0,  4, 0, 0); sca_stencil_set_factor(sten, 15, FD_D2_O8_A4 * inv2_d2); // i1,   i2+4, i3
				sca_stencil_set_location(sten, 16,  0,  0, 0, 0); sca_stencil_set_factor(sten, 16, FD_D2_O8_A0 * inv2_d1 + FD_D2_O8_A0 * inv2_d2); // i1, i2, i3
			}
			else if (fdOrder == 12)
			{
				sca_stencil_set_location(sten,  0, -6,  0, 0, 0); sca_stencil_set_factor(sten,  0, FD_D2_O12_A6 * inv2_d1); // i1-6, i2,   i3
				sca_stencil_set_location(sten,  1, -5,  0, 0, 0); sca_stencil_set_factor(sten,  1, FD_D2_O12_A5 * inv2_d1); // i1-5, i2,   i3
				sca_stencil_set_location(sten,  2, -4,  0, 0, 0); sca_stencil_set_factor(sten,  2, FD_D2_O12_A4 * inv2_d1); // i1-4, i2,   i3
				sca_stencil_set_location(sten,  3, -3,  0, 0, 0); sca_stencil_set_factor(sten,  3, FD_D2_O12_A3 * inv2_d1); // i1-3, i2,   i3
				sca_stencil_set_location(sten,  4, -2,  0, 0, 0); sca_stencil_set_factor(sten,  4, FD_D2_O12_A2 * inv2_d1); // i1-2, i2,   i3
				sca_stencil_set_location(sten,  5, -1,  0, 0, 0); sca_stencil_set_factor(sten,  5, FD_D2_O12_A1 * inv2_d1); // i1-1, i2,   i3
				sca_stencil_set_location(sten,  6,  1,  0, 0, 0); sca_stencil_set_factor(sten,  6, FD_D2_O12_A1 * inv2_d1); // i1+1, i2,   i3
				sca_stencil_set_location(sten,  7,  2,  0, 0, 0); sca_stencil_set_factor(sten,  7, FD_D2_O12_A2 * inv2_d1); // i1+2, i2,   i3
				sca_stencil_set_location(sten,  8,  3,  0, 0, 0); sca_stencil_set_factor(sten,  8, FD_D2_O12_A3 * inv2_d1); // i1+3, i2,   i3
				sca_stencil_set_location(sten,  9,  4,  0, 0, 0); sca_stencil_set_factor(sten,  9, FD_D2_O12_A4 * inv2_d1); // i1+4, i2,   i3
				sca_stencil_set_location(sten, 10,  5,  0, 0, 0); sca_stencil_set_factor(sten, 10, FD_D2_O12_A5 * inv2_d1); // i1+5, i2,   i3
				sca_stencil_set_location(sten, 11,  6,  0, 0, 0); sca_stencil_set_factor(sten, 11, FD_D2_O12_A6 * inv2_d1); // i1+6, i2,   i3
				sca_stencil_set_location(sten, 12,  0, -6, 0, 0); sca_stencil_set_factor(sten, 12, FD_D2_O12_A6 * inv2_d2); // i1,   i2-6, i3
				sca_stencil_set_location(sten, 13,  0, -5, 0, 0); sca_stencil_set_factor(sten, 13, FD_D2_O12_A5 * inv2_d2); // i1,   i2-5, i3
				sca_stencil_set_location(sten, 14,  0, -4, 0, 0); sca_stencil_set_factor(sten, 14, FD_D2_O12_A4 * inv2_d2); // i1,   i2-4, i3
				sca_stencil_set_location(sten, 15,  0, -3, 0, 0); sca_stencil_set_factor(sten, 15, FD_D2_O12_A3 * inv2_d2); // i1,   i2-3, i3
				sca_stencil_set_location(sten, 16,  0, -2, 0, 0); sca_stencil_set_factor(sten, 16, FD_D2_O12_A2 * inv2_d2); // i1,   i2-2, i3
				sca_stencil_set_location(sten, 17,  0, -1, 0, 0); sca_stencil_set_factor(sten, 17, FD_D2_O12_A1 * inv2_d2); // i1,   i2-1, i3
				sca_stencil_set_location(sten, 18,  0,  1, 0, 0); sca_stencil_set_factor(sten, 18, FD_D2_O12_A1 * inv2_d2); // i1,   i2+1, i3
				sca_stencil_set_location(sten, 19,  0,  2, 0, 0); sca_stencil_set_factor(sten, 19, FD_D2_O12_A2 * inv2_d2); // i1,   i2+2, i3
				sca_stencil_set_location(sten, 20,  0,  3, 0, 0); sca_stencil_set_factor(sten, 20, FD_D2_O12_A3 * inv2_d2); // i1,   i2+3, i3
				sca_stencil_set_location(sten, 21,  0,  4, 0, 0); sca_stencil_set_factor(sten, 21, FD_D2_O12_A4 * inv2_d2); // i1,   i2+4, i3
				sca_stencil_set_location(sten, 22,  0,  5, 0, 0); sca_stencil_set_factor(sten, 22, FD_D2_O12_A5 * inv2_d2); // i1,   i2+5, i3
				sca_stencil_set_location(sten, 23,  0,  6, 0, 0); sca_stencil_set_factor(sten, 23, FD_D2_O12_A6 * inv2_d2); // i1,   i2+6, i3
				sca_stencil_set_location(sten, 24,  0,  0, 0, 0); sca_stencil_set_factor(sten, 24, FD_D2_O12_A0 * inv2_d1 + FD_D2_O12_A0 * inv2_d2); // i1, i2, i3
			}
			else if (fdOrder == 16)
			{
				sca_stencil_set_location(sten,  0, -8,  0, 0, 0); sca_stencil_set_factor(sten,  0, FD_D2_O16_A8 * inv2_d1); // i1-8, i2,   i3
				sca_stencil_set_location(sten,  1, -7,  0, 0, 0); sca_stencil_set_factor(sten,  1, FD_D2_O16_A7 * inv2_d1); // i1-7, i2,   i3
				sca_stencil_set_location(sten,  2, -6,  0, 0, 0); sca_stencil_set_factor(sten,  2, FD_D2_O16_A6 * inv2_d1); // i1-6, i2,   i3
				sca_stencil_set_location(sten,  3, -5,  0, 0, 0); sca_stencil_set_factor(sten,  3, FD_D2_O16_A5 * inv2_d1); // i1-5, i2,   i3
				sca_stencil_set_location(sten,  4, -4,  0, 0, 0); sca_stencil_set_factor(sten,  4, FD_D2_O16_A4 * inv2_d1); // i1-4, i2,   i3
				sca_stencil_set_location(sten,  5, -3,  0, 0, 0); sca_stencil_set_factor(sten,  5, FD_D2_O16_A3 * inv2_d1); // i1-3, i2,   i3
				sca_stencil_set_location(sten,  6, -2,  0, 0, 0); sca_stencil_set_factor(sten,  6, FD_D2_O16_A2 * inv2_d1); // i1-2, i2,   i3
				sca_stencil_set_location(sten,  7, -1,  0, 0, 0); sca_stencil_set_factor(sten,  7, FD_D2_O16_A1 * inv2_d1); // i1-1, i2,   i3
				sca_stencil_set_location(sten,  8,  1,  0, 0, 0); sca_stencil_set_factor(sten,  8, FD_D2_O16_A1 * inv2_d1); // i1+1, i2,   i3
				sca_stencil_set_location(sten,  9,  2,  0, 0, 0); sca_stencil_set_factor(sten,  9, FD_D2_O16_A2 * inv2_d1); // i1+2, i2,   i3
				sca_stencil_set_location(sten, 10,  3,  0, 0, 0); sca_stencil_set_factor(sten, 10, FD_D2_O16_A3 * inv2_d1); // i1+3, i2,   i3
				sca_stencil_set_location(sten, 11,  4,  0, 0, 0); sca_stencil_set_factor(sten, 11, FD_D2_O16_A4 * inv2_d1); // i1+4, i2,   i3
				sca_stencil_set_location(sten, 12,  5,  0, 0, 0); sca_stencil_set_factor(sten, 12, FD_D2_O16_A5 * inv2_d1); // i1+5, i2,   i3
				sca_stencil_set_location(sten, 13,  6,  0, 0, 0); sca_stencil_set_factor(sten, 13, FD_D2_O16_A6 * inv2_d1); // i1+6, i2,   i3
				sca_stencil_set_location(sten, 14,  7,  0, 0, 0); sca_stencil_set_factor(sten, 14, FD_D2_O16_A7 * inv2_d1); // i1+7, i2,   i3
				sca_stencil_set_location(sten, 15,  8,  0, 0, 0); sca_stencil_set_factor(sten, 15, FD_D2_O16_A8 * inv2_d1); // i1+8, i2,   i3
				sca_stencil_set_location(sten, 16,  0, -8, 0, 0); sca_stencil_set_factor(sten, 16, FD_D2_O16_A8 * inv2_d2); // i1,   i2-8, i3
				sca_stencil_set_location(sten, 17,  0, -7, 0, 0); sca_stencil_set_factor(sten, 17, FD_D2_O16_A7 * inv2_d2); // i1,   i2-7, i3
				sca_stencil_set_location(sten, 18,  0, -6, 0, 0); sca_stencil_set_factor(sten, 18, FD_D2_O16_A6 * inv2_d2); // i1,   i2-6, i3
				sca_stencil_set_location(sten, 19,  0, -5, 0, 0); sca_stencil_set_factor(sten, 19, FD_D2_O16_A5 * inv2_d2); // i1,   i2-5, i3
				sca_stencil_set_location(sten, 20,  0, -4, 0, 0); sca_stencil_set_factor(sten, 20, FD_D2_O16_A4 * inv2_d2); // i1,   i2-4, i3
				sca_stencil_set_location(sten, 21,  0, -3, 0, 0); sca_stencil_set_factor(sten, 21, FD_D2_O16_A3 * inv2_d2); // i1,   i2-3, i3
				sca_stencil_set_location(sten, 22,  0, -2, 0, 0); sca_stencil_set_factor(sten, 22, FD_D2_O16_A2 * inv2_d2); // i1,   i2-2, i3
				sca_stencil_set_location(sten, 23,  0, -1, 0, 0); sca_stencil_set_factor(sten, 23, FD_D2_O16_A1 * inv2_d2); // i1,   i2-1, i3
				sca_stencil_set_location(sten, 24,  0,  1, 0, 0); sca_stencil_set_factor(sten, 24, FD_D2_O16_A1 * inv2_d2); // i1,   i2+1, i3
				sca_stencil_set_location(sten, 25,  0,  2, 0, 0); sca_stencil_set_factor(sten, 25, FD_D2_O16_A2 * inv2_d2); // i1,   i2+2, i3
				sca_stencil_set_location(sten, 26,  0,  3, 0, 0); sca_stencil_set_factor(sten, 26, FD_D2_O16_A3 * inv2_d2); // i1,   i2+3, i3
				sca_stencil_set_location(sten, 27,  0,  4, 0, 0); sca_stencil_set_factor(sten, 27, FD_D2_O16_A4 * inv2_d2); // i1,   i2+4, i3
				sca_stencil_set_location(sten, 28,  0,  5, 0, 0); sca_stencil_set_factor(sten, 28, FD_D2_O16_A5 * inv2_d2); // i1,   i2+5, i3
				sca_stencil_set_location(sten, 29,  0,  6, 0, 0); sca_stencil_set_factor(sten, 29, FD_D2_O16_A6 * inv2_d2); // i1,   i2+6, i3
				sca_stencil_set_location(sten, 30,  0,  7, 0, 0); sca_stencil_set_factor(sten, 30, FD_D2_O16_A7 * inv2_d2); // i1,   i2+7, i3
				sca_stencil_set_location(sten, 31,  0,  8, 0, 0); sca_stencil_set_factor(sten, 31, FD_D2_O16_A8 * inv2_d2); // i1,   i2+8, i3
				sca_stencil_set_location(sten, 32,  0,  0, 0, 0); sca_stencil_set_factor(sten, 32, FD_D2_O16_A0 * inv2_d1 + FD_D2_O16_A0 * inv2_d2); // i1, i2, i3
			}
		}

		else if (dim == DIM3)
		{
			if (fdOrder == 2)
			{
				sca_stencil_set_location(sten, 0, -1,  0,  0, 0); sca_stencil_set_factor(sten, 0, FD_D2_O2_A1 * inv2_d1); // i1-1, i2,   i3
				sca_stencil_set_location(sten, 1,  1,  0,  0, 0); sca_stencil_set_factor(sten, 1, FD_D2_O2_A1 * inv2_d1); // i1+1, i2,   i3
				sca_stencil_set_location(sten, 2,  0, -1,  0, 0); sca_stencil_set_factor(sten, 2, FD_D2_O2_A1 * inv2_d2); // i1,   i2-1, i3
				sca_stencil_set_location(sten, 3,  0,  1,  0, 0); sca_stencil_set_factor(sten, 3, FD_D2_O2_A1 * inv2_d2); // i1,   i2+1, i3
				sca_stencil_set_location(sten, 4,  0,  0, -1, 0); sca_stencil_set_factor(sten, 4, FD_D2_O2_A1 * inv2_d3); // i1,   i2,   i3-1
				sca_stencil_set_location(sten, 5,  0,  0,  1, 0); sca_stencil_set_factor(sten, 5, FD_D2_O2_A1 * inv2_d3); // i1,   i2,   i3+1
				sca_stencil_set_location(sten, 6,  0,  0,  0, 0); sca_stencil_set_factor(sten, 6, FD_D2_O2_A0 * inv2_d1 + FD_D2_O2_A0 * inv2_d2 + FD_D2_O2_A0 * inv2_d3); // i1, i2, i3
			}
			else if (fdOrder == 4)
			{
				sca_stencil_set_location(sten,  0, -2,  0,  0, 0); sca_stencil_set_factor(sten,  0, FD_D2_O4_A2 * inv2_d1); // i1-2, i2,   i3
				sca_stencil_set_location(sten,  1, -1,  0,  0, 0); sca_stencil_set_factor(sten,  1, FD_D2_O4_A1 * inv2_d1); // i1-1, i2,   i3
				sca_stencil_set_location(sten,  2,  1,  0,  0, 0); sca_stencil_set_factor(sten,  2, FD_D2_O4_A1 * inv2_d1); // i1+1, i2,   i3
				sca_stencil_set_location(sten,  3,  2,  0,  0, 0); sca_stencil_set_factor(sten,  3, FD_D2_O4_A2 * inv2_d1); // i1+2, i2,   i3
				sca_stencil_set_location(sten,  4,  0, -2,  0, 0); sca_stencil_set_factor(sten,  4, FD_D2_O4_A2 * inv2_d2); // i1,   i2-2, i3
				sca_stencil_set_location(sten,  5,  0, -1,  0, 0); sca_stencil_set_factor(sten,  5, FD_D2_O4_A1 * inv2_d2); // i1,   i2-1, i3
				sca_stencil_set_location(sten,  6,  0,  1,  0, 0); sca_stencil_set_factor(sten,  6, FD_D2_O4_A1 * inv2_d2); // i1,   i2+1, i3
				sca_stencil_set_location(sten,  7,  0,  2,  0, 0); sca_stencil_set_factor(sten,  7, FD_D2_O4_A2 * inv2_d2); // i1,   i2+2, i3
				sca_stencil_set_location(sten,  8,  0,  0, -2, 0); sca_stencil_set_factor(sten,  8, FD_D2_O4_A2 * inv2_d3); // i1,   i2,   i3-2
				sca_stencil_set_location(sten,  9,  0,  0, -1, 0); sca_stencil_set_factor(sten,  9, FD_D2_O4_A1 * inv2_d3); // i1,   i2,   i3-1
				sca_stencil_set_location(sten, 10,  0,  0,  1, 0); sca_stencil_set_factor(sten, 10, FD_D2_O4_A1 * inv2_d3); // i1,   i2,   i3+1
				sca_stencil_set_location(sten, 11,  0,  0,  2, 0); sca_stencil_set_factor(sten, 11, FD_D2_O4_A2 * inv2_d3); // i1,   i2,   i3+2
				sca_stencil_set_location(sten, 12,  0,  0,  0, 0); sca_stencil_set_factor(sten, 12, FD_D2_O4_A0 * inv2_d1 + FD_D2_O4_A0 * inv2_d2 + FD_D2_O4_A0 * inv2_d3); // i1, i2, i3
			}
			else if (fdOrder == 8)
			{
				sca_stencil_set_location(sten,  0, -4,  0,  0, 0); sca_stencil_set_factor(sten,  0, FD_D2_O8_A4 * inv2_d1); // i1-4, i2,   i3
				sca_stencil_set_location(sten,  1, -3,  0,  0, 0); sca_stencil_set_factor(sten,  1, FD_D2_O8_A3 * inv2_d1); // i1-3, i2,   i3
				sca_stencil_set_location(sten,  2, -2,  0,  0, 0); sca_stencil_set_factor(sten,  2, FD_D2_O8_A2 * inv2_d1); // i1-2, i2,   i3
				sca_stencil_set_location(sten,  3, -1,  0,  0, 0); sca_stencil_set_factor(sten,  3, FD_D2_O8_A1 * inv2_d1); // i1-1, i2,   i3
				sca_stencil_set_location(sten,  4,  1,  0,  0, 0); sca_stencil_set_factor(sten,  4, FD_D2_O8_A1 * inv2_d1); // i1+1, i2,   i3
				sca_stencil_set_location(sten,  5,  2,  0,  0, 0); sca_stencil_set_factor(sten,  5, FD_D2_O8_A2 * inv2_d1); // i1+2, i2,   i3
				sca_stencil_set_location(sten,  6,  3,  0,  0, 0); sca_stencil_set_factor(sten,  6, FD_D2_O8_A3 * inv2_d1); // i1+3, i2,   i3
				sca_stencil_set_location(sten,  7,  4,  0,  0, 0); sca_stencil_set_factor(sten,  7, FD_D2_O8_A4 * inv2_d1); // i1+4, i2,   i3
				sca_stencil_set_location(sten,  8,  0, -4,  0, 0); sca_stencil_set_factor(sten,  8, FD_D2_O8_A4 * inv2_d2); // i1,   i2-4, i3
				sca_stencil_set_location(sten,  9,  0, -3,  0, 0); sca_stencil_set_factor(sten,  9, FD_D2_O8_A3 * inv2_d2); // i1,   i2-3, i3
				sca_stencil_set_location(sten, 10,  0, -2,  0, 0); sca_stencil_set_factor(sten, 10, FD_D2_O8_A2 * inv2_d2); // i1,   i2-2, i3
				sca_stencil_set_location(sten, 11,  0, -1,  0, 0); sca_stencil_set_factor(sten, 11, FD_D2_O8_A1 * inv2_d2); // i1,   i2-1, i3
				sca_stencil_set_location(sten, 12,  0,  1,  0, 0); sca_stencil_set_factor(sten, 12, FD_D2_O8_A1 * inv2_d2); // i1,   i2+1, i3
				sca_stencil_set_location(sten, 13,  0,  2,  0, 0); sca_stencil_set_factor(sten, 13, FD_D2_O8_A2 * inv2_d2); // i1,   i2+2, i3
				sca_stencil_set_location(sten, 14,  0,  3,  0, 0); sca_stencil_set_factor(sten, 14, FD_D2_O8_A3 * inv2_d2); // i1,   i2+3, i3
				sca_stencil_set_location(sten, 15,  0,  4,  0, 0); sca_stencil_set_factor(sten, 15, FD_D2_O8_A4 * inv2_d2); // i1,   i2+4, i3
				sca_stencil_set_location(sten, 16,  0,  0, -4, 0); sca_stencil_set_factor(sten, 16, FD_D2_O8_A4 * inv2_d3); // i1,   i2,   i3-4
				sca_stencil_set_location(sten, 17,  0,  0, -3, 0); sca_stencil_set_factor(sten, 17, FD_D2_O8_A3 * inv2_d3); // i1,   i2,   i3-3
				sca_stencil_set_location(sten, 18,  0,  0, -2, 0); sca_stencil_set_factor(sten, 18, FD_D2_O8_A2 * inv2_d3); // i1,   i2,   i3-2
				sca_stencil_set_location(sten, 19,  0,  0, -1, 0); sca_stencil_set_factor(sten, 19, FD_D2_O8_A1 * inv2_d3); // i1,   i2,   i3-1
				sca_stencil_set_location(sten, 20,  0,  0,  1, 0); sca_stencil_set_factor(sten, 20, FD_D2_O8_A1 * inv2_d3); // i1,   i2,   i3+1
				sca_stencil_set_location(sten, 21,  0,  0,  2, 0); sca_stencil_set_factor(sten, 21, FD_D2_O8_A2 * inv2_d3); // i1,   i2,   i3+2
				sca_stencil_set_location(sten, 22,  0,  0,  3, 0); sca_stencil_set_factor(sten, 22, FD_D2_O8_A3 * inv2_d3); // i1,   i2,   i3+3
				sca_stencil_set_location(sten, 23,  0,  0,  4, 0); sca_stencil_set_factor(sten, 23, FD_D2_O8_A4 * inv2_d3); // i1,   i2,   i3+4
				sca_stencil_set_location(sten, 24,  0,  0,  0, 0); sca_stencil_set_factor(sten, 24, FD_D2_O8_A0 * inv2_d1 + FD_D2_O8_A0 * inv2_d2 + FD_D2_O8_A0 * inv2_d3); // i1, i2, i3
			}
			else if (fdOrder == 12)
			{
				sca_stencil_set_location(sten,  0, -6,  0,  0, 0); sca_stencil_set_factor(sten,  0, FD_D2_O12_A6 * inv2_d1); // i1-6, i2,   i3
				sca_stencil_set_location(sten,  1, -5,  0,  0, 0); sca_stencil_set_factor(sten,  1, FD_D2_O12_A5 * inv2_d1); // i1-5, i2,   i3
				sca_stencil_set_location(sten,  2, -4,  0,  0, 0); sca_stencil_set_factor(sten,  2, FD_D2_O12_A4 * inv2_d1); // i1-4, i2,   i3
				sca_stencil_set_location(sten,  3, -3,  0,  0, 0); sca_stencil_set_factor(sten,  3, FD_D2_O12_A3 * inv2_d1); // i1-3, i2,   i3
				sca_stencil_set_location(sten,  4, -2,  0,  0, 0); sca_stencil_set_factor(sten,  4, FD_D2_O12_A2 * inv2_d1); // i1-2, i2,   i3
				sca_stencil_set_location(sten,  5, -1,  0,  0, 0); sca_stencil_set_factor(sten,  5, FD_D2_O12_A1 * inv2_d1); // i1-1, i2,   i3
				sca_stencil_set_location(sten,  6,  1,  0,  0, 0); sca_stencil_set_factor(sten,  6, FD_D2_O12_A1 * inv2_d1); // i1+1, i2,   i3
				sca_stencil_set_location(sten,  7,  2,  0,  0, 0); sca_stencil_set_factor(sten,  7, FD_D2_O12_A2 * inv2_d1); // i1+2, i2,   i3
				sca_stencil_set_location(sten,  8,  3,  0,  0, 0); sca_stencil_set_factor(sten,  8, FD_D2_O12_A3 * inv2_d1); // i1+3, i2,   i3
				sca_stencil_set_location(sten,  9,  4,  0,  0, 0); sca_stencil_set_factor(sten,  9, FD_D2_O12_A4 * inv2_d1); // i1+4, i2,   i3
				sca_stencil_set_location(sten, 10,  5,  0,  0, 0); sca_stencil_set_factor(sten, 10, FD_D2_O12_A5 * inv2_d1); // i1+5, i2,   i3
				sca_stencil_set_location(sten, 11,  6,  0,  0, 0); sca_stencil_set_factor(sten, 11, FD_D2_O12_A6 * inv2_d1); // i1+6, i2,   i3
				sca_stencil_set_location(sten, 12,  0, -6,  0, 0); sca_stencil_set_factor(sten, 12, FD_D2_O12_A6 * inv2_d2); // i1,   i2-6, i3
				sca_stencil_set_location(sten, 13,  0, -5,  0, 0); sca_stencil_set_factor(sten, 13, FD_D2_O12_A5 * inv2_d2); // i1,   i2-5, i3
				sca_stencil_set_location(sten, 14,  0, -4,  0, 0); sca_stencil_set_factor(sten, 14, FD_D2_O12_A4 * inv2_d2); // i1,   i2-4, i3
				sca_stencil_set_location(sten, 15,  0, -3,  0, 0); sca_stencil_set_factor(sten, 15, FD_D2_O12_A3 * inv2_d2); // i1,   i2-3, i3
				sca_stencil_set_location(sten, 16,  0, -2,  0, 0); sca_stencil_set_factor(sten, 16, FD_D2_O12_A2 * inv2_d2); // i1,   i2-2, i3
				sca_stencil_set_location(sten, 17,  0, -1,  0, 0); sca_stencil_set_factor(sten, 17, FD_D2_O12_A1 * inv2_d2); // i1,   i2-1, i3
				sca_stencil_set_location(sten, 18,  0,  1,  0, 0); sca_stencil_set_factor(sten, 18, FD_D2_O12_A1 * inv2_d2); // i1,   i2+1, i3
				sca_stencil_set_location(sten, 19,  0,  2,  0, 0); sca_stencil_set_factor(sten, 19, FD_D2_O12_A2 * inv2_d2); // i1,   i2+2, i3
				sca_stencil_set_location(sten, 20,  0,  3,  0, 0); sca_stencil_set_factor(sten, 20, FD_D2_O12_A3 * inv2_d2); // i1,   i2+3, i3
				sca_stencil_set_location(sten, 21,  0,  4,  0, 0); sca_stencil_set_factor(sten, 21, FD_D2_O12_A4 * inv2_d2); // i1,   i2+4, i3
				sca_stencil_set_location(sten, 22,  0,  5,  0, 0); sca_stencil_set_factor(sten, 22, FD_D2_O12_A5 * inv2_d2); // i1,   i2+5, i3
				sca_stencil_set_location(sten, 23,  0,  6,  0, 0); sca_stencil_set_factor(sten, 23, FD_D2_O12_A6 * inv2_d2); // i1,   i2+6, i3
				sca_stencil_set_location(sten, 24,  0,  0, -6, 0); sca_stencil_set_factor(sten, 24, FD_D2_O12_A6 * inv2_d3); // i1,   i2,   i3-6
				sca_stencil_set_location(sten, 25,  0,  0, -5, 0); sca_stencil_set_factor(sten, 25, FD_D2_O12_A5 * inv2_d3); // i1,   i2,   i3-5
				sca_stencil_set_location(sten, 26,  0,  0, -4, 0); sca_stencil_set_factor(sten, 26, FD_D2_O12_A4 * inv2_d3); // i1,   i2,   i3-4
				sca_stencil_set_location(sten, 27,  0,  0, -3, 0); sca_stencil_set_factor(sten, 27, FD_D2_O12_A3 * inv2_d3); // i1,   i2,   i3-3
				sca_stencil_set_location(sten, 28,  0,  0, -2, 0); sca_stencil_set_factor(sten, 28, FD_D2_O12_A2 * inv2_d3); // i1,   i2,   i3-2
				sca_stencil_set_location(sten, 29,  0,  0, -1, 0); sca_stencil_set_factor(sten, 29, FD_D2_O12_A1 * inv2_d3); // i1,   i2,   i3-1
				sca_stencil_set_location(sten, 30,  0,  0,  1, 0); sca_stencil_set_factor(sten, 30, FD_D2_O12_A1 * inv2_d3); // i1,   i2,   i3+1
				sca_stencil_set_location(sten, 31,  0,  0,  2, 0); sca_stencil_set_factor(sten, 31, FD_D2_O12_A2 * inv2_d3); // i1,   i2,   i3+2
				sca_stencil_set_location(sten, 32,  0,  0,  3, 0); sca_stencil_set_factor(sten, 32, FD_D2_O12_A3 * inv2_d3); // i1,   i2,   i3+3
				sca_stencil_set_location(sten, 33,  0,  0,  4, 0); sca_stencil_set_factor(sten, 33, FD_D2_O12_A4 * inv2_d3); // i1,   i2,   i3+4
				sca_stencil_set_location(sten, 34,  0,  0,  5, 0); sca_stencil_set_factor(sten, 34, FD_D2_O12_A5 * inv2_d3); // i1,   i2,   i3+5
				sca_stencil_set_location(sten, 35,  0,  0,  6, 0); sca_stencil_set_factor(sten, 35, FD_D2_O12_A6 * inv2_d3); // i1,   i2,   i3+6
				sca_stencil_set_location(sten, 36,  0,  0,  0, 0); sca_stencil_set_factor(sten, 36, FD_D2_O12_A0 * inv2_d1 + FD_D2_O12_A0 * inv2_d2 + FD_D2_O12_A0 * inv2_d3); // i1, i2, i3
			}
			else if (fdOrder == 16)
			{
				sca_stencil_set_location(sten,  0, -8,  0,  0, 0); sca_stencil_set_factor(sten,  0, FD_D2_O16_A8 * inv2_d1); // i1-8, i2,   i3
				sca_stencil_set_location(sten,  1, -7,  0,  0, 0); sca_stencil_set_factor(sten,  1, FD_D2_O16_A7 * inv2_d1); // i1-7, i2,   i3
				sca_stencil_set_location(sten,  2, -6,  0,  0, 0); sca_stencil_set_factor(sten,  2, FD_D2_O16_A6 * inv2_d1); // i1-6, i2,   i3
				sca_stencil_set_location(sten,  3, -5,  0,  0, 0); sca_stencil_set_factor(sten,  3, FD_D2_O16_A5 * inv2_d1); // i1-5, i2,   i3
				sca_stencil_set_location(sten,  4, -4,  0,  0, 0); sca_stencil_set_factor(sten,  4, FD_D2_O16_A4 * inv2_d1); // i1-4, i2,   i3
				sca_stencil_set_location(sten,  5, -3,  0,  0, 0); sca_stencil_set_factor(sten,  5, FD_D2_O16_A3 * inv2_d1); // i1-3, i2,   i3
				sca_stencil_set_location(sten,  6, -2,  0,  0, 0); sca_stencil_set_factor(sten,  6, FD_D2_O16_A2 * inv2_d1); // i1-2, i2,   i3
				sca_stencil_set_location(sten,  7, -1,  0,  0, 0); sca_stencil_set_factor(sten,  7, FD_D2_O16_A1 * inv2_d1); // i1-1, i2,   i3
				sca_stencil_set_location(sten,  8,  1,  0,  0, 0); sca_stencil_set_factor(sten,  8, FD_D2_O16_A1 * inv2_d1); // i1+1, i2,   i3
				sca_stencil_set_location(sten,  9,  2,  0,  0, 0); sca_stencil_set_factor(sten,  9, FD_D2_O16_A2 * inv2_d1); // i1+2, i2,   i3
				sca_stencil_set_location(sten, 10,  3,  0,  0, 0); sca_stencil_set_factor(sten, 10, FD_D2_O16_A3 * inv2_d1); // i1+3, i2,   i3
				sca_stencil_set_location(sten, 11,  4,  0,  0, 0); sca_stencil_set_factor(sten, 11, FD_D2_O16_A4 * inv2_d1); // i1+4, i2,   i3
				sca_stencil_set_location(sten, 12,  5,  0,  0, 0); sca_stencil_set_factor(sten, 12, FD_D2_O16_A5 * inv2_d1); // i1+5, i2,   i3
				sca_stencil_set_location(sten, 13,  6,  0,  0, 0); sca_stencil_set_factor(sten, 13, FD_D2_O16_A6 * inv2_d1); // i1+6, i2,   i3
				sca_stencil_set_location(sten, 14,  7,  0,  0, 0); sca_stencil_set_factor(sten, 14, FD_D2_O16_A7 * inv2_d1); // i1+7, i2,   i3
				sca_stencil_set_location(sten, 15,  8,  0,  0, 0); sca_stencil_set_factor(sten, 15, FD_D2_O16_A8 * inv2_d1); // i1+8, i2,   i3
				sca_stencil_set_location(sten, 16,  0, -8,  0, 0); sca_stencil_set_factor(sten, 16, FD_D2_O16_A8 * inv2_d2); // i1,   i2-8, i3
				sca_stencil_set_location(sten, 17,  0, -7,  0, 0); sca_stencil_set_factor(sten, 17, FD_D2_O16_A7 * inv2_d2); // i1,   i2-7, i3
				sca_stencil_set_location(sten, 18,  0, -6,  0, 0); sca_stencil_set_factor(sten, 18, FD_D2_O16_A6 * inv2_d2); // i1,   i2-6, i3
				sca_stencil_set_location(sten, 19,  0, -5,  0, 0); sca_stencil_set_factor(sten, 19, FD_D2_O16_A5 * inv2_d2); // i1,   i2-5, i3
				sca_stencil_set_location(sten, 20,  0, -4,  0, 0); sca_stencil_set_factor(sten, 20, FD_D2_O16_A4 * inv2_d2); // i1,   i2-4, i3
				sca_stencil_set_location(sten, 21,  0, -3,  0, 0); sca_stencil_set_factor(sten, 21, FD_D2_O16_A3 * inv2_d2); // i1,   i2-3, i3
				sca_stencil_set_location(sten, 22,  0, -2,  0, 0); sca_stencil_set_factor(sten, 22, FD_D2_O16_A2 * inv2_d2); // i1,   i2-2, i3
				sca_stencil_set_location(sten, 23,  0, -1,  0, 0); sca_stencil_set_factor(sten, 23, FD_D2_O16_A1 * inv2_d2); // i1,   i2-1, i3
				sca_stencil_set_location(sten, 24,  0,  1,  0, 0); sca_stencil_set_factor(sten, 24, FD_D2_O16_A1 * inv2_d2); // i1,   i2+1, i3
				sca_stencil_set_location(sten, 25,  0,  2,  0, 0); sca_stencil_set_factor(sten, 25, FD_D2_O16_A2 * inv2_d2); // i1,   i2+2, i3
				sca_stencil_set_location(sten, 26,  0,  3,  0, 0); sca_stencil_set_factor(sten, 26, FD_D2_O16_A3 * inv2_d2); // i1,   i2+3, i3
				sca_stencil_set_location(sten, 27,  0,  4,  0, 0); sca_stencil_set_factor(sten, 27, FD_D2_O16_A4 * inv2_d2); // i1,   i2+4, i3
				sca_stencil_set_location(sten, 28,  0,  5,  0, 0); sca_stencil_set_factor(sten, 28, FD_D2_O16_A5 * inv2_d2); // i1,   i2+5, i3
				sca_stencil_set_location(sten, 29,  0,  6,  0, 0); sca_stencil_set_factor(sten, 29, FD_D2_O16_A6 * inv2_d2); // i1,   i2+6, i3
				sca_stencil_set_location(sten, 30,  0,  7,  0, 0); sca_stencil_set_factor(sten, 30, FD_D2_O16_A7 * inv2_d2); // i1,   i2+7, i3
				sca_stencil_set_location(sten, 31,  0,  8,  0, 0); sca_stencil_set_factor(sten, 31, FD_D2_O16_A8 * inv2_d2); // i1,   i2+8, i3
				sca_stencil_set_location(sten, 32,  0,  0, -8, 0); sca_stencil_set_factor(sten, 32, FD_D2_O16_A8 * inv2_d3); // i1,   i2,   i3-8
				sca_stencil_set_location(sten, 33,  0,  0, -7, 0); sca_stencil_set_factor(sten, 33, FD_D2_O16_A7 * inv2_d3); // i1,   i2,   i3-7
				sca_stencil_set_location(sten, 34,  0,  0, -6, 0); sca_stencil_set_factor(sten, 34, FD_D2_O16_A6 * inv2_d3); // i1,   i2,   i3-6
				sca_stencil_set_location(sten, 35,  0,  0, -5, 0); sca_stencil_set_factor(sten, 35, FD_D2_O16_A5 * inv2_d3); // i1,   i2,   i3-5
				sca_stencil_set_location(sten, 36,  0,  0, -4, 0); sca_stencil_set_factor(sten, 36, FD_D2_O16_A4 * inv2_d3); // i1,   i2,   i3-4
				sca_stencil_set_location(sten, 37,  0,  0, -3, 0); sca_stencil_set_factor(sten, 37, FD_D2_O16_A3 * inv2_d3); // i1,   i2,   i3-3
				sca_stencil_set_location(sten, 38,  0,  0, -2, 0); sca_stencil_set_factor(sten, 38, FD_D2_O16_A2 * inv2_d3); // i1,   i2,   i3-2
				sca_stencil_set_location(sten, 39,  0,  0, -1, 0); sca_stencil_set_factor(sten, 39, FD_D2_O16_A1 * inv2_d3); // i1,   i2,   i3-1
				sca_stencil_set_location(sten, 40,  0,  0,  1, 0); sca_stencil_set_factor(sten, 40, FD_D2_O16_A1 * inv2_d3); // i1,   i2,   i3+1
				sca_stencil_set_location(sten, 41,  0,  0,  2, 0); sca_stencil_set_factor(sten, 41, FD_D2_O16_A2 * inv2_d3); // i1,   i2,   i3+2
				sca_stencil_set_location(sten, 42,  0,  0,  3, 0); sca_stencil_set_factor(sten, 42, FD_D2_O16_A3 * inv2_d3); // i1,   i2,   i3+3
				sca_stencil_set_location(sten, 43,  0,  0,  4, 0); sca_stencil_set_factor(sten, 43, FD_D2_O16_A4 * inv2_d3); // i1,   i2,   i3+4
				sca_stencil_set_location(sten, 44,  0,  0,  5, 0); sca_stencil_set_factor(sten, 44, FD_D2_O16_A5 * inv2_d3); // i1,   i2,   i3+5
				sca_stencil_set_location(sten, 45,  0,  0,  6, 0); sca_stencil_set_factor(sten, 45, FD_D2_O16_A6 * inv2_d3); // i1,   i2,   i3+6
				sca_stencil_set_location(sten, 46,  0,  0,  7, 0); sca_stencil_set_factor(sten, 46, FD_D2_O16_A7 * inv2_d3); // i1,   i2,   i3+7
				sca_stencil_set_location(sten, 47,  0,  0,  8, 0); sca_stencil_set_factor(sten, 47, FD_D2_O16_A8 * inv2_d3); // i1,   i2,   i3+8
				sca_stencil_set_location(sten, 48,  0,  0,  0, 0); sca_stencil_set_factor(sten, 48, FD_D2_O16_A0 * inv2_d1 + FD_D2_O16_A0 * inv2_d2 + FD_D2_O16_A0 * inv2_d3); // i1, i2, i3
			}
		}

		Myfloat * const w = Wgrid.grid_3d ;
		Myfloat * const u = this->grid_3d ;

		// sca_code_create
		Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
		getGridIndex(pType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End) ;
		for (Myint i = 0; i < nPtPerStencil; i++) sca_stencil_set_input_array(sten, i, 1, n1, n2, n3, &u[i1Start + i2Start*n1 + i3Start*n2*n1]);
		sca_stencil_set_output_array(sten, 1, n1, n2, n3, &w[i1Start + i2Start*n1 + i3Start*n2*n1]);
		sca_code_create(&code_FD_LAPLACIAN, sten, i1End-i1Start+1, i2End-i2Start+1, i3End-i3Start+1, 1);
		sca_stencil_destroy(sten);
	}

	flag_code_FD_LAPLACIAN = true ;
	return(RTN_CODE_OK) ;

#else
	printError("Grid_NEC_SCA::initialize_code_FD_LAPLACIAN, not supported on your platform") ;
	return(RTN_CODE_KO) ;
#endif

	printDebug(MID_DEBUG, "OUT Grid_NEC_SCA::initialize_code_FD_LAPLACIAN");
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_NEC_SCA::applyBoundaryCondition(BoundCond_type boundCondType)
{
	printDebug(FULL_DEBUG, "IN Grid_NEC_SCA::applyBoundaryCondition");

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
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
						// set inner point to 0
						Myint64 iInner1 = i1End+1 ;
						grid_3d[iInner1+i2*n1+i3*n1*n2] = 0.0 ;
#pragma omp simd
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							// set symetrical point to minus inner point value
							grid_3d[i1+i2*n1+i3*n1*n2] = -grid_3d[(iInner1+iInner1-i1)+i2*n1+i3*n1*n2] ;
						}
					}
				}
			}

			// I1HALO2
			if (getNeighbourProc(I1HALO2) == MPI_PROC_NULL)
			{
				getGridIndex(I1HALO2, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End) ;
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i2 = i2Start; i2<= i2End; i2++)
					{
						// set inner point to 0
						Myint64 iInner1 = i1Start-1 ;
						grid_3d[iInner1+i2*n1+i3*n1*n2] = 0.0 ;
#pragma omp simd
						for (Myint64 i1 = i1Start; i1<= i1End; i1++)
						{
							// set symetrical point to minus inner point value
							grid_3d[i1+i2*n1+i3*n1*n2] = -grid_3d[(iInner1-(i1-iInner1))+i2*n1+i3*n1*n2] ;
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
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i1 = i1Start; i1<= i1End; i1++)
					{
						// set inner point to 0
						Myint64 iInner2 = i2End+1 ;
						grid_3d[i1+iInner2*n1+i3*n1*n2] = 0.0 ;
#pragma omp simd
						for (Myint64 i2 = i2Start; i2<= i2End; i2++)
						{
							// set symetrical point to minus inner point value
							grid_3d[i1+i2*n1+i3*n1*n2] = -grid_3d[i1+(iInner2+iInner2-i2)*n1+i3*n1*n2] ;
						}
					}
				}
			}

			// I2HALO2
			if (getNeighbourProc(I2HALO2) == MPI_PROC_NULL)
			{
				getGridIndex(I2HALO2, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End) ;
#pragma omp parallel for collapse(2)
				for (Myint64 i3 = i3Start; i3<= i3End; i3++)
				{
					for (Myint64 i1 = i1Start; i1<= i1End; i1++)
					{
						// set inner point to 0
						Myint64 iInner2 = i2Start-1 ;
						grid_3d[i1+iInner2*n1+i3*n1*n2] = 0.0 ;
#pragma omp simd
						for (Myint64 i2 = i2Start; i2<= i2End; i2++)
						{
							// set symetrical point to minus inner point value
							grid_3d[i1+i2*n1+i3*n1*n2] = -grid_3d[i1+(iInner2-(i2-iInner2))*n1+i3*n1*n2] ;
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
#pragma omp parallel for collapse(2)
				for (Myint64 i2 = i2Start; i2<= i2End; i2++)
				{
					for (Myint64 i1 = i1Start; i1<= i1End; i1++)
					{
						// set inner point to 0
						Myint64 iInner3 = i3End+1 ;
						grid_3d[i1+i2*n1+iInner3*n1*n2] = 0.0 ;
#pragma omp simd
						for (Myint64 i3 = i3Start; i3<= i3End; i3++)
						{
							// set symetrical point to minus inner point value
							grid_3d[i1+i2*n1+i3*n1*n2] = -grid_3d[i1+i2*n1+(iInner3+iInner3-i3)*n1*n2] ;
						}
					}
				}
			}

			// I3HALO2
			if (getNeighbourProc(I3HALO2) == MPI_PROC_NULL)
			{
				getGridIndex(I3HALO2, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End) ;
#pragma omp parallel for collapse(2)
				for (Myint64 i2 = i2Start; i2<= i2End; i2++)
				{
					for (Myint64 i1 = i1Start; i1<= i1End; i1++)
					{
						// set inner point to 0
						Myint64 iInner3 = i3Start-1 ;
						grid_3d[i1+i2*n1+iInner3*n1*n2] = 0.0 ;
#pragma omp simd
						for (Myint64 i3 = i3Start; i3<= i3End; i3++)
						{
							// set symetrical point to minus inner point value
							grid_3d[i1+i2*n1+i3*n1*n2] = -grid_3d[i1+i2*n1+(iInner3-(i3-iInner3))*n1*n2] ;
						}
					}
				}
			}
		}
	}
	else
	{
		printError("IN Grid_NEC_SCA::applyBoundaryCondition, invalid boundCondType", boundCondType) ;
		return(RTN_CODE_KO) ;
	}

	printDebug(FULL_DEBUG, "OUT Grid_NEC_SCA::applyBoundaryCondition");
	return(RTN_CODE_OK) ;
}

} // namespace hpcscan
