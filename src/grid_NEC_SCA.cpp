
//-------------------------------------------------------------------------------------------------------
// This grid is activated with command line option -testMode NEC_SCA
// Derived class from Grid_NEC
// NEC Stencil Code Accelerator (target NEC SX-Aurora TSUBASA Vector Engine)
//-------------------------------------------------------------------------------------------------------

#include "grid_NEC_SCA.h"

#include <cassert>
#include <cfloat>  // for FLT_MAX ;
#include <cmath>   // for fabs
#include <cstddef> // for NULL
#include <fstream>
#include <stdio.h>

#include "mpi.h"
#include <sca.h>

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
#define sca_stencil_append_elements_xyza_ordered sca_stencil_append_elements_xyza_ordered_d
#define sca_stencil_append_elements_xa_ordered sca_stencil_append_elements_xa_ordered_d
#define sca_stencil_append_elements_ya_ordered sca_stencil_append_elements_ya_ordered_d
#define sca_stencil_append_elements_za_ordered sca_stencil_append_elements_za_ordered_d
#define sca_stencil_append_elements_xya_ordered sca_stencil_append_elements_xya_ordered_d
#else
#define sca_stencil_create sca_stencil_create_s
#define sca_stencil_set_factor sca_stencil_set_factor_s
#define sca_stencil_set_input_array sca_stencil_set_input_array_s
#define sca_stencil_set_output_array sca_stencil_set_output_array_s
#define sca_stencil_append_elements_xyza_ordered sca_stencil_append_elements_xyza_ordered_s
#define sca_stencil_append_elements_xa_ordered sca_stencil_append_elements_xa_ordered_s
#define sca_stencil_append_elements_ya_ordered sca_stencil_append_elements_ya_ordered_s
#define sca_stencil_append_elements_za_ordered sca_stencil_append_elements_za_ordered_s
#define sca_stencil_append_elements_xya_ordered sca_stencil_append_elements_xya_ordered_s
#endif

using namespace std;

namespace hpcscan {

const float FD_O2_1D[17][17] = {
		{},
		{},
		{FD_D2_O2_A1, FD_D2_O2_A0, FD_D2_O2_A1}, // fdOrder 2
		{},
		{FD_D2_O4_A2, FD_D2_O4_A1, FD_D2_O4_A0, FD_D2_O4_A1, FD_D2_O4_A2}, // fdOrder 4
		{},
		{FD_D2_O6_A3, FD_D2_O6_A2, FD_D2_O6_A1, FD_D2_O6_A0, FD_D2_O6_A1, FD_D2_O6_A2, FD_D2_O6_A3}, // fdOrder 6
		{},
		{FD_D2_O8_A4, FD_D2_O8_A3, FD_D2_O8_A2, FD_D2_O8_A1, FD_D2_O8_A0, FD_D2_O8_A1, FD_D2_O8_A2, FD_D2_O8_A3, FD_D2_O8_A4}, // fdOrder 8
		{},
		{FD_D2_O10_A5, FD_D2_O10_A4, FD_D2_O10_A3, FD_D2_O10_A2, FD_D2_O10_A1, FD_D2_O10_A0, FD_D2_O10_A1, FD_D2_O10_A2, FD_D2_O10_A3, FD_D2_O10_A4, FD_D2_O10_A5}, // fdOrder 10
		{},
		{FD_D2_O12_A6, FD_D2_O12_A5, FD_D2_O12_A4, FD_D2_O12_A3, FD_D2_O12_A2, FD_D2_O12_A1, FD_D2_O12_A0, FD_D2_O12_A1, FD_D2_O12_A2, FD_D2_O12_A3, FD_D2_O12_A4, FD_D2_O12_A5, FD_D2_O12_A6}, // fdOrder 12
		{},
		{FD_D2_O14_A7, FD_D2_O14_A6, FD_D2_O14_A5, FD_D2_O14_A4, FD_D2_O14_A3, FD_D2_O14_A2, FD_D2_O14_A1, FD_D2_O14_A0, FD_D2_O14_A1, FD_D2_O14_A2, FD_D2_O14_A3, FD_D2_O14_A4, FD_D2_O14_A5, FD_D2_O14_A6, FD_D2_O14_A7}, // fdOrder 14
		{},
		{FD_D2_O16_A8, FD_D2_O16_A7, FD_D2_O16_A6, FD_D2_O16_A5, FD_D2_O16_A4, FD_D2_O16_A3, FD_D2_O16_A2, FD_D2_O16_A1, FD_D2_O16_A0, FD_D2_O16_A1, FD_D2_O16_A2, FD_D2_O16_A3, FD_D2_O16_A4, FD_D2_O16_A5, FD_D2_O16_A6, FD_D2_O16_A7, FD_D2_O16_A8}, // fdOrder 16
};

//-------------------------------------------------------------------------------------------------------

Grid_NEC_SCA::Grid_NEC_SCA(Grid_type gridTypeIn) : Grid_NEC(gridTypeIn)
{
	printDebug(MID_DEBUG, "IN Grid_NEC_SCA::Grid_NEC_SCA");

	gridMode = GRID_MODE_NEC_SCA ;

	flag_code_FD_D2_N1  = false ;
	flag_code_FD_D2_N2  = false ;
	flag_code_FD_D2_N3  = false ;
	flag_code_FD_LAPLACIAN = false ;

	printDebug(MID_DEBUG, "OUT Grid_NEC_SCA::Grid_NEC_SCA");
}

//-------------------------------------------------------------------------------------------------------

Grid_NEC_SCA::Grid_NEC_SCA(Grid_type gridTypeIn, Dim_type dimIn,
		Myint64 n1InnerIn, Myint64 n2InnerIn, Myint64 n3InnerIn) : Grid_NEC(gridTypeIn, dimIn,
				n1InnerIn, n2InnerIn, n3InnerIn)
{
	printDebug(MID_DEBUG, "IN Grid_NEC_SCA::Grid_NEC_SCA");

	gridMode = GRID_MODE_NEC_SCA ;

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
	if (flag_code_FD_D2_N1) sca_code_destroy(code_FD_D2_N1);
	if (flag_code_FD_D2_N2) sca_code_destroy(code_FD_D2_N2);
	if (flag_code_FD_D2_N3) sca_code_destroy(code_FD_D2_N3);
	if (flag_code_FD_LAPLACIAN) sca_code_destroy(code_FD_LAPLACIAN);

	printDebug(MID_DEBUG, "OUT Grid_NEC_SCA::~Grid_NEC_SCA");
}

//-------------------------------------------------------------------------------------------------------

void Grid_NEC_SCA::info(void)
{
	printDebug(FULL_DEBUG, "IN Grid_NEC_SCA::info");

	// parent class info
	Grid_NEC::info() ;

	// additional info
	//printInfo(MASTER, " TO BE COMPLETED") ;

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
	sca_code_execute(code_FD_D2_N1);

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
	sca_code_execute(code_FD_D2_N2);

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
	sca_code_execute(code_FD_D2_N3);

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
	sca_code_execute(code_FD_LAPLACIAN);

	printDebug(MID_DEBUG, "OUT Grid_NEC_SCA::FD_LAPLACIAN");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_NEC_SCA::initializeGrid(void)
{
	printDebug(MID_DEBUG, "IN Grid_NEC_SCA::initializeGrid");

	Grid_NEC::initializeGrid() ;

	printDebug(MID_DEBUG, "OUT Grid_NEC_SCA::initializeGrid");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_NEC_SCA::initialize_code_FD_D2_N1(Point_type pType, const Grid& Wgrid, Myint fdOrder)
{
	printDebug(MID_DEBUG, "IN Grid_NEC_SCA::initialize_code_FD_D2_N1");

	if (!flag_code_FD_D2_N1)
	{
		Myint nPtPerStencil = getPtPerStencilFD_D2(fdOrder) ;
		Myint fdOrderHalf = fdOrder / 2;

		const Myfloat inv_d1  = Myfloat(1.0) / d1 ;
		const Myfloat inv_d2  = Myfloat(1.0) / d2 ;
		const Myfloat inv_d3  = Myfloat(1.0) / d3 ;

		const Myfloat inv2_d1 = inv_d1 * inv_d1 ;
		const Myfloat inv2_d2 = inv_d2 * inv_d2 ;
		const Myfloat inv2_d3 = inv_d3 * inv_d3 ;

		// check if fdOrder is supported
		Myint64 nFdOrder = sizeof(FD_O2_1D) / sizeof(FD_O2_1D[0]) ;
		if (fdOrder >= nFdOrder) {
			printError("Grid_NEC_SCA::initialize_code_FD_D2_N1, fdOrder not supported") ;
			return(RTN_CODE_KO) ;
		}

		Myfloat coefs[fdOrder+1];
		for (Myint i = 0; i < fdOrder + 1; i++) coefs[i] = FD_O2_1D[fdOrder][i] * inv2_d1;

		Myfloat * const w = Wgrid.grid_3d ;
		Myfloat * const u = this->grid_3d ;

		Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
		getGridIndex(pType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End) ;

		// sca_code_create
		sca_stencil_t sten ;
		sca_stencil_create(&sten);
		sca_stencil_append_elements_xa_ordered(sten, fdOrderHalf, fdOrderHalf, 1, n1, n2, 0, &u[i1Start + i2Start*n1 + i3Start*n2*n1], coefs);
		sca_stencil_set_output_array(sten, 1, n1, n2, 0, &w[i1Start + i2Start*n1 + i3Start*n2*n1]);
		sca_code_create(&code_FD_D2_N1, sten, i1End-i1Start+1, i2End-i2Start+1, i3End-i3Start+1, 1);
		sca_stencil_destroy(sten);
	}

	flag_code_FD_D2_N1 = true ;

	printDebug(MID_DEBUG, "OUT Grid_NEC_SCA::initialize_code_FD_D2_N1");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_NEC_SCA::initialize_code_FD_D2_N2(Point_type pType, const Grid& Wgrid, Myint fdOrder)
{
	printDebug(MID_DEBUG, "IN Grid_NEC_SCA::initialize_code_FD_D2_N2");

	if (!flag_code_FD_D2_N2)
	{
		Myint nPtPerStencil = getPtPerStencilFD_D2(fdOrder) ;
		Myint fdOrderHalf = fdOrder / 2;

		const Myfloat inv_d1  = Myfloat(1.0) / d1 ;
		const Myfloat inv_d2  = Myfloat(1.0) / d2 ;
		const Myfloat inv_d3  = Myfloat(1.0) / d3 ;

		const Myfloat inv2_d1 = inv_d1 * inv_d1 ;
		const Myfloat inv2_d2 = inv_d2 * inv_d2 ;
		const Myfloat inv2_d3 = inv_d3 * inv_d3 ;

		// check if fdOrder is supported
		Myint64 nFdOrder = sizeof(FD_O2_1D) / sizeof(FD_O2_1D[0]) ;
		if (fdOrder >= nFdOrder) {
			printError("Grid_NEC_SCA::initialize_code_FD_D2_N2, fdOrder not supported") ;
			return(RTN_CODE_KO) ;
		}

		Myfloat coefs[fdOrder+1];
		for (Myint i = 0; i < fdOrder + 1; i++) coefs[i] = FD_O2_1D[fdOrder][i] * inv2_d2;

		Myfloat * const w = Wgrid.grid_3d ;
		Myfloat * const u = this->grid_3d ;

		Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
		getGridIndex(pType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End) ;

		// sca_code_create
		sca_stencil_t sten ;
		sca_stencil_create(&sten);
		sca_stencil_append_elements_ya_ordered(sten, fdOrderHalf, fdOrderHalf, 1, n1, n2, 0, &u[i1Start + i2Start*n1 + i3Start*n2*n1], coefs);
		sca_stencil_set_output_array(sten, 1, n1, n2, 0, &w[i1Start + i2Start*n1 + i3Start*n2*n1]);
		sca_code_create(&code_FD_D2_N2, sten, i1End-i1Start+1, i2End-i2Start+1, i3End-i3Start+1, 1);
		sca_stencil_destroy(sten);
	}

	flag_code_FD_D2_N2 = true ;

	printDebug(MID_DEBUG, "OUT Grid_NEC_SCA::initialize_code_FD_D2_N2");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_NEC_SCA::initialize_code_FD_D2_N3(Point_type pType, const Grid& Wgrid, Myint fdOrder)
{
	printDebug(MID_DEBUG, "IN Grid_NEC_SCA::initialize_code_FD_D2_N3");

	if (!flag_code_FD_D2_N3)
	{
		Myint nPtPerStencil = getPtPerStencilFD_D2(fdOrder) ;
		Myint fdOrderHalf = fdOrder / 2;

		const Myfloat inv_d1  = Myfloat(1.0) / d1 ;
		const Myfloat inv_d2  = Myfloat(1.0) / d2 ;
		const Myfloat inv_d3  = Myfloat(1.0) / d3 ;

		const Myfloat inv2_d1 = inv_d1 * inv_d1 ;
		const Myfloat inv2_d2 = inv_d2 * inv_d2 ;
		const Myfloat inv2_d3 = inv_d3 * inv_d3 ;

		// check if fdOrder is supported
		Myint64 nFdOrder = sizeof(FD_O2_1D) / sizeof(FD_O2_1D[0]) ;
		if (fdOrder >= nFdOrder) {
			printError("Grid_NEC_SCA::initialize_code_FD_D2_N3, fdOrder not supported") ;
			return(RTN_CODE_KO) ;
		}

		Myfloat coefs[fdOrder+1];
		for (Myint i = 0; i < fdOrder + 1; i++) coefs[i] = FD_O2_1D[fdOrder][i] * inv2_d3;

		Myfloat * const w = Wgrid.grid_3d ;
		Myfloat * const u = this->grid_3d ;

		Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
		getGridIndex(pType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End) ;

		// sca_code_create
		sca_stencil_t sten ;
		sca_stencil_create(&sten);
		sca_stencil_append_elements_za_ordered(sten, fdOrderHalf, fdOrderHalf, 1, n1, n2, 0, &u[i1Start + i2Start*n1 + i3Start*n2*n1], coefs);
		sca_stencil_set_output_array(sten, 1, n1, n2, 0, &w[i1Start + i2Start*n1 + i3Start*n2*n1]);
		sca_code_create(&code_FD_D2_N3, sten, i1End-i1Start+1, i2End-i2Start+1, i3End-i3Start+1, 1);
		sca_stencil_destroy(sten);
	}

	flag_code_FD_D2_N3 = true ;

	printDebug(MID_DEBUG, "OUT Grid_NEC_SCA::initialize_code_FD_D2_N3");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid_NEC_SCA::initialize_code_FD_LAPLACIAN(Point_type pType, const Grid& Wgrid, Myint fdOrder)
{
	printDebug(MID_DEBUG, "IN Grid_NEC_SCA::initialize_code_FD_LAPLACIAN");

	if (!flag_code_FD_LAPLACIAN)
	{
		Myint nPtPerStencil = getPtPerStencilFD_LAPLACIAN(fdOrder) ;

		const Myfloat inv_d1  = Myfloat(1.0) / d1 ;
		const Myfloat inv_d2  = Myfloat(1.0) / d2 ;
		const Myfloat inv_d3  = Myfloat(1.0) / d3 ;

		const Myfloat inv2_d1 = inv_d1 * inv_d1 ;
		const Myfloat inv2_d2 = inv_d2 * inv_d2 ;
		const Myfloat inv2_d3 = inv_d3 * inv_d3 ;

		// check if fdOrder is supported
		Myint64 nFdOrder = sizeof(FD_O2_1D) / sizeof(FD_O2_1D[0]) ;
		if (fdOrder >= nFdOrder) {
			printError("Grid_NEC_SCA::initialize_code_FD_LAPLACIAN, fdOrder not supported") ;
			return(RTN_CODE_KO) ;
		}

		Myfloat * const w = Wgrid.grid_3d ;
		Myfloat * const u = this->grid_3d ;

		Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
		getGridIndex(pType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End) ;

		Myint fdOrderHalf = fdOrder / 2;

		if (dim == DIM1)
		{
			Myfloat coefs[fdOrder+1];

			for (Myint i = 0; i < fdOrder + 1; i++) coefs[i] = FD_O2_1D[fdOrder][i] * inv2_d1;

			// sca_code_create
			sca_stencil_t sten ;
			sca_stencil_create(&sten);
			sca_stencil_append_elements_xa_ordered(sten, fdOrderHalf, fdOrderHalf, 1, n1, n2, 0, &u[i1Start + i2Start*n1 + i3Start*n2*n1], coefs);
			sca_stencil_set_output_array(sten, 1, n1, n2, 0, &w[i1Start + i2Start*n1 + i3Start*n2*n1]);
			sca_code_create(&code_FD_LAPLACIAN, sten, i1End-i1Start+1, i2End-i2Start+1, i3End-i3Start+1, 1);
			sca_stencil_destroy(sten);
		}

		else if (dim == DIM2)
		{
			Myint n = 0;
			Myfloat coefs[2*fdOrder+1];

			for (Myint i = 0; i < fdOrderHalf; i++) { coefs[n] = FD_O2_1D[fdOrder][i] * inv2_d2; n++; }  // i2-n,..,i2-1
			for (Myint i = 0; i < fdOrderHalf; i++) { coefs[n] = FD_O2_1D[fdOrder][i] * inv2_d1; n++; }  // i1-n,..,i1-1
			coefs[n] = FD_O2_1D[fdOrder][fdOrderHalf] * inv2_d1 + FD_O2_1D[fdOrder][fdOrderHalf] * inv2_d2; n++; // i1,i2
			for (Myint i = 0; i < fdOrderHalf; i++) { coefs[n] = FD_O2_1D[fdOrder][i+fdOrderHalf+1] * inv2_d1; n++; }  // i1+1,..,i1+n
			for (Myint i = 0; i < fdOrderHalf; i++) { coefs[n] = FD_O2_1D[fdOrder][i+fdOrderHalf+1] * inv2_d2; n++; }  // i2+1,..,i2+n

			// sca_code_create
			sca_stencil_t sten ;
			sca_stencil_create(&sten);
			sca_stencil_append_elements_xya_ordered(sten, fdOrderHalf, fdOrderHalf, fdOrderHalf, fdOrderHalf, 1, n1, n2, 0, &u[i1Start + i2Start*n1 + i3Start*n2*n1], coefs);
			sca_stencil_set_output_array(sten, 1, n1, n2, 0, &w[i1Start + i2Start*n1 + i3Start*n2*n1]);
			sca_code_create(&code_FD_LAPLACIAN, sten, i1End-i1Start+1, i2End-i2Start+1, i3End-i3Start+1, 1);
			sca_stencil_destroy(sten);
		}

		else if (dim == DIM3)
		{
			Myint n = 0;
			Myfloat coefs[3*fdOrder+1];

			for (Myint i = 0; i < fdOrderHalf; i++) { coefs[n] = FD_O2_1D[fdOrder][i] * inv2_d3; n++; }  // i3-n,..,i3-1
			for (Myint i = 0; i < fdOrderHalf; i++) { coefs[n] = FD_O2_1D[fdOrder][i] * inv2_d2; n++; }  // i2-n,..,i2-1
			for (Myint i = 0; i < fdOrderHalf; i++) { coefs[n] = FD_O2_1D[fdOrder][i] * inv2_d1; n++; }  // i1-n,..,i1-1
			coefs[n] = FD_O2_1D[fdOrder][fdOrderHalf] * inv2_d1 + FD_O2_1D[fdOrder][fdOrderHalf] * inv2_d2 + FD_O2_1D[fdOrder][fdOrderHalf] * inv2_d3; n++; // i1,i2,i3
			for (Myint i = 0; i < fdOrderHalf; i++) { coefs[n] = FD_O2_1D[fdOrder][i+fdOrderHalf+1] * inv2_d1; n++; }  // i1+1,..,i1+n
			for (Myint i = 0; i < fdOrderHalf; i++) { coefs[n] = FD_O2_1D[fdOrder][i+fdOrderHalf+1] * inv2_d2; n++; }  // i2+1,..,i2+n
			for (Myint i = 0; i < fdOrderHalf; i++) { coefs[n] = FD_O2_1D[fdOrder][i+fdOrderHalf+1] * inv2_d3; n++; }  // i3+1,..,i3+n

			// sca_code_create
			sca_stencil_t sten ;
			sca_stencil_create(&sten);
			sca_stencil_append_elements_xyza_ordered(sten, fdOrderHalf, fdOrderHalf, fdOrderHalf, fdOrderHalf,
					fdOrderHalf, fdOrderHalf, 1, n1, n2, 0, &u[i1Start + i2Start*n1 + i3Start*n2*n1], coefs);
			sca_stencil_set_output_array(sten, 1, n1, n2, 0, &w[i1Start + i2Start*n1 + i3Start*n2*n1]);
			sca_code_create(&code_FD_LAPLACIAN, sten, i1End-i1Start+1, i2End-i2Start+1, i3End-i3Start+1, 1);
			sca_stencil_destroy(sten);
		}
	}

	flag_code_FD_LAPLACIAN = true ;

	printDebug(MID_DEBUG, "OUT Grid_NEC_SCA::initialize_code_FD_LAPLACIAN");
	return(RTN_CODE_OK) ;
}

} // namespace hpcscan
