
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
	//Grid::initializeGrid() ;

	Myint nlayer = Config::Instance()->nlayer ;

	if (Config::Instance()->fdOrder == 2)
	{
		haloWidth = FD_D2_O2_HSTL ;
	}
	else if (Config::Instance()->fdOrder == 4)
	{
		haloWidth = FD_D2_O4_HSTL ;
	}
	else if (Config::Instance()->fdOrder == 8)
	{
		haloWidth = FD_D2_O8_HSTL ;
	}
	else if (Config::Instance()->fdOrder == 12)
	{
		haloWidth = FD_D2_O12_HSTL ;
	}
	else if (Config::Instance()->fdOrder == 16)
	{
		haloWidth = FD_D2_O16_HSTL ;
	}
	else
	{
		haloWidth = 0 ;
	}

	// get global grid size (Inner points)
	Myint n1InnerLoc = n1Inner ;
	Myint n2InnerLoc = n2Inner ;
	Myint n3InnerLoc = n3Inner ;

	Myint subIdx1 = 0 ;
	Myint subIdx2 = 0 ;
	Myint subIdx3 = 0 ;

	i1ProcIdStart = MPI_PROC_NULL ;
	i1ProcIdEnd   = MPI_PROC_NULL ;
	i2ProcIdStart = MPI_PROC_NULL ;
	i2ProcIdEnd   = MPI_PROC_NULL ;
	i3ProcIdStart = MPI_PROC_NULL ;
	i3ProcIdEnd   = MPI_PROC_NULL ;

	// divide grid in subdomains
	//--------------------------
	if (gridType == GRID_LOCAL)
	{
		// retrieve index of subdomain in global grid
		//-------------------------------------------
		subIdx3 = myid_world / (Config::Instance()->nsub1*Config::Instance()->nsub2) ;
		subIdx2 = (myid_world - subIdx3 *
				(Config::Instance()->nsub1*Config::Instance()->nsub2)) / Config::Instance()->nsub1 ;
		subIdx1 = myid_world
				- subIdx3 * (Config::Instance()->nsub1*Config::Instance()->nsub2)
				- subIdx2 * (Config::Instance()->nsub1) ;

		printDebug(LIGHT_DEBUG, "subIdx1", subIdx1) ;
		printDebug(LIGHT_DEBUG, "subIdx2", subIdx2) ;
		printDebug(LIGHT_DEBUG, "subIdx3", subIdx3) ;

		// retrieve number of points in inner domain
		//------------------------------------------
		Myint n1InnerLocTmp = n1InnerLoc / Config::Instance()->nsub1 ;
		if (subIdx1 < (Config::Instance()->nsub1-1))
		{
			n1InnerLoc = n1InnerLocTmp ;
		}
		else
		{
			n1InnerLoc = n1InnerLoc - (Config::Instance()->nsub1-1)*n1InnerLocTmp ;
		}

		Myint n2InnerLocTmp = n2InnerLoc / Config::Instance()->nsub2 ;
		if (subIdx2 < (Config::Instance()->nsub2-1))
		{
			n2InnerLoc = n2InnerLocTmp ;
		}
		else
		{
			n2InnerLoc = n2InnerLoc - (Config::Instance()->nsub2-1)*n2InnerLocTmp ;
		}

		Myint n3InnerLocTmp = n3InnerLoc / Config::Instance()->nsub3 ;
		if (subIdx3 < (Config::Instance()->nsub3-1))
		{
			n3InnerLoc = n3InnerLocTmp ;
		}
		else
		{
			n3InnerLoc = n3InnerLoc - (Config::Instance()->nsub3-1)*n3InnerLocTmp ;
		}

		// retrieve MPI rank of neighbours
		//--------------------------------
		if (Config::Instance()->nsub1 == 1)
		{
			i1ProcIdStart = MPI_PROC_NULL ;
			i1ProcIdEnd   = MPI_PROC_NULL ;
		}
		else
		{
			if (subIdx1 == 0)
			{
				i1ProcIdStart = MPI_PROC_NULL ;
			}
			else
			{
				i1ProcIdStart = myid_world - 1 ;
			}

			if (subIdx1 == Config::Instance()->nsub1-1)
			{
				i1ProcIdEnd = MPI_PROC_NULL ;
			}
			else
			{
				i1ProcIdEnd = myid_world + 1 ;
			}
		}

		if (Config::Instance()->nsub2 == 1)
		{
			i2ProcIdStart = MPI_PROC_NULL ;
			i2ProcIdEnd   = MPI_PROC_NULL ;
		}
		else
		{
			if (subIdx2 == 0)
			{
				i2ProcIdStart = MPI_PROC_NULL ;
			}
			else
			{
				i2ProcIdStart = myid_world - Config::Instance()->nsub1 ;
			}

			if (subIdx2 == Config::Instance()->nsub2-1)
			{
				i2ProcIdEnd = MPI_PROC_NULL ;
			}
			else
			{
				i2ProcIdEnd = myid_world + Config::Instance()->nsub1 ;
			}
		}

		if (Config::Instance()->nsub3 == 1)
		{
			i3ProcIdStart = MPI_PROC_NULL ;
			i3ProcIdEnd   = MPI_PROC_NULL ;
		}
		else
		{
			if (subIdx3 == 0)
			{
				i3ProcIdStart = MPI_PROC_NULL ;
			}
			else
			{
				i3ProcIdStart = myid_world -
						(Config::Instance()->nsub1 * Config::Instance()->nsub2);
			}

			if (subIdx3 == Config::Instance()->nsub3-1)
			{
				i3ProcIdEnd = MPI_PROC_NULL ;
			}
			else
			{
				i3ProcIdEnd = myid_world +
						(Config::Instance()->nsub1 * Config::Instance()->nsub2);
			}
		}
	}

	//================================================================================
	// initialize grid index axis 1
	//================================================================================
	i1Halo1Start = 0 ;
	i1Halo1End   = i1Halo1Start + haloWidth - 1 ;
	if (nlayer == 0)
	{
		i1Layer1Start = i1Halo1End + 1 ;
		i1Layer1End   = i1Layer1Start ;
		i1InnerStart  = i1Layer1Start ;
	}
	else
	{
		i1Layer1Start = i1Halo1End + 1 ;
		i1Layer1End   = i1Layer1Start + nlayer - 1 ;
		i1InnerStart  = i1Layer1End + 1 ;
	}
	i1InnerEnd   = i1InnerStart + n1InnerLoc - 1 ;
	if (nlayer == 0)
	{
		i1Layer2Start = i1InnerEnd ;
		i1Layer2End   = i1InnerEnd ;
	}
	else
	{
		i1Layer2Start = i1InnerEnd + 1 ;
		i1Layer2End   = i1Layer2Start + nlayer - 1 ;
	}
	i1Halo2Start = i1Layer2End + 1 ;
	i1Halo2End   = i1Halo2Start + haloWidth - 1  ;

	// Padding
	padGridn1() ;

	// total n1
	n1 = i1PadEnd + 1 ;

	// spatial sampling
	d1 = Config::Instance()->h ;

	// global coordinate of the grid
	i1OffsetGlob = subIdx1 * (Config::Instance()->n1 / Config::Instance()->nsub1) ;
	printDebug(LIGHT_DEBUG, "i1OffsetGlob", i1OffsetGlob) ;
	Orig1 = (i1OffsetGlob - i1InnerStart) * d1 ;

	//================================================================================
	// initialize grid index axis 2
	//================================================================================
	if (dim >= DIM2)
	{
		i2Halo1Start = 0 ;
		i2Halo1End   = i2Halo1Start + haloWidth - 1 ;
		if (nlayer == 0)
		{
			i2Layer1Start = i2Halo1End + 1 ;
			i2Layer1End   = i2Layer1Start ;
			i2InnerStart  = i2Layer1Start ;
		}
		else
		{
			i2Layer1Start = i2Halo1End + 1 ;
			i2Layer1End   = i2Layer1Start + nlayer - 1 ;
			i2InnerStart  = i2Layer1End + 1 ;
		}
		i2InnerEnd   = i2InnerStart + n2InnerLoc - 1 ;
		if (nlayer == 0)
		{
			i2Layer2Start = i2InnerEnd ;
			i2Layer2End   = i2InnerEnd ;
		}
		else
		{
			i2Layer2Start = i2InnerEnd + 1 ;
			i2Layer2End   = i2Layer2Start + nlayer - 1 ;
		}
		i2Halo2Start = i2Layer2End + 1 ;
		i2Halo2End   = i2Halo2Start + haloWidth - 1  ;

		// Padding
		padGridn2() ;

		// total n2
		n2 = i2PadEnd + 1 ;

		// spatial sampling
		d2 = Config::Instance()->h ;

		// global coordinate of the grid
		i2OffsetGlob = subIdx2 * (Config::Instance()->n2 / Config::Instance()->nsub2) ;
		printDebug(LIGHT_DEBUG, "i2OffsetGlob", i2OffsetGlob) ;
		Orig2 = (i2OffsetGlob - i2InnerStart) * d2 ;
	}
	else
	{
		i2Halo1Start  = 0 ;
		i2Halo1End    = 0 ;
		i2Layer1Start = 0 ;
		i2Layer1End   = 0 ;
		i2InnerStart  = 0 ;
		i2InnerEnd    = 0 ;
		i2Layer2Start = 0 ;
		i2Layer2End   = 0 ;
		i2Halo2Start  = 0 ;
		i2Halo2End    = 0 ;
		i2PadStart    = 0 ;
		i2PadEnd      = 0 ;
		i2ProcIdStart = MPI_PROC_NULL ;
		i2ProcIdEnd   = MPI_PROC_NULL ;
		d2            = 0 ;
		i2OffsetGlob  = 0 ;
		Orig2         = 0 ;
		n2            = 1 ;
		n2InnerLoc    = 1 ;
	}

	//================================================================================
	// initialize grid index axis 3
	//================================================================================
	if (dim >= DIM3)
	{
		i3Halo1Start = 0 ;
		i3Halo1End   = i3Halo1Start + haloWidth - 1 ;
		if (nlayer == 0)
		{
			i3Layer1Start = i3Halo1End + 1 ;
			i3Layer1End   = i3Layer1Start ;
			i3InnerStart  = i3Layer1Start ;
		}
		else
		{
			i3Layer1Start = i3Halo1End + 1 ;
			i3Layer1End   = i3Layer1Start + nlayer - 1 ;
			i3InnerStart  = i3Layer1End + 1 ;
		}
		i3InnerEnd   = i3InnerStart + n3InnerLoc - 1 ;
		if (nlayer == 0)
		{
			i3Layer2Start = i3InnerEnd ;
			i3Layer2End   = i3InnerEnd ;
		}
		else
		{
			i3Layer2Start = i3InnerEnd + 1 ;
			i3Layer2End   = i3Layer2Start + nlayer - 1 ;
		}
		i3Halo2Start = i3Layer2End + 1 ;
		i3Halo2End   = i3Halo2Start + haloWidth - 1  ;

		// Padding
		padGridn3() ;

		// total n3
		n3 = i3PadEnd + 1 ;

		// spatial sampling
		d3 = Config::Instance()->h ;

		// global coordinate of the grid
		i3OffsetGlob = subIdx3 * (Config::Instance()->n3 / Config::Instance()->nsub3) ;
		printDebug(LIGHT_DEBUG, "i3OffsetGlob", i3OffsetGlob) ;
		Orig3 = (i3OffsetGlob - i3InnerStart) * d3 ;
	}
	else
	{
		i3Halo1Start  = 0 ;
		i3Halo1End    = 0 ;
		i3Layer1Start = 0 ;
		i3Layer1End   = 0 ;
		i3InnerStart  = 0 ;
		i3InnerEnd    = 0 ;
		i3Layer2Start = 0 ;
		i3Layer2End   = 0 ;
		i3Halo2Start  = 0 ;
		i3Halo2End    = 0 ;
		i3PadStart    = 0 ;
		i3PadEnd      = 0 ;
		i3ProcIdStart = MPI_PROC_NULL ;
		i3ProcIdEnd   = MPI_PROC_NULL ;
		d3            = 0 ;
		i3OffsetGlob  = 0 ;
		Orig3         = 0 ;

		n3            = 1 ;
		n3InnerLoc    = 1 ;
	}

	//-------------------------------------------------------------------------------------
	// create MPI communication type using:
	//
	// MPI_Type_contiguous - Creates a contiguous datatype.
	// int MPI_Type_contiguous(int count, MPI_Datatype oldtype,
	//     MPI_Datatype *newtype)
	//
	// MPI_Type_vector - Creates a vector (strided) datatype.
	// int MPI_Type_vector(int count, int blocklength, int stride,
	//     MPI_Datatype oldtype, MPI_Datatype *newtype)
	//
	// MPI_Type_create_hvector - Creates a vector (strided) data type with offset in bytes.
	// int MPI_Type_create_hvector(int count, int blocklength,
	//     MPI_Aint stride, MPI_Datatype oldtype, MPI_Datatype *newtype)
	//
	// MPI_Type_commit - Commits a data type.
	// int MPI_Type_commit(MPI_Datatype *datatype)
	//-------------------------------------------------------------------------------------

	{
		// create i1HaloDataType used for I1HALO1 & I1HALO2
		MPI_Datatype dataType1 ;
		MPI_Type_contiguous(haloWidth, MPI_MYFLOAT, &dataType1) ;

		MPI_Datatype dataType2 ;
		MPI_Aint stride2 = n1*sizeof(Myfloat);
		MPI_Type_create_hvector(n2InnerLoc, 1, stride2, dataType1, &dataType2) ;

		MPI_Aint stride3 = n1*n2*sizeof(Myfloat);
		MPI_Type_create_hvector(n3InnerLoc, 1, stride3, dataType2, &i1HaloDataType) ;
		MPI_Type_commit(&i1HaloDataType) ;
	}

	{
		// create i2HaloDataType used for I2HALO1 & I2HALO2
		MPI_Datatype dataType1 ;
		MPI_Type_contiguous(n1InnerLoc, MPI_MYFLOAT, &dataType1) ;

		MPI_Datatype dataType2 ;
		MPI_Aint stride2 = n1*sizeof(Myfloat);
		MPI_Type_create_hvector(haloWidth, 1, stride2, dataType1, &dataType2) ;

		MPI_Aint stride3 = n1*n2*sizeof(Myfloat);
		MPI_Type_create_hvector(n3InnerLoc, 1, stride3, dataType2, &i2HaloDataType) ;
		MPI_Type_commit(&i2HaloDataType) ;
	}

	{
		// create i3HaloDataType used for I3HALO1 & I3HALO2
		MPI_Datatype dataType1 ;
		MPI_Type_contiguous(n1InnerLoc, MPI_MYFLOAT, &dataType1) ;

		MPI_Datatype dataType2 ;
		MPI_Aint stride2 = n1*sizeof(Myfloat);
		MPI_Type_create_hvector(n2InnerLoc, 1, stride2, dataType1, &dataType2) ;

		MPI_Aint stride3 = n1*n2*sizeof(Myfloat);
		MPI_Type_create_hvector(haloWidth, 1, stride3, dataType2, &i3HaloDataType) ;
		MPI_Type_commit(&i3HaloDataType) ;
	}

	// allocate grid
	//--------------
	npoint = n1 * n2 * n3 ;
	//grid_3d = new Myfloat[npoint] ;
#pragma acc enter data create(grid_3d[0:(n1*n2*n3)-1])

	printDebug(FULL_DEBUG, "Out Grid_GPU2::initializeGrid") ;
}

//-------------------------------------------------------------------------------------------------------
void Grid_GPU2::fill(Point_type pointType, Myfloat val)
{
	printDebug(FULL_DEBUG, "In Grid_GPU2::fill") ;

	Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
	getGridIndex(pointType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End) ;
	//Myfloat *grid_3d_2;
#pragma acc enter data create(grid_3d[0:(n1*n2*n3)-1])
#pragma acc parallel loop collapse(2)// present(grid_3d_2)
	for (Myint64 i3 = i3Start; i3<= i3End; i3++)
	{
		for (Myint64 i2 = i2Start; i2<= i2End; i2++)
		{
			for (Myint64 i1 = i1Start; i1<= i1End; i1++)
			{
				Myint64 ii = i1 + i2*n1 + i3*n2*n1 ;
				grid_3d[ii] = val ;
			}
		}
	}
	//#pragma acc exit data delete(grid_3d)
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
