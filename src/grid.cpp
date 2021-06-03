
//-------------------------------------------------------------------------------------------------------
//
// This grid is activated with command line option -testMode Baseline
//
// ## THIS IS THE MAIN CLASS OF HPCSCAN AND THE REFERENCE IMPLEMENTATION ##
// ## ALL GRID IMPLEMENTATIONS SHOULD DERIVE FROM THIS CLASS             ##
//
// Handle all grid data in hpcscan and all operations performed on grids
//
// ####################################################################################
// ##                                                                                ##
// ##              IMPORTANT: THIS FILE SHOULD NOT BE MODIFIED                       ##
// ## To implement specialization of some functions, you may create a new grid class ##
// ## that derives from this one. See for example: grid_CacheBlk.cpp                 ##
// ##                                                                                ##
// ####################################################################################
//
// Unique definition for 1D, 2D and 3D grids
// For all axis, grid index is as follows:
//
// <-----><------><-----><------><-----><--------> --> Axis N (N=1,2 or 3)
//  Halo1  Layer1  Inner  Layer2  Halo2  Pad
//
// Origin is when Inner starts
// All area have start and end index
// For example, axis 1, halo1:
//  start at point index: i1Halo1Start
//  end at point index  : i1PadEnd
//
// For all dimension (1,2 or 3)
// grid values are stored in the 1D array grid_3d
// Points can be accessed by converting coordinates i1, i2, i3 as
// Myint64 ii = i1 + i2*n1 + i3*n2*n1 ;
//
// where n1, n2, n3 are the dimension of the 3D grid
// For 1D grids, only n1 is meaning full (and n2=n3=1)
// For 2D grids, only n1 and n2 are meaning full (and n3=1)
//-------------------------------------------------------------------------------------------------------

#include "grid.h"

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

Grid::Grid(Grid_type gridTypeIn)
				{
	printDebug(MID_DEBUG, "IN Grid::Grid");

	gridMode = GRID_MODE_BASELINE ;
	n1Inner  = Config::Instance()->n1 ;
	n2Inner  = Config::Instance()->n2 ;
	n3Inner  = Config::Instance()->n3 ;
	gridType = gridTypeIn ;
	dim      = Config::Instance()->dim ;
	grid_3d  = NULL;

	printDebug(MID_DEBUG, "OUT Grid::Grid");
				}

//-------------------------------------------------------------------------------------------------------

Grid::Grid(Grid_type gridTypeIn, Dim_type dimTypeIn,
		Myint64 n1InnerIn, Myint64 n2InnerIn, Myint64 n3InnerIn)
{
	printDebug(MID_DEBUG, "IN Grid::Grid");

	gridMode = GRID_MODE_BASELINE ;
	n1Inner  = n1InnerIn ;
	n2Inner  = n2InnerIn ;
	n3Inner  = n3InnerIn ;
	gridType = gridTypeIn ;
	dim      = dimTypeIn ;
	grid_3d  = NULL;

	printDebug(MID_DEBUG, "OUT Grid::Grid");
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid::initializeGrid(void)
{
	printDebug(MID_DEBUG, "IN Grid::initializeGrid");

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
		subIdx3 = myMpiRank / (Config::Instance()->nsub1*Config::Instance()->nsub2) ;
		subIdx2 = (myMpiRank - subIdx3 *
				(Config::Instance()->nsub1*Config::Instance()->nsub2)) / Config::Instance()->nsub1 ;
		subIdx1 = myMpiRank
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
				i1ProcIdStart = myMpiRank - 1 ;
			}

			if (subIdx1 == Config::Instance()->nsub1-1)
			{
				i1ProcIdEnd = MPI_PROC_NULL ;
			}
			else
			{
				i1ProcIdEnd = myMpiRank + 1 ;
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
				i2ProcIdStart = myMpiRank - Config::Instance()->nsub1 ;
			}

			if (subIdx2 == Config::Instance()->nsub2-1)
			{
				i2ProcIdEnd = MPI_PROC_NULL ;
			}
			else
			{
				i2ProcIdEnd = myMpiRank + Config::Instance()->nsub1 ;
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
				i3ProcIdStart = myMpiRank -
						(Config::Instance()->nsub1 * Config::Instance()->nsub2);
			}

			if (subIdx3 == Config::Instance()->nsub3-1)
			{
				i3ProcIdEnd = MPI_PROC_NULL ;
			}
			else
			{
				i3ProcIdEnd = myMpiRank +
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
	grid_3d = new Myfloat[npoint] ;

	printDebug(MID_DEBUG, "OUT Grid::initializeGrid");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

void Grid::padGridn1(void)
{
	printDebug(MID_DEBUG, "IN Grid::padGridn1");

	if (Config::Instance()->autoPad == true)
	{
		printWarning("Grid::padGridn1, autoPad is not available") ;
		i1PadStart = i1Halo2End ;
		i1PadEnd   = i1Halo2End ;
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

	printDebug(MID_DEBUG, "OUT Grid::padGridn1");
	return ;
}

//-------------------------------------------------------------------------------------------------------

void Grid::padGridn2(void)
{
	printDebug(MID_DEBUG, "IN Grid::padGridn2");

	if (Config::Instance()->autoPad == true)
	{
		printWarning("Grid::padGridn2, autoPad is not available") ;
		i2PadStart = i2Halo2End ;
		i2PadEnd   = i2Halo2End ;
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

	printDebug(MID_DEBUG, "OUT Grid::padGridn2");
}

//-------------------------------------------------------------------------------------------------------

void Grid::padGridn3(void)
{
	printDebug(MID_DEBUG, "IN Grid::padGridn3");

	if (Config::Instance()->autoPad == true)
	{
		printWarning("Grid::padGridn3, autoPad is not available") ;
		i3PadStart = i3Halo2End ;
		i3PadEnd   = i3Halo2End ;
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

	printDebug(MID_DEBUG, "OUT Grid::padGridn3");
}

//-------------------------------------------------------------------------------------------------------

Grid::~Grid(void)
{
	printDebug(MID_DEBUG, "IN Grid::~Grid");

	delete[] grid_3d ;

	printDebug(MID_DEBUG, "OUT Grid::~Grid");
}

//-------------------------------------------------------------------------------------------------------

void Grid::info(void)
{
	printDebug(FULL_DEBUG, "IN Grid::info");

	print_blank() ;
	printInfo(MASTER, " * Grid parameters *") ;

	if (gridType == GRID_LOCAL)
	{
		printInfo(MASTER, " Grid type\t", "LOCAL") ;
	}
	else if (gridType == GRID_GLOBAL)
	{
		printInfo(MASTER, " Grid type\t", "GLOBAL") ;
	}

	printInfo(MASTER, " Grid mode\t", gridMode) ;


	if (dim == DIM1)
	{
		printInfo(MASTER, " Grid dim.\t", "1D") ;
		printInfo(MASTER, " Inner n1\t", i1InnerEnd-i1InnerStart+1) ;
		printInfo(MASTER, " Total n1\t", n1) ;
		printDebug(LIGHT_DEBUG, " Grid sampl. d1\t", d1) ;
		printDebug(LIGHT_DEBUG, " Min. Inner coord1", getMinCoord(AXIS1)) ;
		printDebug(LIGHT_DEBUG, " Max. Inner coord1", getMaxCoord(AXIS1)) ;
	}
	else if (dim == DIM2)
	{
		printInfo(MASTER, " Grid dim.\t", "2D") ;
		printInfo(MASTER, " Total n1\t", n1) ;
		printInfo(MASTER, " Total n2\t", n2) ;
		printInfo(MASTER, " Inner n1\t", i1InnerEnd-i1InnerStart+1) ;
		printInfo(MASTER, " Inner n2\t", i2InnerEnd-i2InnerStart+1) ;
		printDebug(LIGHT_DEBUG, " Grid sampl. d1\t", d1) ;
		printDebug(LIGHT_DEBUG, " Grid sampl. d2\t", d2) ;
		printDebug(LIGHT_DEBUG, " Min. Inner coord1", getMinCoord(AXIS1)) ;
		printDebug(LIGHT_DEBUG, " Max. Inner coord1", getMaxCoord(AXIS1)) ;
		printDebug(LIGHT_DEBUG, " Min. Inner coord2", getMinCoord(AXIS2)) ;
		printDebug(LIGHT_DEBUG, " Max. Inner coord2", getMaxCoord(AXIS2)) ;
	}
	else if (dim == DIM3)
	{
		printInfo(MASTER, " Grid dim.\t", "3D") ;
		printInfo(MASTER, " Total n1\t", n1) ;
		printInfo(MASTER, " Total n2\t", n2) ;
		printInfo(MASTER, " Total n3\t", n3) ;
		printInfo(MASTER, " Inner n1\t", i1InnerEnd-i1InnerStart+1) ;
		printInfo(MASTER, " Inner n2\t", i2InnerEnd-i2InnerStart+1) ;
		printInfo(MASTER, " Inner n3\t", i3InnerEnd-i3InnerStart+1) ;
		printDebug(LIGHT_DEBUG, " Grid sampl. d1\t", d1) ;
		printDebug(LIGHT_DEBUG, " Grid sampl. d2\t", d2) ;
		printDebug(LIGHT_DEBUG, " Grid sampl. d3\t", d3) ;
		printDebug(LIGHT_DEBUG, " Min. Inner coord1", getMinCoord(AXIS1)) ;
		printDebug(LIGHT_DEBUG, " Max. Inner coord1", getMaxCoord(AXIS1)) ;
		printDebug(LIGHT_DEBUG, " Min. Inner coord2", getMinCoord(AXIS2)) ;
		printDebug(LIGHT_DEBUG, " Max. Inner coord2", getMaxCoord(AXIS2)) ;
		printDebug(LIGHT_DEBUG, " Min. Inner coord3", getMinCoord(AXIS3)) ;
		printDebug(LIGHT_DEBUG, " Max. Inner coord3", getMaxCoord(AXIS3)) ;
	}

	Myint nPad = i1PadEnd-i1Halo2End ;
	if (nPad > 0) printInfo(MASTER, " n1 padding\t", nPad) ;
	nPad = i2PadEnd-i2Halo2End ;
	if (nPad > 0) printInfo(MASTER, " n2 padding\t", nPad) ;
	nPad = i3PadEnd-i3Halo2End ;
	if (nPad > 0) printInfo(MASTER, " n3 padding\t", nPad) ;

	Myint64 gridSize = npoint * sizeof(Myfloat) ;
	if (gridSize < 1e3)
	{
		printInfo(MASTER, " Grid size (Byte)", gridSize) ;
	}
	else if (gridSize < 1e6)
	{
		printInfo(MASTER, " Grid size (KB)\t", Myfloat(gridSize/1e3)) ;
	}
	else if (gridSize < 1e9)
	{
		printInfo(MASTER, " Grid size (MB)\t", Myfloat(gridSize/1e6)) ;
	}
	else
	{
		printInfo(MASTER, " Grid size (GB)\t", Myfloat(gridSize/1e9)) ;
	}

	printDebug(LIGHT_DEBUG, "i1Halo1Start" , i1Halo1Start);
	printDebug(LIGHT_DEBUG, "i1Halo1End  " , i1Halo1End);
	printDebug(LIGHT_DEBUG, "i1Layer1Start", i1Layer1Start);
	printDebug(LIGHT_DEBUG, "i1Layer1End  ", i1Layer1End);
	printDebug(LIGHT_DEBUG, "i1InnerStart" , i1InnerStart);
	printDebug(LIGHT_DEBUG, "i1InnerEnd  " , i1InnerEnd);
	printDebug(LIGHT_DEBUG, "i1Layer2Start", i1Layer2Start);
	printDebug(LIGHT_DEBUG, "i1Layer2End  ", i1Layer2End);
	printDebug(LIGHT_DEBUG, "i1Halo2Start" , i1Halo2Start);
	printDebug(LIGHT_DEBUG, "i1Halo2End  " , i1Halo2End);
	printDebug(LIGHT_DEBUG, "i1PadStart"   , i1PadStart);
	printDebug(LIGHT_DEBUG, "i1padEnd  "   , i1PadEnd);

	printDebug(LIGHT_DEBUG, "i2Halo1Start" , i2Halo1Start);
	printDebug(LIGHT_DEBUG, "i2Halo1End  " , i2Halo1End);
	printDebug(LIGHT_DEBUG, "i2Layer1Start", i2Layer1Start);
	printDebug(LIGHT_DEBUG, "i2Layer1End  ", i2Layer1End);
	printDebug(LIGHT_DEBUG, "i2InnerStart" , i2InnerStart);
	printDebug(LIGHT_DEBUG, "i2InnerEnd  " , i2InnerEnd);
	printDebug(LIGHT_DEBUG, "i2Layer2Start", i2Layer2Start);
	printDebug(LIGHT_DEBUG, "i2Layer2End  ", i2Layer2End);
	printDebug(LIGHT_DEBUG, "i2Halo2Start" , i2Halo2Start);
	printDebug(LIGHT_DEBUG, "i2Halo2End  " , i2Halo2End);
	printDebug(LIGHT_DEBUG, "i2PadStart"   , i2PadStart);
	printDebug(LIGHT_DEBUG, "i2padEnd  "   , i2PadEnd);

	printDebug(LIGHT_DEBUG, "i3Halo1Start" , i3Halo1Start);
	printDebug(LIGHT_DEBUG, "i3Halo1End  " , i3Halo1End);
	printDebug(LIGHT_DEBUG, "i3Layer1Start", i3Layer1Start);
	printDebug(LIGHT_DEBUG, "i3Layer1End  ", i3Layer1End);
	printDebug(LIGHT_DEBUG, "i3InnerStart" , i3InnerStart);
	printDebug(LIGHT_DEBUG, "i3InnerEnd  " , i3InnerEnd);
	printDebug(LIGHT_DEBUG, "i3Layer2Start", i3Layer2Start);
	printDebug(LIGHT_DEBUG, "i3Layer2End  ", i3Layer2End);
	printDebug(LIGHT_DEBUG, "i3Halo2Start" , i3Halo2Start);
	printDebug(LIGHT_DEBUG, "i3Halo2End  " , i3Halo2End);
	printDebug(LIGHT_DEBUG, "i3PadStart"   , i3PadStart);
	printDebug(LIGHT_DEBUG, "i3padEnd  "   , i3PadEnd);

	printDebug(LIGHT_DEBUG, "i1ProcIdStart", i1ProcIdStart) ;
	printDebug(LIGHT_DEBUG, "i1ProcIdEnd\t",   i1ProcIdEnd) ;
	printDebug(LIGHT_DEBUG, "i2ProcIdStart", i2ProcIdStart) ;
	printDebug(LIGHT_DEBUG, "i2ProcIdEnd\t",   i2ProcIdEnd) ;
	printDebug(LIGHT_DEBUG, "i3ProcIdStart", i3ProcIdStart) ;
	printDebug(LIGHT_DEBUG, "i3ProcIdEnd\t",   i3ProcIdEnd) ;

	printDebug(FULL_DEBUG, "OUT Grid::info");
}

//-------------------------------------------------------------------------------------------------------

void Grid::write(string file_name)
{
	printDebug(LIGHT_DEBUG, "IN Grid::write");

	// each proc write is own file

	if (Config::Instance()->writeGrid)
	{
		// write grid in binary format
		string file_name_bin = file_name + ".proc" + to_string(myMpiRank) + ".grid.bin" ;
		printInfo(MASTER, "* write on disk\t", file_name_bin.c_str());
		printDebug(LIGHT_DEBUG, "* write on disk\t", file_name_bin.c_str());
		ofstream out_file ;
		out_file.open(file_name_bin.c_str(), ios::binary | ios::app) ;
		assert(out_file.is_open());
		out_file.write((char*) grid_3d, npoint * sizeof(Myfloat)) ;
		out_file.close() ;

		// write grid info
		string file_name_info = file_name + ".proc" + to_string(myMpiRank) + ".grid.info" ;
		printDebug(LIGHT_DEBUG, "* write on disk\t", file_name_info.c_str());
		ofstream out_file2(file_name_info) ;
		out_file2 << this->n1 << "\n" ;
		out_file2 << this->n2 << "\n" ;
		out_file2 << this->n3 << "\n" ;
		out_file2.close() ;
	}

	printDebug(LIGHT_DEBUG, "OUT Grid::write");
}

//-------------------------------------------------------------------------------------------------------

void Grid::fill(Point_type pointType, Func_type t1,  Func_type t2, Func_type t3,
		Myfloat64 param1, Myfloat64 param2, Myfloat64 param3, Myfloat64 amp)
{
	printDebug(MID_DEBUG, "IN Grid::fill");

	Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End, s1, s2, s3 ;
	getGridIndex(pointType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End) ;

	s1 = i1End - i1Start + 1;
	s2 = i2End - i2Start + 1;
	s3 = i3End - i3Start + 1;

	Myfloat64 *val1 = new Myfloat64[s1] ;
	Myfloat64 *val2 = new Myfloat64[s2] ;
	Myfloat64 *val3 = new Myfloat64[s3] ;

#pragma omp simd
	for (Myint64 i1 = i1Start; i1<= i1End; i1++)
	{
		Myfloat64 coord1 = Myfloat64(Orig1 + i1 * d1) ;

		if (dim >= DIM1)
		{
			if (t1 == FUNC_SINE)
			{
				val1[i1-i1Start] = sin(coord1 * param1) ;
			}
			else if (t1 == FUNC_COSINE)
			{
				val1[i1-i1Start] = cos(coord1 * param1) ;
			}
			else if (t1 == FUNC_LINEAR)
			{
				val1[i1-i1Start] = coord1 ;
			}
			else if (t1 == FUNC_CONST)
			{
				val1[i1-i1Start] = 1.0 ;
			}
			else if (t1 == FUNC_RANDOM)
			{
				val1[i1-i1Start] = Myfloat(rand()) / RAND_MAX ;
			}
			else
			{
				val1[i1-i1Start] = 1.0;
			}
		}
		else
		{
			val1[i1-i1Start] = 1.0 ;
		}
	}

#pragma omp simd
	for (Myint64 i2 = i2Start; i2<= i2End; i2++)
	{
		Myfloat64 coord2 = Myfloat64(Orig2 + i2 * d2) ;
		if (dim >= DIM2)
		{
			if (t2 == FUNC_SINE)
			{
				val2[i2-i2Start] = sin(coord2 * param2) ;
			}
			else if (t2 == FUNC_COSINE)
			{
				val2[i2-i2Start] = cos(coord2 * param2) ;
			}
			else if (t2 == FUNC_LINEAR)
			{
				val2[i2-i2Start] = coord2 ;
			}
			else if (t2 == FUNC_CONST)
			{
				val2[i2-i2Start] = 1.0 ;
			}
			else if (t2 == FUNC_RANDOM)
			{
				val2[i2-i2Start] = Myfloat(rand()) / RAND_MAX ;
			}
			else
			{
				val2[i2-i2Start] = 1.0;
			}
		}
		else
		{
			val2[i2-i2Start] = 1.0 ;
		}
	}

#pragma omp simd
	for (Myint64 i3 = i3Start; i3<= i3End; i3++)
	{
		Myfloat64 coord3 = Myfloat64(Orig3 + i3 * d3) ;
		if (dim >= DIM3)
		{
			if (t3 == FUNC_SINE)
			{
				val3[i3-i3Start] = sin(coord3 * param3) ;
			}
			else if (t3 == FUNC_COSINE)
			{
				val3[i3-i3Start] = cos(coord3 * param3) ;
			}
			else if (t3 == FUNC_LINEAR)
			{
				val3[i3-i3Start] = coord3 ;
			}
			else if (t3 == FUNC_CONST)
			{
				val3[i3-i3Start] = 1.0 ;
			}
			else if (t3 == FUNC_RANDOM)
			{
				val3[i3-i3Start] = Myfloat(rand()) / RAND_MAX ;
			}
			else
			{
				val3[i3-i3Start] = 1.0;
			}
		}
		else
		{
			val3[i3-i3Start] = 1.0 ;
		}
	}

#pragma omp parallel for collapse(2)
	for (Myint64 i3 = i3Start; i3<= i3End; i3++)
	{
		for (Myint64 i2 = i2Start; i2<= i2End; i2++)
		{
#pragma omp simd
			for (Myint64 i1 = i1Start; i1<= i1End; i1++)
			{
				// assign value at current grid point
				grid_3d[i1+i2*n1+i3*n1*n2] = amp * val1[i1-i1Start] * val2[i2-i2Start] * val3[i3-i3Start] ;

			}
		}
	}

	delete[] val1 ;
	delete[] val2 ;
	delete[] val3 ;

	printDebug(MID_DEBUG, "OUT Grid::fill") ;
}

//-------------------------------------------------------------------------------------------------------

void Grid::fill(Point_type pointType, Myfloat val)
{
	printDebug(MID_DEBUG, "IN Grid::fill");

	Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
	getGridIndex(pointType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End) ;

#pragma omp parallel for collapse(2)
	for (Myint64 i3 = i3Start; i3<= i3End; i3++)
	{
		for (Myint64 i2 = i2Start; i2<= i2End; i2++)
		{
#pragma omp simd
			for (Myint64 i1 = i1Start; i1<= i1End; i1++)
			{
				Myint64 ii = i1 + i2*n1 + i3*n2*n1 ;
				grid_3d[ii] = val ;
			}
		}
	}

	printDebug(MID_DEBUG, "OUT Grid::fill");
}

//-------------------------------------------------------------------------------------------------------

void Grid::fillArray(Myfloat val)
{
	printDebug(MID_DEBUG, "IN Grid::fillArray");

	Myfloat * const w = grid_3d ;

#pragma omp parallel for
	for (Myint64 ii=0; ii<npoint; ii++)
	{
		w[ii] = val ;
	}

	printDebug(MID_DEBUG, "OUT Grid::fillArray");
}

//-------------------------------------------------------------------------------------------------------

void Grid::copyArray(const Grid& gridIn)
{
	printDebug(MID_DEBUG, "IN Grid::copyArray");

	Myfloat * const w = grid_3d ;
	Myfloat * const u = gridIn.grid_3d ;

#pragma omp parallel for
	for (Myint64 ii=0; ii<npoint; ii++)
	{
		w[ii] = u[ii] ;
	}

	printDebug(MID_DEBUG, "OUT Grid::copyArray");
}

//-------------------------------------------------------------------------------------------------------

void Grid::addArray(const Grid& gridIn1, const Grid& gridIn2)
{
	printDebug(MID_DEBUG, "IN Grid::addArray");

	Myfloat * const w = grid_3d ;
	Myfloat * const u = gridIn1.grid_3d ;
	Myfloat * const v = gridIn2.grid_3d ;

#pragma omp parallel for
	for (Myint64 ii=0; ii<npoint; ii++)
	{
		w[ii] = u[ii] + v[ii];
	}

	printDebug(MID_DEBUG, "OUT Grid::addArray");
}

//-------------------------------------------------------------------------------------------------------

void Grid::multiplyArray(const Grid& gridIn1, const Grid& gridIn2)
{
	printDebug(MID_DEBUG, "IN Grid::multiplyArray");

	Myfloat * const w = grid_3d ;
	Myfloat * const u = gridIn1.grid_3d ;
	Myfloat * const v = gridIn2.grid_3d ;

#pragma omp parallel for
	for (Myint64 ii=0; ii<npoint; ii++)
	{
		w[ii] = u[ii] * v[ii];
	}

	printDebug(MID_DEBUG, "OUT Grid::multiplyArray");
}

//-------------------------------------------------------------------------------------------------------

void Grid::addUpdateArray(const Grid& gridIn)
{
	printDebug(MID_DEBUG, "IN Grid::addUpdateArray");

	Myfloat * const w = grid_3d ;
	Myfloat * const u = gridIn.grid_3d ;

#pragma omp parallel for
	for (Myint64 ii=0; ii<npoint; ii++)
	{
		w[ii] = w[ii] + u[ii];
	}

	printDebug(MID_DEBUG, "OUT Grid::addUpdateArray");
}

//-------------------------------------------------------------------------------------------------------

// define a unit grid, with for the 3 axis:
// min. coordinate of Inner domain = 0 (in Global grid)
// max. coordinate of Inner domain = 1 (in Global grid)
// set d1, d2, d3, Orig1, Orig2, Orig3

void Grid::defineUnitGrid(void)
{
	printDebug(FULL_DEBUG, "IN Grid::defineUnitGrid");

	// save previous sampling
	Myfloat64 d1Old = d1 ;
	Myfloat64 d2Old = d2 ;
	Myfloat64 d3Old = d3 ;

	// define new spatial sampling
	Myfloat64 extentGrid = 1.0 ;
	d1 = extentGrid / (Config::Instance()->n1-1) ;
	if (dim >= DIM2) d2 = extentGrid / (Config::Instance()->n2-1) ;
	if (dim >= DIM3) d3 = extentGrid / (Config::Instance()->n3-1) ;

	// reset origin of grid
	Orig1 = Orig1 * d1/d1Old ;
	if (dim >= DIM2) Orig2 = Orig2 * d2/d2Old ;
	if (dim >= DIM3) Orig3 = Orig3 * d3/d3Old ;

	printDebug(LIGHT_DEBUG, "d1" , d1);
	printDebug(LIGHT_DEBUG, "d2" , d2);
	printDebug(LIGHT_DEBUG, "d3" , d3);

	printDebug(LIGHT_DEBUG, "Orig1" , Orig1);
	printDebug(LIGHT_DEBUG, "Orig2" , Orig2);
	printDebug(LIGHT_DEBUG, "Orig3" , Orig3);

	printDebug(FULL_DEBUG, "OUT Grid::defineUnitGrid");
}

//-------------------------------------------------------------------------------------------------------

Myfloat Grid::maxErr(Point_type pointType, const Grid& gridIn) const
{
	printDebug(FULL_DEBUG, "IN Grid::maxErr");

	// check grids have same size
	if (!(this->sameSize(gridIn)))
	{
		printError("Grid::maxErr, grids have different size") ;
		return(-1.0) ;
	}

	Myfloat err = -FLT_MAX, err2 = 0.0 ;

	Myfloat* u1 = this->grid_3d ;
	Myfloat* u2 = gridIn.grid_3d ;

	Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
	getGridIndex(pointType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End) ;

#pragma omp parallel for collapse(2) reduction(max:err)
	for (Myint64 i3 = i3Start; i3<= i3End; i3++)
	{
		for (Myint64 i2 = i2Start; i2<= i2End; i2++)
		{
			for (Myint64 i1 = i1Start; i1<= i1End; i1++)
			{
				Myint64 ii = i1 + i2*n1 + i3*n2*n1 ;

				// prevent divide by 0
				if (fabs(u2[ii]) < MAX_ERR_FLOAT)
				{
					err2 = fabs(u1[ii] - u2[ii]) ;
				}
				else
				{
					err2 = fabs(u1[ii] - u2[ii]) / fabs(u2[ii]) ;
				}

				if (err2 > err)
				{
					err = err2 ;
				}
			}
		}
	}

	printDebug(FULL_DEBUG, "OUT Grid::maxErr");
	return(err) ;
}

//-------------------------------------------------------------------------------------------------------

Myfloat Grid::L1Err(Point_type pointType, const Grid& gridIn) const
{
	printDebug(LIGHT_DEBUG, "IN Grid::L1Err");

	// check grids have same size
	if (!(this->sameSize(gridIn)))
	{
		printError("Grid::L1Err, grids have different size") ;
		return(-1.0) ;
	}

	Myfloat64 sum1 = 0.0 ;
	Myfloat64 sum2 = 0.0 ;

	Myfloat* u1 = this->grid_3d ;
	Myfloat* u2 = gridIn.grid_3d ;

	Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
	getGridIndex(pointType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End) ;


#pragma omp parallel for collapse(2) reduction(+:sum1,sum2)
	for (Myint64 i3 = i3Start; i3<= i3End; i3++)
	{
		for (Myint64 i2 = i2Start; i2<= i2End; i2++)
		{
			for (Myint64 i1 = i1Start; i1<= i1End; i1++)
			{
				Myint64 ii = i1 + i2*n1 + i3*n2*n1 ;

				sum1 += fabs(u1[ii] - u2[ii]) ;
				sum2 += fabs(u2[ii]) ;
			}
		}
	}

	// prevent divide by zero
	if (sum2 < MAX_ERR_FLOAT) sum2 = 1.0 * npoint ;
	Myfloat err = sum1 / sum2 ;

	printDebug(LIGHT_DEBUG, "sum1", sum1) ;
	printDebug(LIGHT_DEBUG, "sum2", sum2) ;
	printDebug(LIGHT_DEBUG, "err", err) ;

	if (std::isnan(err))
	{
		printError("In Grid::L1Err, std::isnan(err)") ;
	}

	printDebug(LIGHT_DEBUG, "OUT Grid::L1Err");
	return(err) ;
}

//-------------------------------------------------------------------------------------------------------

Myfloat Grid::allProcL1Err(Point_type pointType, const Grid& gridIn) const
{
	printDebug(LIGHT_DEBUG, "IN Grid::allProcL1Err");

	// check grids have same size
	if (!(this->sameSize(gridIn)))
	{
		printError("Grid::allProcL1Err, grids have different size") ;
		return(-1.0) ;
	}

	Myfloat64 sum1Loc = 0.0 ;
	Myfloat64 sum2Loc = 0.0 ;
	Myfloat64 sum1 = 0.0 ;
	Myfloat64 sum2 = 0.0 ;

	Myfloat* u1 = this->grid_3d ;
	Myfloat* u2 = gridIn.grid_3d ;

	Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
	getGridIndex(pointType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End) ;

#pragma omp parallel for collapse(2) reduction(+:sum1Loc,sum2Loc)
	for (Myint64 i3 = i3Start; i3<= i3End; i3++)
	{
		for (Myint64 i2 = i2Start; i2<= i2End; i2++)
		{
			for (Myint64 i1 = i1Start; i1<= i1End; i1++)
			{
				Myint64 ii = i1 + i2*n1 + i3*n2*n1 ;

				sum1Loc += fabs(u1[ii] - u2[ii]) ;
				sum2Loc += fabs(u2[ii]) ;
			}
		}
	}

	// MPI reduction
	MPI_Reduce(&sum1Loc, &sum1, 1, MPI_MYFLOAT64, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&sum2Loc, &sum2, 1, MPI_MYFLOAT64, MPI_SUM, 0, MPI_COMM_WORLD);

	// prevent divide by zero
	if (sum2 == 0.0) sum2 = 1.0 * npoint ;
	Myfloat err = sum1 / sum2 ;

	printDebug(LIGHT_DEBUG, "sum1", sum1) ;
	printDebug(LIGHT_DEBUG, "sum2", sum2) ;
	printDebug(LIGHT_DEBUG, "err", err) ;

	if (std::isnan(err))
	{
		printError("In Grid::allProcL1Err, std::isnan(err)") ;
	}

	printDebug(LIGHT_DEBUG, "OUT Grid::allProcL1Err");
	return(err) ;
}

//-------------------------------------------------------------------------------------------------------

Myfloat Grid::getSumAbsDiff(Point_type pointType, const Grid& gridIn) const
{
	printDebug(LIGHT_DEBUG, "IN Grid::getSumAbsDiff");

	// check grids have same size
	if (!(this->sameSize(gridIn)))
	{
		printError("Grid::getSumAbsDiff, grids have different size") ;
		return(-1.0) ;
	}

	Myfloat64 sum1Loc = 0.0 ;
	Myfloat64 sum1 = 0.0 ;

	Myfloat* u1 = this->grid_3d ;
	Myfloat* u2 = gridIn.grid_3d ;

	Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
	getGridIndex(pointType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End) ;

#pragma omp parallel for collapse(2) reduction(+:sum1Loc)
	for (Myint64 i3 = i3Start; i3<= i3End; i3++)
	{
		for (Myint64 i2 = i2Start; i2<= i2End; i2++)
		{
			for (Myint64 i1 = i1Start; i1<= i1End; i1++)
			{
				Myint64 ii = i1 + i2*n1 + i3*n2*n1 ;

				sum1Loc += fabs(u1[ii] - u2[ii]) ;
			}
		}
	}

	// reduction
	if (gridType == GRID_LOCAL)
	{
		MPI_Reduce(&sum1Loc, &sum1, 1, MPI_MYFLOAT64, MPI_SUM, 0, MPI_COMM_WORLD);
	}
	else
	{
		sum1 = sum1Loc ;
	}

	printDebug(LIGHT_DEBUG, "sum1", sum1) ;

	if (std::isnan(sum1))
	{
		printError("In Grid::getSumAbsDiff, std::isnan(sum1)") ;
	}

	printDebug(LIGHT_DEBUG, "OUT Grid::getSumAbsDiff");
	return(sum1) ;
}

//-------------------------------------------------------------------------------------------------------

Myfloat Grid::getSumAbs(Point_type pointType) const
{
	printDebug(LIGHT_DEBUG, "IN Grid::getSumAbs");

	Myfloat64 sum2Loc = 0.0 ;
	Myfloat64 sum2 = 0.0 ;

	Myfloat* u1 = this->grid_3d ;

	Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
	getGridIndex(pointType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End) ;


#pragma omp parallel for collapse(2) reduction(+:sum2Loc)
	for (Myint64 i3 = i3Start; i3<= i3End; i3++)
	{
		for (Myint64 i2 = i2Start; i2<= i2End; i2++)
		{
			for (Myint64 i1 = i1Start; i1<= i1End; i1++)
			{
				Myint64 ii = i1 + i2*n1 + i3*n2*n1 ;

				sum2Loc += fabs(u1[ii]) ;
			}
		}
	}

	// reduction
	if (gridType == GRID_LOCAL)
	{
		MPI_Reduce(&sum2Loc, &sum2, 1, MPI_MYFLOAT64, MPI_SUM, 0, MPI_COMM_WORLD);
	}
	else
	{
		sum2 = sum2Loc ;
	}

	printDebug(LIGHT_DEBUG, "sum2", sum2) ;

	if (std::isnan(sum2))
	{
		printError("In Grid::getSumAbs, std::isnan(sum2)") ;
	}

	printDebug(LIGHT_DEBUG, "OUT Grid::getSumAbs");
	return(sum2) ;
}

//-------------------------------------------------------------------------------------------------------

Myfloat Grid::getMin(Point_type pointType)
{
	printDebug(LIGHT_DEBUG, "IN Grid::getMin");

	Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
	getGridIndex(pointType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End) ;
	Myfloat val = FLT_MAX ;

#pragma omp parallel for collapse(2) reduction(min:val)
	for (Myint64 i3 = i3Start; i3<= i3End; i3++)
	{
		for (Myint64 i2 = i2Start; i2<= i2End; i2++)
		{
			for (Myint64 i1 = i1Start; i1<= i1End; i1++)
			{
				Myint64 ii = i1 + i2*n1 + i3*n2*n1 ;
				if (grid_3d[ii] < val) val = grid_3d[ii] ;
			}
		}
	}

	printDebug(LIGHT_DEBUG, "Min val", val);

	printDebug(LIGHT_DEBUG, "OUT Grid::getMin");

	return(val) ;
}

//-------------------------------------------------------------------------------------------------------

Myfloat Grid::getMax(Point_type pointType)
{
	printDebug(LIGHT_DEBUG, "IN Grid::getMax");

	Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
	getGridIndex(pointType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End) ;

	Myfloat val = -FLT_MAX ;

#pragma omp parallel for collapse(2) reduction(max:val)
	for (Myint64 i3 = i3Start; i3<= i3End; i3++)
	{
		for (Myint64 i2 = i2Start; i2<= i2End; i2++)
		{
			for (Myint64 i1 = i1Start; i1<= i1End; i1++)
			{
				Myint64 ii = i1 + i2*n1 + i3*n2*n1 ;
				if (grid_3d[ii] > val) val = grid_3d[ii] ;
			}
		}
	}

	printDebug(LIGHT_DEBUG, "Max val", val);

	printDebug(LIGHT_DEBUG, "OUT Grid::getMax");
	return(val) ;
}

//-------------------------------------------------------------------------------------------------------

Myint Grid::getFlopPerPtFD_D2(Myint fdOrder)
{
	printDebug(FULL_DEBUG, "IN Grid::getFlopPerPtFD_D2");

	Myint nOpDer=0, val=0 ;

	if (fdOrder == 2)
	{
		nOpDer = FD_D2_O2_NOP ;
	}
	else if (fdOrder == 4)
	{
		nOpDer = FD_D2_O4_NOP ;
	}
	else if (fdOrder == 8)
	{
		nOpDer = FD_D2_O8_NOP ;
	}
	else if (fdOrder == 12)
	{
		nOpDer = FD_D2_O12_NOP ;
	}
	else if (fdOrder == 16)
	{
		nOpDer = FD_D2_O16_NOP ;
	}

	if (dim ==DIM1) val = nOpDer ;
	else if (dim == DIM2) val = nOpDer ;
	else if (dim == DIM3) val = nOpDer ;

	printDebug(FULL_DEBUG, "OUT Grid::getFlopPerPtFD_D2");
	return(val) ;
}

//-------------------------------------------------------------------------------------------------------

Myint Grid::getFlopPerPtFD_LAPLACIAN(Myint fdOrder)
{
	printDebug(FULL_DEBUG, "IN Grid::getFlopPerPtFD_LAPLACIAN");

	Myint nOpDer=0, val=0 ;

	nOpDer = getFlopPerPtFD_D2(fdOrder) ;

	if (dim ==DIM1) val = nOpDer ;
	else if (dim == DIM2) val = 2*nOpDer + 1 ; // + 1 ADD
	else if (dim == DIM3) val = 3*nOpDer + 2 ; // + 2 ADD

	printDebug(FULL_DEBUG, "OUT Grid::getFlopPerPtFD_LAPLACIAN");
	return(val) ;
}

//-------------------------------------------------------------------------------------------------------

Myint Grid::getPtPerStencilFD_D2(Myint fdOrder)
{
	printDebug(FULL_DEBUG, "IN Grid::getPtPerStencilFD_D2");

	Myint nHsDer = fdOrder / 2 ; // half stencil length
	Myint val = nHsDer * 2 + 1 ;

	printDebug(FULL_DEBUG, "OUT Grid::getPtPerStencilFD_D2");
	return(val) ;
}

//-------------------------------------------------------------------------------------------------------

Myint Grid::getPtPerStencilFD_LAPLACIAN(Myint fdOrder)
{
	printDebug(FULL_DEBUG, "IN Grid::getPtPerStencilFD_LAPLACIAN");

	Myint nHsDer = fdOrder / 2 ; // half stencil length
	Myint val = 0 ;

	if (dim ==DIM1) val = nHsDer * 2 + 1 ;
	else if (dim == DIM2) val = nHsDer * 4 + 1 ;
	else if (dim == DIM3) val = nHsDer * 6 + 1 ;

	printDebug(FULL_DEBUG, "OUT Grid::getPtPerStencilFD_LAPLACIAN");
	return(val) ;
}


//-------------------------------------------------------------------------------------------------------

Myint Grid::getNumberOfNeighbour(void)
{
	printDebug(FULL_DEBUG, "IN Grid::getNumberOfNeighbour");

	Myint val = 0 ;

	if (i1ProcIdStart != MPI_PROC_NULL) val++ ;
	if (i1ProcIdEnd   != MPI_PROC_NULL) val++ ;
	if (i2ProcIdStart != MPI_PROC_NULL) val++ ;
	if (i2ProcIdEnd   != MPI_PROC_NULL) val++ ;
	if (i3ProcIdStart != MPI_PROC_NULL) val++ ;
	if (i3ProcIdEnd   != MPI_PROC_NULL) val++ ;

	printDebug(FULL_DEBUG, "OUT Grid::getNumberOfNeighbour");
	return(val) ;
}

//-------------------------------------------------------------------------------------------------------

Myfloat Grid::getMinCoord(Axis_type axisType)
{
	printDebug(MID_DEBUG, "IN Grid::getMinCoord");
	Myfloat retVal = -1 ;
	if (axisType == AXIS1) retVal = Orig1 + d1*i1InnerStart ;
	else if (axisType == AXIS2) retVal = Orig2 + d2*i2InnerStart ;
	else if (axisType == AXIS3) retVal = Orig3 + d3*i3InnerStart ;
	printDebug(LIGHT_DEBUG, "retVal", retVal) ;
	printDebug(MID_DEBUG, "OUT Grid::getMinCoord");
	return(retVal) ;
}

//-------------------------------------------------------------------------------------------------------

Myfloat Grid::getMaxCoord(Axis_type axisType)
{
	printDebug(MID_DEBUG, "IN Grid::getMaxCoord");
	Myfloat retVal = -1 ;
	if (axisType == AXIS1) retVal = Orig1 + d1*i1InnerEnd ;
	else if (axisType == AXIS2) retVal = Orig2 + d2*i2InnerEnd ;
	else if (axisType == AXIS3) retVal = Orig3 + d3*i3InnerEnd ;
	printDebug(LIGHT_DEBUG, "retVal", retVal) ;
	printDebug(MID_DEBUG, "OUT Grid::getMaxCoord");
	return(retVal) ;
}

//-------------------------------------------------------------------------------------------------------

Myfloat Grid::getMinSpaceSampling()
{
	printDebug(MID_DEBUG, "IN Grid::getMinSpaceSampling");

	Myfloat minSpaceSampling = d1 ;
	if (dim >= DIM2) minSpaceSampling = min((double) minSpaceSampling, d2) ;
	if (dim >= DIM3) minSpaceSampling = min((double) minSpaceSampling, d3) ;

	return(minSpaceSampling) ;
}

//-------------------------------------------------------------------------------------------------------

Myint Grid::getNeighbourProc(Point_type pointType)
{
	printDebug(MID_DEBUG, "IN Grid::getNeighbourProc");
	Myint retVal = MPI_PROC_NULL;
	if (pointType == I1HALO1) retVal = i1ProcIdStart ;
	else if (pointType == I1HALO2) retVal = i1ProcIdEnd ;
	else if (pointType == I2HALO1) retVal = i2ProcIdStart ;
	else if (pointType == I2HALO2) retVal = i2ProcIdEnd ;
	else if (pointType == I3HALO1) retVal = i3ProcIdStart ;
	else if (pointType == I3HALO2) retVal = i3ProcIdEnd ;
	else
	{
		printError("IN Grid::getNeighbourProc, invalid pointType") ;
	}
	printDebug(MID_DEBUG, "IN Grid::getNeighbourProc");
	return(retVal) ;
}

//-------------------------------------------------------------------------------------------------------

Myint64 Grid::getNumberOfGridPointCommHalo(Grid_type gridTypeIn)
{
	printDebug(FULL_DEBUG, "IN Grid::getNumberOfGridPointCommHalo");

	Myint64 nGridPointLoc = 0 ;
	if (i1ProcIdStart != MPI_PROC_NULL) nGridPointLoc += getNumberOfGridPoint(GRID_LOCAL, I1HALO1) ;
	if (i1ProcIdEnd   != MPI_PROC_NULL) nGridPointLoc += getNumberOfGridPoint(GRID_LOCAL, I1HALO2) ;
	if (i2ProcIdStart != MPI_PROC_NULL) nGridPointLoc += getNumberOfGridPoint(GRID_LOCAL, I2HALO1) ;
	if (i2ProcIdEnd   != MPI_PROC_NULL) nGridPointLoc += getNumberOfGridPoint(GRID_LOCAL, I2HALO2) ;
	if (i3ProcIdStart != MPI_PROC_NULL) nGridPointLoc += getNumberOfGridPoint(GRID_LOCAL, I3HALO1) ;
	if (i3ProcIdEnd   != MPI_PROC_NULL) nGridPointLoc += getNumberOfGridPoint(GRID_LOCAL, I3HALO2) ;

	Myint64 nGridPoint = 0 ;
	if ((gridTypeIn == GRID_GLOBAL) && (gridType == GRID_LOCAL))
	{
		//MPI_Allreduce(&nGridPointLoc, &nGridPoint, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
		MPI_Reduce(&nGridPointLoc, &nGridPoint, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	}
	else
	{
		nGridPoint = nGridPointLoc ;
	}

	printDebug(FULL_DEBUG, "OUT Grid::getNumberOfGridPointCommHalo");
	return(nGridPoint) ;
}

//-------------------------------------------------------------------------------------------------------

Myint64 Grid::getNumberOfGridPointBoundaryCondition(BoundCond_type boundCondType)
{
	printDebug(FULL_DEBUG, "IN Grid::getNumberOfGridPointBoundaryCondition");
	Myint64 retVal = 0 ;

	if (boundCondType == NO_BOUND_COND)
	{
		// no point
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
				retVal += (iInner1-i1Start+1)*(i2End-i2Start+1)*(i3End-i3Start+1) ;
			}

			// I1HALO2
			if (getNeighbourProc(I1HALO2) == MPI_PROC_NULL)
			{
				getGridIndex(I1HALO2, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End) ;
				Myint64 iInner1 = i1Start-1 ;
				retVal += (i1End-iInner1+1)*(i2End-i2Start+1)*(i3End-i3Start+1) ;
			}
		}

		if (dim >= DIM2)
		{
			// I2HALO1
			if (getNeighbourProc(I2HALO1) == MPI_PROC_NULL)
			{
				getGridIndex(I2HALO1, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End) ;
				Myint64 iInner2 = i2End+1 ;
				retVal += (i1End-i1Start+1)*(iInner2-i2Start+1)*(i3End-i3Start+1) ;
			}

			// I2HALO2
			if (getNeighbourProc(I2HALO2) == MPI_PROC_NULL)
			{
				getGridIndex(I2HALO2, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End) ;
				Myint64 iInner2 = i2Start-1 ;
				retVal += (i1End-i1Start+1)*(i2End-iInner2+1)*(i3End-i3Start+1) ;
			}
		}

		if (dim >= DIM3)
		{
			// I3HALO1
			if (getNeighbourProc(I3HALO1) == MPI_PROC_NULL)
			{
				getGridIndex(I3HALO1, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End) ;
				Myint64 iInner3 = i3End+1 ;
				retVal += (i1End-i1Start+1)*(i2End-i2Start+1)*(iInner3-i3Start+1) ;
			}

			// I3HALO2
			if (getNeighbourProc(I3HALO2) == MPI_PROC_NULL)
			{
				getGridIndex(I3HALO2, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End) ;
				Myint64 iInner3 = i3Start-1 ;
				retVal += (i1End-i1Start+1)*(i2End-i2Start+1)*(i3End-iInner3+1) ;
			}
		}
	}
	else
	{
		printError("IN Grid::getNumberOfGridPointBoundaryCondition", boundCondType) ;
		retVal = -1 ;
	}

	printDebug(FULL_DEBUG, "OUT Grid::getNumberOfGridPointBoundaryCondition");
	return(retVal) ;
}

//-------------------------------------------------------------------------------------------------------
Rtn_code Grid::getGridIndex(Point_type pointType, Myint64* i1Start, Myint64* i1End,
		Myint64* i2Start, Myint64* i2End, Myint64* i3Start, Myint64* i3End) const
{
	printDebug(FULL_DEBUG, "IN Grid::getGridIndex");

	if (pointType == INNER_POINTS)
	{
		*i1Start = i1InnerStart ;
		*i1End   = i1InnerEnd ;
		*i2Start = i2InnerStart ;
		*i2End   = i2InnerEnd ;
		*i3Start = i3InnerStart ;
		*i3End   = i3InnerEnd ;
	}
	else if (pointType == ALL_POINTS)
	{
		*i1Start = i1Halo1Start ;
		*i1End   = i1Halo2End ;
		*i2Start = i2Halo1Start ;
		*i2End   = i2Halo2End ;
		*i3Start = i3Halo1Start ;
		*i3End   = i3Halo2End ;
	}
	else if (pointType == MIDDLE_POINT)
	{
		*i1Start = n1 / 2;
		*i1End   = *i1Start ;
		*i2Start = n2 / 2  ;
		*i2End   = *i2Start ;
		*i3Start = n3 / 2 ;
		*i3End   = *i3Start ;
	}
	else if (pointType == I1HALO1)
	{
		*i1Start = i1Halo1Start ;
		*i1End   = i1Halo1End ;
		*i2Start = i2InnerStart ;
		*i2End   = i2InnerEnd ;
		*i3Start = i3InnerStart ;
		*i3End   = i3InnerEnd ;
	}
	else if (pointType == I1HALO2)
	{
		*i1Start = i1Halo2Start ;
		*i1End   = i1Halo2End ;
		*i2Start = i2InnerStart ;
		*i2End   = i2InnerEnd ;
		*i3Start = i3InnerStart ;
		*i3End   = i3InnerEnd ;
	}
	else if (pointType == I2HALO1)
	{
		*i1Start = i1InnerStart ;
		*i1End   = i1InnerEnd ;
		*i2Start = i2Halo1Start ;
		*i2End   = i2Halo1End ;
		*i3Start = i3InnerStart ;
		*i3End   = i3InnerEnd ;
	}
	else if (pointType == I2HALO2)
	{
		*i1Start = i1InnerStart ;
		*i1End   = i1InnerEnd ;
		*i2Start = i2Halo2Start ;
		*i2End   = i2Halo2End ;
		*i3Start = i3InnerStart ;
		*i3End   = i3InnerEnd ;
	}
	else if (pointType == I3HALO1)
	{
		*i1Start = i1InnerStart ;
		*i1End   = i1InnerEnd ;
		*i2Start = i2InnerStart ;
		*i2End   = i2InnerEnd ;
		*i3Start = i3Halo1Start ;
		*i3End   = i3Halo1End ;
	}
	else if (pointType == I3HALO2)
	{
		*i1Start = i1InnerStart ;
		*i1End   = i1InnerEnd ;
		*i2Start = i2InnerStart ;
		*i2End   = i2InnerEnd ;
		*i3Start = i3Halo2Start ;
		*i3End   = i3Halo2End ;
	}
	else if (pointType == I1INNERHALO1)
	{
		*i1Start = i1InnerStart ;
		*i1End   = *i1Start + haloWidth ;
		*i2Start = i2InnerStart ;
		*i2End   = i2InnerEnd ;
		*i3Start = i3InnerStart ;
		*i3End   = i3InnerEnd ;
	}
	else if (pointType == I1INNERHALO2)
	{
		*i1Start = i1InnerEnd - haloWidth ;
		*i1End   = *i1Start + haloWidth ;
		*i2Start = i2InnerStart ;
		*i2End   = i2InnerEnd ;
		*i3Start = i3InnerStart ;
		*i3End   = i3InnerEnd ;
	}
	else if (pointType == I2INNERHALO1)
	{
		*i1Start = i1InnerStart ;
		*i1End   = i1InnerEnd ;
		*i2Start = i2InnerStart ;
		*i2End   = *i2Start + haloWidth ;
		*i3Start = i3InnerStart ;
		*i3End   = i3InnerEnd ;
	}
	else if (pointType == I2INNERHALO2)
	{
		*i1Start = i1InnerStart ;
		*i1End   = i1InnerEnd ;
		*i2Start = i2InnerEnd - haloWidth ;
		*i2End   = *i2Start + haloWidth ;
		*i3Start = i3InnerStart ;
		*i3End   = i3InnerEnd ;
	}
	else if (pointType == I3INNERHALO1)
	{
		*i1Start = i1InnerStart ;
		*i1End   = i1InnerEnd ;
		*i2Start = i2InnerStart ;
		*i2End   = i2InnerEnd ;
		*i3Start = i3InnerStart ;
		*i3End   = *i3Start + haloWidth ;
	}
	else if (pointType == I3INNERHALO2)
	{
		*i1Start = i1InnerStart ;
		*i1End   = i1InnerEnd ;
		*i2Start = i2InnerStart ;
		*i2End   = i2InnerEnd ;
		*i3Start = i3InnerEnd - haloWidth ;
		*i3End   = *i3Start + haloWidth ;
	}
	else
	{
		printError("IN Grid::getGridIndex, Invalid pointType") ;
	}

	printDebug(MID_DEBUG, "i1Start", *i1Start);
	printDebug(MID_DEBUG, "i1End",   *i1End);
	printDebug(MID_DEBUG, "i2Start", *i2Start);
	printDebug(MID_DEBUG, "i2End",   *i2End);
	printDebug(MID_DEBUG, "i3Start", *i3Start);
	printDebug(MID_DEBUG, "i3End",   *i3End);

	printDebug(FULL_DEBUG, "OUT Grid::getGridIndex");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Myint64 Grid::getNumberOfGridPoint(Grid_type gridTypeIn, Point_type pointType)
{
	printDebug(FULL_DEBUG, "IN Grid::getNumberOfGridPoint");

	Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End, retVal ;
	getGridIndex(pointType, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End) ;
	Myint64 nGridPointLoc = (i1End-i1Start+1)*(i2End-i2Start+1)*(i3End-i3Start+1) ;

	Myint64 nGridPoint = 0 ;
	if ((gridTypeIn == GRID_GLOBAL) && (gridType == GRID_LOCAL))
	{
		MPI_Allreduce(&nGridPointLoc, &nGridPoint, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
	}
	else
	{
		nGridPoint = nGridPointLoc ;
	}

	printDebug(FULL_DEBUG, "OUT Grid::getNumberOfGridPoint");
	return(nGridPoint) ;
}

//-------------------------------------------------------------------------------------------------------

bool Grid::sameSize(const Grid& gridIn) const
{
	printDebug(FULL_DEBUG, "IN Grid::sameSize");

	bool retVal = true ;
	if (this->npoint != gridIn.npoint) retVal = false ;
	if (this->n1 != gridIn.n1) retVal = false ;
	if (this->n2 != gridIn.n2) retVal = false ;
	if (this->n3 != gridIn.n3) retVal = false ;

	printDebug(FULL_DEBUG, "OUT Grid::sameSize");
	return(retVal) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid::exchangeHalos(MPI_comm_mode_type commMode)
{
	printDebug(FULL_DEBUG, "IN Grid::exchangeHalos");

	if (commMode == MPI_COMM_MODE_SENDRECV)
	{
		printDebug(MID_DEBUG, "MPI_COMM_MODE_SENDRECV") ;

		// exchange the 6 halos
		if (exchangeHalo(MPI_COMM_MODE_SENDRECV, I1HALO1) != RTN_CODE_OK) return (RTN_CODE_KO) ;
		if (exchangeHalo(MPI_COMM_MODE_SENDRECV, I1HALO2) != RTN_CODE_OK) return (RTN_CODE_KO) ;
		if (exchangeHalo(MPI_COMM_MODE_SENDRECV, I2HALO1) != RTN_CODE_OK) return (RTN_CODE_KO) ;
		if (exchangeHalo(MPI_COMM_MODE_SENDRECV, I2HALO2) != RTN_CODE_OK) return (RTN_CODE_KO) ;
		if (exchangeHalo(MPI_COMM_MODE_SENDRECV, I3HALO1) != RTN_CODE_OK) return (RTN_CODE_KO) ;
		if (exchangeHalo(MPI_COMM_MODE_SENDRECV, I3HALO2) != RTN_CODE_OK) return (RTN_CODE_KO) ;

	}
	else
	{
		printError("IN Grid::exchangeHalos, invalid commMode", commMode) ;
		return(RTN_CODE_KO) ;
	}

	printDebug(FULL_DEBUG, "OUT Grid::exchangeHalos");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid::exchangeHalo(MPI_comm_mode_type commMode, Point_type pointType)
{
	printDebug(FULL_DEBUG, "IN Grid::exchangeHalo");

	if (commMode == MPI_COMM_MODE_SENDRECV)
	{
		printDebug(FULL_DEBUG, "MPI_COMM_MODE_SENDRECV") ;

		MPI_Status status ;
		Myint rankSend, rankRecv ;
		MPI_Datatype typeSend, typeRecv ;
		Myfloat *bufSend, *bufRecv ;

		if (pointType == I1HALO1)
		{
			printDebug(FULL_DEBUG, "I1HALO1") ;
			Myint64 i1Send = i1InnerEnd - haloWidth + 1 ;
			Myint64 i2Send = i2InnerStart ;
			Myint64 i3Send = i3InnerStart ;
			Myint64 i1Recv = i1Halo1Start ;
			Myint64 i2Recv = i2InnerStart ;
			Myint64 i3Recv = i3InnerStart ;
			Myint64 idxSend = i1Send+i2Send*n1+i3Send*n1*n2 ;
			Myint64 idxRecv = i1Recv+i2Recv*n1+i3Recv*n1*n2 ;
			bufSend = &(grid_3d[idxSend]) ;
			bufRecv = &(grid_3d[idxRecv]) ;
			typeSend = i1HaloDataType ;
			typeRecv = i1HaloDataType ;
			rankSend = i1ProcIdEnd ;
			rankRecv = i1ProcIdStart ;
		}

		else if (pointType == I1HALO2)
		{
			printDebug(FULL_DEBUG, "I1HALO2") ;
			Myint64 i1Send = i1InnerStart ;
			Myint64 i2Send = i2InnerStart ;
			Myint64 i3Send = i3InnerStart ;
			Myint64 i1Recv = i1Halo2Start ;
			Myint64 i2Recv = i2InnerStart ;
			Myint64 i3Recv = i3InnerStart ;
			Myint64 idxSend = i1Send+i2Send*n1+i3Send*n1*n2 ;
			Myint64 idxRecv = i1Recv+i2Recv*n1+i3Recv*n1*n2 ;
			bufSend = &(grid_3d[idxSend]) ;
			bufRecv = &(grid_3d[idxRecv]) ;
			typeSend = i1HaloDataType ;
			typeRecv = i1HaloDataType ;
			rankSend = i1ProcIdStart ;
			rankRecv = i1ProcIdEnd ;
		}

		else if (pointType == I2HALO1)
		{
			printDebug(FULL_DEBUG, "I2HALO1") ;
			Myint64 i1Send = i1InnerStart ;
			Myint64 i2Send = i2InnerEnd - haloWidth + 1 ;
			Myint64 i3Send = i3InnerStart ;
			Myint64 i1Recv = i1InnerStart ;
			Myint64 i2Recv = i2Halo1Start ;
			Myint64 i3Recv = i3InnerStart ;
			Myint64 idxSend = i1Send+i2Send*n1+i3Send*n1*n2 ;
			Myint64 idxRecv = i1Recv+i2Recv*n1+i3Recv*n1*n2 ;
			bufSend = &(grid_3d[idxSend]) ;
			bufRecv = &(grid_3d[idxRecv]) ;
			typeSend = i2HaloDataType ;
			typeRecv = i2HaloDataType ;
			rankSend = i2ProcIdEnd ;
			rankRecv = i2ProcIdStart ;
		}

		else if (pointType == I2HALO2)
		{
			printDebug(FULL_DEBUG, "I2HALO2") ;
			Myint64 i1Send = i1InnerStart ;
			Myint64 i2Send = i2InnerStart ;
			Myint64 i3Send = i3InnerStart ;
			Myint64 i1Recv = i1InnerStart ;
			Myint64 i2Recv = i2Halo2Start ;
			Myint64 i3Recv = i3InnerStart ;
			Myint64 idxSend = i1Send+i2Send*n1+i3Send*n1*n2 ;
			Myint64 idxRecv = i1Recv+i2Recv*n1+i3Recv*n1*n2 ;
			bufSend = &(grid_3d[idxSend]) ;
			bufRecv = &(grid_3d[idxRecv]) ;
			typeSend = i2HaloDataType ;
			typeRecv = i2HaloDataType ;
			rankSend = i2ProcIdStart ;
			rankRecv = i2ProcIdEnd ;
		}

		else if (pointType == I3HALO1)
		{
			printDebug(FULL_DEBUG, "I3HALO1") ;
			Myint64 i1Send = i1InnerStart ;
			Myint64 i2Send = i2InnerStart ;
			Myint64 i3Send = i3InnerEnd - haloWidth + 1 ;
			Myint64 i1Recv = i1InnerStart ;
			Myint64 i2Recv = i2InnerStart ;
			Myint64 i3Recv = i3Halo1Start ;
			Myint64 idxSend = i1Send+i2Send*n1+i3Send*n1*n2 ;
			Myint64 idxRecv = i1Recv+i2Recv*n1+i3Recv*n1*n2 ;
			bufSend = &(grid_3d[idxSend]) ;
			bufRecv = &(grid_3d[idxRecv]) ;
			typeSend = i3HaloDataType ;
			typeRecv = i3HaloDataType ;
			rankSend = i3ProcIdEnd ;
			rankRecv = i3ProcIdStart ;
		}

		else if (pointType == I3HALO2)
		{
			printDebug(FULL_DEBUG, "I3HALO2") ;
			Myint64 i1Send = i1InnerStart ;
			Myint64 i2Send = i2InnerStart ;
			Myint64 i3Send = i3InnerStart ;
			Myint64 i1Recv = i1InnerStart ;
			Myint64 i2Recv = i2InnerStart ;
			Myint64 i3Recv = i3Halo2Start ;
			Myint64 idxSend = i1Send+i2Send*n1+i3Send*n1*n2 ;
			Myint64 idxRecv = i1Recv+i2Recv*n1+i3Recv*n1*n2 ;
			bufSend = &(grid_3d[idxSend]) ;
			bufRecv = &(grid_3d[idxRecv]) ;
			typeSend = i3HaloDataType ;
			typeRecv = i3HaloDataType ;
			rankSend = i3ProcIdStart ;
			rankRecv = i3ProcIdEnd ;
		}

		else
		{
			printError("IN Grid::exchangeHalo, invalid pointType", pointType) ;
			return(RTN_CODE_KO) ;
		}

		// call MPI_Sendrecv
		printDebug(FULL_DEBUG, "MPI_Sendrecv", rankSend, rankRecv) ;
		MPI_Sendrecv(bufSend, 1, typeSend, rankSend, 0,
				bufRecv, 1, typeRecv, rankRecv, 0,
				MPI_COMM_WORLD, &status);
	}
	else
	{
		printError("IN Grid::exchangeHalo, invalid commMode", commMode) ;
		return(RTN_CODE_KO) ;
	}

	printDebug(FULL_DEBUG, "OUT Grid::exchangeHalo");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid::applyBoundaryCondition(BoundCond_type boundCondType)
{
	printDebug(FULL_DEBUG, "IN Grid::applyBoundaryCondition");

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
		printError("IN Grid::applyBoundaryCondition, invalid boundCondType", boundCondType) ;
		return(RTN_CODE_KO) ;
	}

	printDebug(FULL_DEBUG, "OUT Grid::applyBoundaryCondition");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid::FD_D2_N1(Point_type pType, const Grid& Wgrid, Myint fdOrder)
{
	printDebug(MID_DEBUG, "IN Grid::FD_D2_N1");

	// check grids are same size
	if (this->sameSize(Wgrid) != true)
	{
		printError("Grid::FD_D2_N1, grids have not same size") ;
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

	// compute FD along N1
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
					w[i1+i2*n1+i3*n1*n2] =
							FD_D2_O2_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
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
					w[i1+i2*n1+i3*n1*n2] =
							FD_D2_O4_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
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
					w[i1+i2*n1+i3*n1*n2] =
							FD_D2_O8_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
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
					w[i1+i2*n1+i3*n1*n2] =
							FD_D2_O12_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
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
					w[i1+i2*n1+i3*n1*n2] =
							FD_D2_O16_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
				}

			}
		}
	}

	printDebug(MID_DEBUG, "OUT Grid::FD_D2_N1");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid::FD_D2_N2(Point_type pType, const Grid& Wgrid, Myint fdOrder)
{
	printDebug(MID_DEBUG, "IN Grid::FD_D2_N2");

	// check grids are same size
	if (this->sameSize(Wgrid) != true)
	{
		printError("Grid::FD_D2_N2, grids have not same size") ;
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
					w[i1+i2*n1+i3*n1*n2] =
							FD_D2_O2_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
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
					w[i1+i2*n1+i3*n1*n2] =
							FD_D2_O4_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
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
					w[i1+i2*n1+i3*n1*n2] =
							FD_D2_O8_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
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
					w[i1+i2*n1+i3*n1*n2] =
							FD_D2_O12_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
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
					w[i1+i2*n1+i3*n1*n2] =
							FD_D2_O16_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
				}

			}
		}
	}

	printDebug(MID_DEBUG, "OUT Grid::FD_D2_N2");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid::FD_D2_N3(Point_type pType, const Grid& Wgrid, Myint fdOrder)
{
	printDebug(MID_DEBUG, "IN Grid::FD_D2_N3");

	// check grids are same size
	if (this->sameSize(Wgrid) != true)
	{
		printError("Grid::FD_D2_N3, grids have not same size") ;
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
					w[i1+i2*n1+i3*n1*n2] =
							FD_D2_O2_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
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
					w[i1+i2*n1+i3*n1*n2] =
							FD_D2_O4_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
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
					w[i1+i2*n1+i3*n1*n2] =
							FD_D2_O8_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
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
					w[i1+i2*n1+i3*n1*n2] =
							FD_D2_O12_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
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
					w[i1+i2*n1+i3*n1*n2] =
							FD_D2_O16_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
				}

			}
		}
	}

	printDebug(MID_DEBUG, "OUT Grid::FD_D2_N3");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid::FD_LAPLACIAN(Point_type pType, const Grid& Wgrid, Myint fdOrder)
{
	printDebug(MID_DEBUG, "IN Grid::FD_LAPLACIAN");

	// check grids are same size
	if (this->sameSize(Wgrid) != true)
	{
		printError("Grid::FD_LAPLACIAN, grids have not same size") ;
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
#pragma omp simd
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
#pragma omp simd
					for (Myint64 i1 = i1Start; i1<= i1End; i1++)
					{
						w[i1+i2*n1+i3*n1*n2] =
								FD_D2_O4_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
								+ FD_D2_O4_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
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
						w[i1+i2*n1+i3*n1*n2] =
								FD_D2_O8_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
								+ FD_D2_O8_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
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
						w[i1+i2*n1+i3*n1*n2] =
								FD_D2_O12_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
								+ FD_D2_O12_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
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
#pragma omp simd
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
#pragma omp simd
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
						w[i1+i2*n1+i3*n1*n2] =
								FD_D2_O8_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
								+ FD_D2_O8_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
								+ FD_D2_O8_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
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
						w[i1+i2*n1+i3*n1*n2] =
								FD_D2_O12_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
								+ FD_D2_O12_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
								+ FD_D2_O12_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
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
						w[i1+i2*n1+i3*n1*n2] =
								FD_D2_O16_N1(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
								+ FD_D2_O16_N2(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3)
								+ FD_D2_O16_N3(u, i1, i2, i3, inv2_d1, inv2_d2, inv2_d3, n1, n2, n3) ;
					}

				}
			}
		}
	}

	printDebug(MID_DEBUG, "OUT Grid::FD_LAPLACIAN");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid::updatePressure(Point_type pType, const Grid& prcGrid,
		const Grid& coefGrid, const Grid& laplaGrid)
{
	printDebug(MID_DEBUG, "IN Grid::updatePressure");

	Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
	getGridIndex(INNER_POINTS, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End) ;

	Myfloat * prn   = this->grid_3d ;
	Myfloat * prc   = prcGrid.grid_3d ;
	Myfloat * lapla = laplaGrid.grid_3d ;
	Myfloat * coef  = coefGrid.grid_3d ;

#pragma omp parallel for collapse(2)
	for (Myint64 i3 = i3Start; i3<= i3End; i3++)
	{
		for (Myint64 i2 = i2Start; i2<= i2End; i2++)
		{
#pragma omp simd
			for (Myint64 i1 = i1Start; i1<= i1End; i1++)
			{
				prn[i1+i2*n1+i3*n1*n2] = TWO * prc[i1+i2*n1+i3*n1*n2] - prn[i1+i2*n1+i3*n1*n2] +
						coef[i1+i2*n1+i3*n1*n2] * lapla[i1+i2*n1+i3*n1*n2] ;
			}
		}
	}

	printDebug(MID_DEBUG, "OUT Grid::updatePressure");
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid::computePressureWithFD(Grid& prcGridIn, Grid& coefGridIn, Myint fdOrder)
{

	printDebug(FULL_DEBUG, "In Grid::computePressureWithFD") ;

	// check grids are same size
	if (this->sameSize(prcGridIn) != true)
	{
		printError("In Grid::computePressureWithFD, grids have not same size") ;
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

	printDebug(FULL_DEBUG, "Out Grid::computePressureWithFD") ;
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid::sendWithMPI(Myint64 nGridPoint, Myint procDestId)
{

	printDebug(FULL_DEBUG, "In Grid::sendWithMPI") ;

	MPI_Send(grid_3d, nGridPoint, MPI_MYFLOAT, procDestId, 0, MPI_COMM_WORLD) ;

	printDebug(FULL_DEBUG, "Out Grid::sendWithMPI") ;
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid::recvWithMPI(Myint64 nGridPoint, Myint procSrcId)
{

	printDebug(FULL_DEBUG, "In Grid::recvWithMPI") ;

	MPI_Status status ;
	MPI_Recv(grid_3d, nGridPoint, MPI_MYFLOAT, procSrcId, 0, MPI_COMM_WORLD, &status) ;

	if (status.MPI_ERROR != MPI_SUCCESS)
	{
		//printError("MPI ERROR", status.MPI_ERROR) ;
	}

	printDebug(FULL_DEBUG, "Out Grid::recvWithMPI") ;
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Grid::sendRecvWithMPI(const Grid& gridDest, Myint idSend, Myint idRecv, Myint64 nGridPoint)
{

	printDebug(FULL_DEBUG, "In Grid::sendRecvWithMPI") ;

	Myfloat *bufSend = grid_3d ;
	Myfloat *bufRecv = gridDest.grid_3d ;

	MPI_Status status ;
	MPI_Sendrecv(bufSend, nGridPoint, MPI_MYFLOAT, idSend, 0,
			bufRecv, nGridPoint, MPI_MYFLOAT, idRecv, 0,
			MPI_COMM_WORLD, &status) ;
	if (status.MPI_ERROR != MPI_SUCCESS)
	{
		//printError("MPI ERROR", status.MPI_ERROR) ;
	}

	printDebug(FULL_DEBUG, "Out Grid::sendRecvWithMPI") ;
	return(RTN_CODE_OK) ;
}

} // namespace hpcscan
