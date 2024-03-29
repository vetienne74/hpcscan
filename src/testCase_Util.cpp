
//-------------------------------------------------------------------------------------------------------
// test cases for internal use by developers
// these cases are used to check the good behavior of internal funcions
// and not for the purpose of performance measurements
//-------------------------------------------------------------------------------------------------------

#include "testCase_Util.h"

#include <algorithm>
#include <cassert>
#include <cfloat>

#include <mpi.h>
#include <omp.h>

#include "config.h"
#include "constant.h"
#include "data_Acquisition_Factory.h"
#include "fdm.h"
#include "global.h"
#include "grid_Factory.h"
#include "output_report.h"

namespace hpcscan {

TestCase_Util::TestCase_Util(void)
		{
	testCaseName    = "Util" ;
	testCaseVersion = "Standard implementation" ;
		}

//-------------------------------------------------------------------------------------------------------

Rtn_code TestCase_Util::run(void)
{

	printDebug(MID_DEBUG, "In TestCase_Util::run") ;

	if (this->initialize() == RTN_CODE_KO) return (RTN_CODE_KO) ;

	const string gridMode = Config::Instance()->testMode ;

	// display grid info
	{
		auto gridLoc2 = Grid_Factory::create(gridMode, GRID_LOCAL) ;
		Grid &gridLoc = *gridLoc2 ;
		gridLoc.initializeGrid() ;
		gridLoc.info() ;
	}

	// create a block for each test
	{
		//------------------------------------
		// check 2 global grids have same size
		//------------------------------------
		print_blank() ;
		string caseName = testCaseName + "GridSameSize" ;
		printInfo(MASTER, " * Case", caseName) ;

		auto gridGlob12 = Grid_Factory::create(gridMode, GRID_GLOBAL) ;
		Grid &gridGlob1 = *gridGlob12 ;
		gridGlob1.initializeGrid() ;

		auto gridGlob22 = Grid_Factory::create(gridMode, GRID_GLOBAL) ;
		Grid &gridGlob2 = *gridGlob22 ;
		gridGlob2.initializeGrid() ;

		bool sameSize = gridGlob1.sameSize(gridGlob2) ;
		checkBoolDiff(sameSize, true) ;
	}

	{
		//-----------------------------------------
		// check size between global and local grid
		// are different when nMpiProc > 1
		//-----------------------------------------
		print_blank() ;
		string caseName = testCaseName + "GridSizeGlobalLocal" ;
		printInfo(MASTER, " * Case", caseName) ;

		auto gridGlob12 = Grid_Factory::create(gridMode, GRID_GLOBAL) ;
		Grid &gridGlob1 = *gridGlob12 ;
		gridGlob1.initializeGrid() ;

		auto gridLoc12 = Grid_Factory::create(gridMode, GRID_LOCAL) ;
		Grid &gridLoc1 = *gridLoc12 ;
		gridLoc1.initializeGrid() ;

		bool sameSize2 = gridGlob1.sameSize(gridLoc1) ;
		if (nMpiProc == 1)
		{
			checkBoolDiff(sameSize2, true) ;
		}
		else
		{
			checkBoolDiff(sameSize2, false) ;
		}
	}

	{
		//-----------------------------------------------------------
		// check number of Inner points between global and local grid
		// are different when nMpiProc > 1
		//------------------------------------------------------------
		print_blank() ;
		string caseName = testCaseName + "GridSizeGlobalLocal" ;
		printInfo(MASTER, " * Case", caseName) ;

		auto gridGlob12 = Grid_Factory::create(gridMode, GRID_GLOBAL) ;
		Grid &gridGlob1 = *gridGlob12 ;
		gridGlob1.initializeGrid() ;

		auto gridLoc12 = Grid_Factory::create(gridMode, GRID_LOCAL) ;
		Grid &gridLoc1 = *gridLoc12 ;
		gridLoc1.initializeGrid() ;

		Myint64 localSize = gridLoc1.getNumberOfGridPoint(GRID_LOCAL, INNER_POINTS) ;
		Myint64 globalSize = gridGlob1.getNumberOfGridPoint(GRID_GLOBAL, INNER_POINTS) ;
		Myint64 sumLocalSize ;
		MPI_Reduce(&localSize, &sumLocalSize, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
		if (checkIntegerDiff(globalSize, sumLocalSize) != RTN_CODE_OK)
		{
			printInfo(MASTER, "globalSize", globalSize) ;
			printInfo(MASTER, "sumSize", sumLocalSize) ;
		}
	}

	{
		//---------------------------------------------
		// check number of Halos comm points is correct
		//---------------------------------------------
		print_blank() ;
		string caseName = testCaseName + "NumberOfHaloCommPoint" ;
		printInfo(MASTER, " * Case", caseName) ;

		auto gridLoc12 = Grid_Factory::create(gridMode, GRID_LOCAL) ;
		Grid &gridLoc1 = *gridLoc12 ;
		gridLoc1.initializeGrid() ;

		Myint64 nGridPointHaloGlob = gridLoc1.getNumberOfGridPointCommHalo(GRID_GLOBAL) ;

		// compute number of halo comm point
		auto dim     = gridLoc1.dim ;
		auto nsub1   = Config::Instance()->nsub1 ;
		auto nsub2   = Config::Instance()->nsub2 ;
		auto nsub3   = Config::Instance()->nsub3 ;
		auto n1      = Config::Instance()->n1 ;
		auto n2      = Config::Instance()->n2 ;
		auto n3      = Config::Instance()->n3 ;
		auto fdOrder = Config::Instance()->fdOrder ;

		Myint64 nGridPointHaloGlob2 = 0 ;

		if (dim ==DIM1)
		{
			nGridPointHaloGlob2 = (nsub1 - 1) * fdOrder ;
		}
		else if (dim ==DIM2)
		{
			nGridPointHaloGlob2 = (nsub1 - 1) * n2*fdOrder + (nsub2 - 1) * n1*fdOrder ;
		}
		else if (dim ==DIM3)
		{
			nGridPointHaloGlob2 = (nsub1 - 1) * n2*n3*fdOrder + (nsub2 - 1) * n1*n3*fdOrder + (nsub3 - 1) * n1*n2*fdOrder ;
		}
		if (checkIntegerDiff(nGridPointHaloGlob, nGridPointHaloGlob2) != RTN_CODE_OK)
		{
			printInfo(MASTER, "nGridPointHaloGlob", nGridPointHaloGlob) ;
			printInfo(MASTER, "nGridPointHaloGlob2", nGridPointHaloGlob2) ;
		}
	}

	{
		//--------------------------------
		// check number of proc neighbours
		//--------------------------------
		print_blank() ;
		string caseName = testCaseName + "NumberOfNeighbours" ;
		printInfo(MASTER, " * Case", caseName) ;

		auto gridLoc2 = Grid_Factory::create(gridMode, GRID_LOCAL) ;
		Grid &gridLoc = *gridLoc2 ;
		gridLoc.initializeGrid() ;

		Myint nNeighLoc  = gridLoc.getNumberOfNeighbour() ;
		Myint nNeighGlobMax = 0 ;
		Myint nNeighGlobMin = 0 ;

		// check min and max number of neighbours
		MPI_Reduce(&nNeighLoc, &nNeighGlobMin, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
		MPI_Reduce(&nNeighLoc, &nNeighGlobMax, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

		if (nMpiProc == 1)
		{
			bool checkTest = ((nNeighGlobMin == 0) && (nNeighGlobMax == 0)) ;
			if (checkBoolDiff(checkTest, true) != RTN_CODE_OK)
			{
				printInfo(MASTER, " nNeighGlobMin", nNeighGlobMin) ;
				printInfo(MASTER, " nNeighGlobMax", nNeighGlobMax) ;
			}
		}
		else
		{
			bool checkTest = ((nNeighGlobMin > 0) && (nNeighGlobMax <= gridLoc.dim*2)) ;
			if (checkBoolDiff(checkTest, true) != RTN_CODE_OK)
			{
				printInfo(MASTER, " nNeighGlobMin", nNeighGlobMin) ;
				printInfo(MASTER, " nNeighGlobMax", nNeighGlobMax) ;
			}
		}
	}

	auto gridLoc2 = Grid_Factory::create(gridMode, GRID_LOCAL) ;
	Grid &gridLoc = *gridLoc2 ;
	gridLoc.initializeGrid() ;

	auto gridRef2 = Grid_Factory::create(gridMode, GRID_LOCAL) ;
	Grid &gridRef = *gridRef2 ;
	gridRef.initializeGrid() ;

	{
		//--------------------------------
		// check fill grid with constant
		//--------------------------------
		print_blank() ;
		string caseName = testCaseName + "FillWithConstant" ;
		printInfo(MASTER, " * Case", caseName) ;

		Myfloat valConst = 999.0 ;
		gridLoc.fill(ALL_POINTS, valConst) ;
		gridLoc.write(caseName+"Loc") ;

		Myfloat gridLocMin = gridLoc.getMin(ALL_POINTS) ;
		Myfloat gridLocGlobMin = 0 ;
		MPI_Reduce(&gridLocMin, &gridLocGlobMin, 1, MPI_MYFLOAT, MPI_MIN, 0, MPI_COMM_WORLD);
		bool checkMin = relErr(gridLocGlobMin, valConst) < MAX_ERR_FLOAT ;

		Myfloat gridLocMax = gridLoc.getMax(ALL_POINTS) ;
		Myfloat gridLocGlobMax = 0 ;
		MPI_Reduce(&gridLocMax, &gridLocGlobMax, 1, MPI_MYFLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
		bool checkMax = relErr(gridLocGlobMax, valConst) < MAX_ERR_FLOAT ;

		checkBoolDiff((checkMin && checkMax), true) ;

	}

	{
		//--------------------------------
		// check fill grid with sine
		//--------------------------------
		print_blank() ;
		string caseName = testCaseName + "FillWithSine" ;
		printInfo(MASTER, " * Case", caseName) ;

		Myfloat valConst = 999.0 ;
		gridLoc.fill(ALL_POINTS, FUNC_SINE, FUNC_SINE, FUNC_SINE, 1.0, 1.0, 1.0, valConst) ;
		gridLoc.write(caseName+"Loc") ;

		Myfloat gridLocMin = gridLoc.getMin(ALL_POINTS) ;
		Myfloat gridLocGlobMin = 0 ;
		MPI_Reduce(&gridLocMin, &gridLocGlobMin, 1, MPI_MYFLOAT, MPI_MIN, 0, MPI_COMM_WORLD);
		bool checkMin = relErr(gridLocGlobMin, -valConst) < 0.1 ;

		Myfloat gridLocMax = gridLoc.getMax(ALL_POINTS) ;
		Myfloat gridLocGlobMax = 0 ;
		MPI_Reduce(&gridLocMax, &gridLocGlobMax, 1, MPI_MYFLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
		bool checkMax = relErr(gridLocGlobMax, valConst) < 0.1 ;

		checkBoolDiff((checkMin && checkMax), true) ;
	}

	{
		//---------------------
		// check UpdateI1HALO1
		//---------------------
		print_blank() ;
		string caseName = testCaseName + "UpdateI1HALO1" ;
		printInfo(MASTER, " * Case", caseName) ;

		// fill gridLoc with proc rank
		gridLoc.fill(ALL_POINTS, myMpiRank) ;
		gridLoc.exchangeHalo(MPI_COMM_MODE_SENDRECV, I1HALO1) ;
		gridLoc.write(caseName+"Loc") ;

		// fill gridRef with expected values
		gridRef.fill(ALL_POINTS, myMpiRank) ;
		// set values in halo
		Myint rankRecv = gridLoc.getNeighbourProc(I1HALO1) ;
		if (rankRecv != MPI_PROC_NULL) gridRef.fill(I1HALO1, rankRecv) ;
		gridRef.write(caseName+"Ref") ;

		checkAllProcGridL1Err(ALL_POINTS, gridLoc, gridRef, MAX_ERR_FLOAT) ;
	}

	{
		//---------------------
		// check UpdateI1HALO2
		//---------------------
		print_blank() ;
		string caseName = testCaseName + "UpdateI1HALO2" ;
		printInfo(MASTER, " * Case", caseName) ;

		// fill gridLoc with proc rank
		gridLoc.fill(ALL_POINTS, myMpiRank) ;
		gridLoc.exchangeHalo(MPI_COMM_MODE_SENDRECV, I1HALO2) ;
		gridLoc.write(caseName+"Loc") ;

		// fill gridRef with expected values
		gridRef.fill(ALL_POINTS, myMpiRank) ;
		// set values in halo
		Myint rankRecv = gridLoc.getNeighbourProc(I1HALO2) ;
		if (rankRecv != MPI_PROC_NULL) gridRef.fill(I1HALO2, rankRecv) ;
		gridRef.write(caseName+"Ref") ;

		checkAllProcGridL1Err(ALL_POINTS, gridLoc, gridRef, MAX_ERR_FLOAT) ;
	}

	{
		//---------------------
		// check UpdateI2HALO1
		//---------------------
		print_blank() ;
		string caseName = testCaseName + "UpdateI2HALO1" ;
		printInfo(MASTER, " * Case", caseName) ;

		// fill gridLoc with proc rank
		gridLoc.fill(ALL_POINTS, myMpiRank) ;
		gridLoc.exchangeHalo(MPI_COMM_MODE_SENDRECV, I2HALO1) ;
		gridLoc.write(caseName+"Loc") ;

		// fill gridRef with expected values
		gridRef.fill(ALL_POINTS, myMpiRank) ;
		// set values in halo
		Myint rankRecv = gridLoc.getNeighbourProc(I2HALO1) ;
		if (rankRecv != MPI_PROC_NULL) gridRef.fill(I2HALO1, rankRecv) ;
		gridRef.write(caseName+"Ref") ;

		checkAllProcGridL1Err(ALL_POINTS, gridLoc, gridRef, MAX_ERR_FLOAT) ;
	}

	{
		//---------------------
		// check UpdateI2HALO2
		//---------------------
		print_blank() ;
		string caseName = testCaseName + "UpdateI2HALO2" ;
		printInfo(MASTER, " * Case", caseName) ;

		// fill gridLoc with proc rank
		gridLoc.fill(ALL_POINTS, myMpiRank) ;
		gridLoc.exchangeHalo(MPI_COMM_MODE_SENDRECV, I2HALO2) ;
		gridLoc.write(caseName+"Loc") ;

		// fill gridRef with expected values
		gridRef.fill(ALL_POINTS, myMpiRank) ;
		// set values in halo
		Myint rankRecv = gridLoc.getNeighbourProc(I2HALO2) ;
		if (rankRecv != MPI_PROC_NULL) gridRef.fill(I2HALO2, rankRecv) ;
		gridRef.write(caseName+"Ref") ;

		checkAllProcGridL1Err(ALL_POINTS, gridLoc, gridRef, MAX_ERR_FLOAT) ;
	}

	{
		//---------------------
		// check UpdateI3HALO1
		//---------------------
		print_blank() ;
		string caseName = testCaseName + "UpdateI3HALO1" ;
		printInfo(MASTER, " * Case", caseName) ;

		// fill gridLoc with proc rank
		gridLoc.fill(ALL_POINTS, myMpiRank) ;
		gridLoc.exchangeHalo(MPI_COMM_MODE_SENDRECV, I3HALO1) ;
		gridLoc.write(caseName+"Loc") ;

		// fill gridRef with expected values
		gridRef.fill(ALL_POINTS, myMpiRank) ;
		// set values in halo
		Myint rankRecv = gridLoc.getNeighbourProc(I3HALO1) ;
		if (rankRecv != MPI_PROC_NULL) gridRef.fill(I3HALO1, rankRecv) ;
		gridRef.write(caseName+"Ref") ;

		checkAllProcGridL1Err(ALL_POINTS, gridLoc, gridRef, MAX_ERR_FLOAT) ;
	}

	{
		//---------------------
		// check UpdateI3HALO2
		//---------------------
		print_blank() ;
		string caseName = testCaseName + "UpdateI3HALO2" ;
		printInfo(MASTER, " * Case", caseName) ;

		// fill gridLoc with proc rank
		gridLoc.fill(ALL_POINTS, myMpiRank) ;
		gridLoc.exchangeHalo(MPI_COMM_MODE_SENDRECV, I3HALO2) ;
		gridLoc.write(caseName+"Loc") ;

		// fill gridRef with expected values
		gridRef.fill(ALL_POINTS, myMpiRank) ;
		// set values in halo
		Myint rankRecv = gridLoc.getNeighbourProc(I3HALO2) ;
		if (rankRecv != MPI_PROC_NULL) gridRef.fill(I3HALO2, rankRecv) ;
		gridRef.write(caseName+"Ref") ;

		checkAllProcGridL1Err(ALL_POINTS, gridLoc, gridRef, MAX_ERR_FLOAT) ;
	}

	const Myfloat64 a1 = Config::Instance()->param1 ;
	const Myfloat64 a2 = Config::Instance()->param2 ;
	const Myfloat64 a3 = Config::Instance()->param3 ;
	const Myfloat64 a4 = Config::Instance()->param4 ;

	{
		//---------------------
		// check exchangeHalos
		//---------------------
		print_blank() ;
		string caseName = testCaseName + "ExchangeHalos" ;
		printInfo(MASTER, " * Case", caseName) ;

		// fill gridLoc with linear function
		// note: any function can be used (linear is easy to debug)
		gridLoc.fill(ALL_POINTS, FUNC_LINEAR, FUNC_LINEAR, FUNC_LINEAR, a1, a2, a3, a4) ;

		// set all halos to zeros
		if (gridLoc.getNeighbourProc(I1HALO1) != MPI_PROC_NULL) gridLoc.fill(I1HALO1, -1) ;
		if (gridLoc.getNeighbourProc(I1HALO2) != MPI_PROC_NULL) gridLoc.fill(I1HALO2, -1) ;
		if (gridLoc.getNeighbourProc(I2HALO1) != MPI_PROC_NULL) gridLoc.fill(I2HALO1, -1) ;
		if (gridLoc.getNeighbourProc(I2HALO2) != MPI_PROC_NULL) gridLoc.fill(I2HALO2, -1) ;
		if (gridLoc.getNeighbourProc(I3HALO1) != MPI_PROC_NULL) gridLoc.fill(I3HALO1, -1) ;
		if (gridLoc.getNeighbourProc(I3HALO2) != MPI_PROC_NULL) gridLoc.fill(I3HALO2, -1) ;

		// update all halos
		gridLoc.exchangeHalos(MPI_COMM_MODE_SENDRECV) ;
		gridLoc.write(caseName+"Loc") ;

		// fill gridRef with expected values
		gridRef.fill(ALL_POINTS, FUNC_LINEAR, FUNC_LINEAR, FUNC_LINEAR, a1, a2, a3, a4) ;
		gridRef.write(caseName+"Ref") ;

		checkAllProcGridL1Err(ALL_POINTS, gridLoc, gridRef, MAX_ERR_FLOAT) ;
	}

	{
		//-----------------------
		// check grid coordinates
		//-----------------------
		print_blank() ;
		string caseName = testCaseName + "GridCoord" ;
		printInfo(MASTER, " * Case", caseName) ;

		Myint nFailed = 0 ;

		Myfloat minCoord1 = gridLoc.getMinCoord(AXIS1) ;
		Myfloat minCoord1Glob ;
		Myfloat minCoord1Ref = 0.0 ;
		MPI_Reduce(&minCoord1, &minCoord1Glob, 1, MPI_MYFLOAT, MPI_MIN, 0, MPI_COMM_WORLD);
		if (relErr(minCoord1Glob, minCoord1Ref) > MAX_ERR_FLOAT) nFailed++ ;

		Myfloat maxCoord1 = gridLoc.getMaxCoord(AXIS1) ;
		Myfloat maxCoord1Glob ;
		Myfloat maxCoord1Ref = minCoord1Ref + (Config::Instance()->n1-1) *gridLoc.d1 ;
		MPI_Reduce(&maxCoord1, &maxCoord1Glob, 1, MPI_MYFLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
		if (relErr(maxCoord1Glob, maxCoord1Ref) > MAX_ERR_FLOAT) nFailed++ ;

		if (gridLoc.dim >= DIM2)
		{
			Myfloat minCoord2 = gridLoc.getMinCoord(AXIS2);
			Myfloat minCoord2Glob ;
			Myfloat minCoord2Ref = 0.0 ;
			MPI_Reduce(&minCoord2, &minCoord2Glob, 1, MPI_MYFLOAT, MPI_MIN, 0, MPI_COMM_WORLD);
			if (relErr(minCoord2Glob, minCoord2Ref) > MAX_ERR_FLOAT) nFailed++ ;

			Myfloat maxCoord2 = gridLoc.getMaxCoord(AXIS2) ;
			Myfloat maxCoord2Glob ;
			Myfloat maxCoord2Ref = minCoord2Ref + (Config::Instance()->n2-1) *gridLoc.d2 ;
			MPI_Reduce(&maxCoord2, &maxCoord2Glob, 1, MPI_MYFLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
			if (relErr(maxCoord2Glob, maxCoord2Ref) > MAX_ERR_FLOAT) nFailed++ ;
		}

		if (gridLoc.dim >= DIM3)
		{
			Myfloat minCoord3 = gridLoc.getMinCoord(AXIS3);
			Myfloat minCoord3Glob ;
			Myfloat minCoord3Ref = 0.0 ;
			MPI_Reduce(&minCoord3, &minCoord3Glob, 1, MPI_MYFLOAT, MPI_MIN, 0, MPI_COMM_WORLD);
			if (relErr(minCoord3Glob, minCoord3Ref) > MAX_ERR_FLOAT) nFailed++ ;

			Myfloat maxCoord3 = gridLoc.getMaxCoord(AXIS3) ;
			Myfloat maxCoord3Glob ;
			Myfloat maxCoord3Ref = minCoord3Ref + (Config::Instance()->n3-1) *gridLoc.d3 ;
			MPI_Reduce(&maxCoord3, &maxCoord3Glob, 1, MPI_MYFLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
			if (relErr(maxCoord3Glob, maxCoord3Ref) > MAX_ERR_FLOAT) nFailed++ ;
		}

		checkIntegerDiff(nFailed, 0) ;
	}


	{
		//------------------------------------
		// check grid coordinates in unit grid
		//------------------------------------
		print_blank() ;
		string caseName = testCaseName + "UnitGridCoord" ;
		printInfo(MASTER, " * Case", caseName) ;

		Myint nFailed = 0 ;

		// define unit grid
		gridLoc.defineUnitGrid() ;

		Myfloat minCoord1 = gridLoc.getMinCoord(AXIS1) ;
		Myfloat minCoord1Glob ;
		Myfloat minCoord1Ref = 0.0 ;
		MPI_Reduce(&minCoord1, &minCoord1Glob, 1, MPI_MYFLOAT, MPI_MIN, 0, MPI_COMM_WORLD);
		if (relErr(minCoord1Glob, minCoord1Ref) > MAX_ERR_FLOAT) nFailed++ ;

		Myfloat maxCoord1 = gridLoc.getMaxCoord(AXIS1) ;
		Myfloat maxCoord1Glob ;
		Myfloat maxCoord1Ref = 1.0 ;
		MPI_Reduce(&maxCoord1, &maxCoord1Glob, 1, MPI_MYFLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
		if (relErr(maxCoord1Glob, maxCoord1Ref) > MAX_ERR_FLOAT) nFailed++ ;

		if (gridLoc.dim >= DIM2)
		{
			Myfloat minCoord2 = gridLoc.getMinCoord(AXIS2);
			Myfloat minCoord2Glob ;
			Myfloat minCoord2Ref = 0.0 ;
			MPI_Reduce(&minCoord2, &minCoord2Glob, 1, MPI_MYFLOAT, MPI_MIN, 0, MPI_COMM_WORLD);
			if (relErr(minCoord2Glob, minCoord2Ref) > MAX_ERR_FLOAT) nFailed++ ;

			Myfloat maxCoord2 = gridLoc.getMaxCoord(AXIS2) ;
			Myfloat maxCoord2Glob ;
			Myfloat maxCoord2Ref = 1.0 ;
			MPI_Reduce(&maxCoord2, &maxCoord2Glob, 1, MPI_MYFLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
			if (relErr(maxCoord2Glob, maxCoord2Ref) > MAX_ERR_FLOAT) nFailed++ ;
		}

		if (gridLoc.dim >= DIM3)
		{
			Myfloat minCoord3 = gridLoc.getMinCoord(AXIS3);
			Myfloat minCoord3Glob ;
			Myfloat minCoord3Ref = 0.0 ;
			MPI_Reduce(&minCoord3, &minCoord3Glob, 1, MPI_MYFLOAT, MPI_MIN, 0, MPI_COMM_WORLD);
			if (relErr(minCoord3Glob, minCoord3Ref) > MAX_ERR_FLOAT) nFailed++ ;

			Myfloat maxCoord3 = gridLoc.getMaxCoord(AXIS3) ;
			Myfloat maxCoord3Glob ;
			Myfloat maxCoord3Ref = 1.0 ;
			MPI_Reduce(&maxCoord3, &maxCoord3Glob, 1, MPI_MYFLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
			if (relErr(maxCoord3Glob, maxCoord3Ref) > MAX_ERR_FLOAT) nFailed++ ;
		}

		checkIntegerDiff(nFailed, 0) ;
	}

	{
		//---------------------------------
		// check exchangeHalos in unit grid
		//---------------------------------
		print_blank() ;
		string caseName = testCaseName + "UnitGridExchangeHalos" ;
		printInfo(MASTER, " * Case", caseName) ;

		// set unit grid
		gridLoc.defineUnitGrid() ;
		gridRef.defineUnitGrid() ;

		// fill gridLoc with linear function
		// note: any function can be used (linear is easy to debug)
		gridLoc.fill(ALL_POINTS, FUNC_LINEAR, FUNC_LINEAR, FUNC_LINEAR, a1, a2, a3, a4) ;

		// set all communication halos to -1
		if (gridLoc.getNeighbourProc(I1HALO1) != MPI_PROC_NULL) gridLoc.fill(I1HALO1, -1) ;
		if (gridLoc.getNeighbourProc(I1HALO2) != MPI_PROC_NULL) gridLoc.fill(I1HALO2, -1) ;
		if (gridLoc.getNeighbourProc(I2HALO1) != MPI_PROC_NULL) gridLoc.fill(I2HALO1, -1) ;
		if (gridLoc.getNeighbourProc(I2HALO2) != MPI_PROC_NULL) gridLoc.fill(I2HALO2, -1) ;
		if (gridLoc.getNeighbourProc(I3HALO1) != MPI_PROC_NULL) gridLoc.fill(I3HALO1, -1) ;
		if (gridLoc.getNeighbourProc(I3HALO2) != MPI_PROC_NULL) gridLoc.fill(I3HALO2, -1) ;

		// update all halos
		gridLoc.exchangeHalos(MPI_COMM_MODE_SENDRECV) ;
		gridLoc.write(caseName+"Loc") ;

		// fill gridRef with expected values
		gridRef.fill(ALL_POINTS, FUNC_LINEAR, FUNC_LINEAR, FUNC_LINEAR, a1, a2, a3, a4) ;
		gridRef.write(caseName+"Ref") ;

		checkAllProcGridL1Err(ALL_POINTS, gridLoc, gridRef, MAX_ERR_FLOAT) ;
	}

	{
		//----------------------------------
		// check FD coefficients
		// Sum of coef should be equal to 0
		//----------------------------------
		print_blank() ;
		string caseName = testCaseName + "FDcoef" ;
		printInfo(MASTER, " * Case", caseName) ;

		auto fdOrder = Config::Instance()->fdOrder ;
		checkFloatDiff(0.0, getSumFD_D2Coef(fdOrder), MAX_ERR_FLOAT) ;
	}

	{
		//-------------------------------------
		// check write and read grid from disk
		// write GRID_GLOBAL ALL_POINTS
		// read GRID_GLOBAL ALL_POINTS
		//-------------------------------------
		if (Config::Instance()->writeGrid)
		{
			print_blank();
			string caseName = testCaseName + "WriteReadAllPointsGlobalGrid";
			printInfo(MASTER, " * Case", caseName);

			// create global grid and write on disk
			auto gridGlobWrite2 = Grid_Factory::create(gridMode, GRID_GLOBAL);
			Grid &gridGlobWrite = *gridGlobWrite2;
			gridGlobWrite.initializeGrid();
			gridGlobWrite.defineUnitGrid();
			Myfloat amp = 999.0;
			gridGlobWrite.fill(ALL_POINTS, FUNC_LINEAR, FUNC_LINEAR, FUNC_LINEAR, 0.0, 0.0, 0.0, amp);
			gridGlobWrite.write(ALL_POINTS, caseName + "Out");

			// read grid from disk
			auto gridGlobRead2 = Grid_Factory::create(gridMode, GRID_GLOBAL);
			Grid &gridGlobRead = *gridGlobRead2;
			gridGlobRead.initializeGrid();
			gridGlobRead.defineUnitGrid();
			MPI_Barrier(MPI_COMM_WORLD) ;
			gridGlobRead.read(ALL_POINTS, "UtilWriteReadAllPointsGlobalGridOut.proc0.grid.bin");

			// check grids are identical
			//checkAllProcGridL1Err(ALL_POINTS, gridGlobWrite, gridGlobRead, MAX_ERR_FLOAT);
			checkGridL1Err(ALL_POINTS, gridGlobWrite, gridGlobRead, MAX_ERR_FLOAT);
		}
	}

	{
		//-------------------------------------
		// check write and read grid from disk
		// write GRID_GLOBAL INNER_POINTS
		// read GRID_GLOBAL INNER_POINTS
		//-------------------------------------
		if (Config::Instance()->writeGrid)
		{
			print_blank();
			string caseName = testCaseName + "WriteReadInnerPointsGlobalGrid";
			printInfo(MASTER, " * Case", caseName);

			// create global grid and write on disk
			auto gridGlobWrite2 = Grid_Factory::create(gridMode, GRID_GLOBAL);
			Grid &gridGlobWrite = *gridGlobWrite2;
			gridGlobWrite.initializeGrid();
			gridGlobWrite.defineUnitGrid();
			Myfloat amp = 999.0;
			gridGlobWrite.fill(INNER_POINTS, FUNC_LINEAR, FUNC_LINEAR, FUNC_LINEAR, 0.0, 0.0, 0.0, amp);
			gridGlobWrite.write(INNER_POINTS, caseName + "Out");

			// read grid from disk
			auto gridGlobRead2 = Grid_Factory::create(gridMode, GRID_GLOBAL);
			Grid &gridGlobRead = *gridGlobRead2;
			gridGlobRead.initializeGrid();
			gridGlobRead.defineUnitGrid();
			MPI_Barrier(MPI_COMM_WORLD) ;
			gridGlobRead.read(INNER_POINTS, "UtilWriteReadInnerPointsGlobalGridOut.proc0.grid.bin");

			// check grids are identical
			//checkAllProcGridL1Err(INNER_POINTS, gridGlobWrite, gridGlobRead, MAX_ERR_FLOAT);
			checkGridL1Err(INNER_POINTS, gridGlobWrite, gridGlobRead, MAX_ERR_FLOAT);
		}
	}

	{
		//-------------------------------------
		// check read grid from disk		
		// read GRID_LOCAL INNER_POINTS
		//-------------------------------------
		if (Config::Instance()->writeGrid)
		{
			print_blank();
			string caseName = testCaseName + "ReadInnerPointsLocalGrid";
			printInfo(MASTER, " * Case", caseName);

			// create local grid (reference grid)
			auto gridLocRef2 = Grid_Factory::create(gridMode, GRID_LOCAL);
			Grid &gridLocRef = *gridLocRef2;
			gridLocRef.initializeGrid();
			gridLocRef.defineUnitGrid();
			Myfloat amp = 999.0;
			gridLocRef.fill(INNER_POINTS, FUNC_LINEAR, FUNC_LINEAR, FUNC_LINEAR, 0.0, 0.0, 0.0, amp);
			gridLocRef.write(INNER_POINTS, caseName + "Out");

			// read grid from disk (output from previous test case)
			auto gridLocRead2 = Grid_Factory::create(gridMode, GRID_LOCAL);
			Grid &gridLocRead = *gridLocRead2;
			gridLocRead.initializeGrid();
			gridLocRead.defineUnitGrid();
			MPI_Barrier(MPI_COMM_WORLD) ;
			gridLocRead.read(INNER_POINTS, "UtilWriteReadInnerPointsGlobalGridOut.proc0.grid.bin");

			// check grids are identical
			checkAllProcGridL1Err(INNER_POINTS, gridLocRef, gridLocRead, MAX_ERR_FLOAT);
		}
	}

	{
		//-------------------------------------
		// check write grid from disk		
		// write GRID_LOCAL INNER_POINTS
		//-------------------------------------
		if (Config::Instance()->writeGrid)
		{
			print_blank();
			string caseName = testCaseName + "WriteInnerPointsLocalGrid";
			printInfo(MASTER, " * Case", caseName);

			// create local grid
			auto gridLocWrite2 = Grid_Factory::create(gridMode, GRID_LOCAL);
			Grid &gridLocWrite = *gridLocWrite2;
			gridLocWrite.initializeGrid();
			gridLocWrite.defineUnitGrid();
			Myfloat amp = 999.0;
			gridLocWrite.fill(INNER_POINTS, FUNC_LINEAR, FUNC_LINEAR, FUNC_LINEAR, 0.0, 0.0, 0.0, amp);
			
			// write global grid
			gridLocWrite.writeGlobal(INNER_POINTS, caseName + "Out", 0);

			// read grid from disk
			auto gridLocRead2 = Grid_Factory::create(gridMode, GRID_LOCAL);
			Grid &gridLocRead = *gridLocRead2;
			gridLocRead.initializeGrid();
			gridLocRead.defineUnitGrid();
			MPI_Barrier(MPI_COMM_WORLD) ;
			gridLocRead.read(INNER_POINTS, "UtilWriteInnerPointsLocalGridOut.global.grid.bin");

			// check grids are identical
			checkAllProcGridL1Err(INNER_POINTS, gridLocWrite, gridLocRead, MAX_ERR_FLOAT);
		}
	}

	//-------------------
	// check Acquisition
	//-------------------
	if (Config::Instance()->writeGrid)
	{		
		print_blank();
		string caseName = testCaseName + "Acquisition";
		printInfo(MASTER, " * Case", caseName);

		// create grid
		auto grid2 = Grid_Factory::create(gridMode, GRID_LOCAL);
		Grid &grid = *grid2;
		grid.initializeGrid();
		grid.defineUnitGrid();
		Myfloat amp = 999.0;
		grid.fill(INNER_POINTS, amp);

		// create acquisition
		auto acqui2 = DataAcquisition_Factory::create(gridMode) ;		
		DataAcquisition &acqui = *acqui2 ;
		Myint nt = Config::Instance()->nt ;
		acqui.initialize(ACQUI_BUILTIN, grid, nt) ;

		// loop over nt
		for (Myint it=0; it<nt; it++)
		{
			// record trace
			acqui.recordTrace(grid, FORWARD, it, 1.0) ;
		}

		// write traces
		acqui.writeTrace(caseName);
		MPI_Barrier(MPI_COMM_WORLD);

		//------------------
		// check file size
		//------------------

		if (myMpiRank == 0)
		{

			// get size of trace file
			ifstream in_file;
			string file_in = caseName + ".trace.bin";
			in_file.open(file_in.c_str(), ios::binary);
			assert(in_file.is_open());
			in_file.seekg(0, ios_base::end);
			Myint64 inFileSize = in_file.tellg();
			

			// read header
			ifstream in_file2(caseName + ".trace.info");
			assert(in_file2.is_open());
			Myint n1, n2, n3;
			in_file2 >> n1;
			in_file2 >> n2;
			in_file2 >> n3;
			in_file2.close();
			printDebug(MID_DEBUG, " Trace Header n1", n1);
			printDebug(MID_DEBUG, " Trace Header n2", n2);
			printDebug(MID_DEBUG, " Trace Header n3", n3);
			Myint64 expectedFileSize = n1 * n2 * sizeof(Myfloat);

			printInfo(MASTER, " File size on disk", inFileSize);
			printInfo(MASTER, " File size from header", expectedFileSize);
			checkIntegerDiff(inFileSize, expectedFileSize);

			//---------------------------
			// check content of the file
			//---------------------------
			in_file.seekg(0, ios_base::beg);
			Myfloat* valInFile = new Myfloat[n1*n2*n3] ;
			in_file.read((char *)valInFile, sizeof(Myfloat) * n1*n2*n3) ;
			in_file.close();

			// get max diff between values in file and expected values

			Myfloat maxAbsErr = 0 ;
			for (Myint64 i3 = 0; i3 < n3; i3++)
			{
				for (Myint64 i2 = 0; i2 < n2; i2++)
				{
					for (Myint64 i1 = 0; i1 < n1; i1++)
					{
						Myint64 ii = i1 + i2*n1 + i3*(n1*n2) ;
						if (fabs(valInFile[ii] - amp) > maxAbsErr) maxAbsErr = fabs(valInFile[ii] - amp) ;
					}
				}
			}
			printDebug(MID_DEBUG, "maxErr", maxAbsErr) ;
			checkFloatDiff(maxAbsErr, 0.0, MAX_ERR_FLOAT) ;
		}
	}

	this->finalize() ;

	printDebug(MID_DEBUG, "Out TestCase_Util::run") ;
	return(RTN_CODE_OK) ;
}

} // namespace hpcscan

