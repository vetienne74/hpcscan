
//-------------------------------------------------------------------------------------------------------
// test cases for internal use by developers
// these cases are used to check the good behavior of internal funcions
// and not for the purpose of performance measurements
//-------------------------------------------------------------------------------------------------------

#include "testCase_Util.h"

#include <algorithm>
#include <cfloat>

#include <mpi.h>
#include <omp.h>

#include "config.h"
#include "constant.h"
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
		if (nproc_world == 1)
		{
			checkBoolDiff(sameSize2, true) ;
		}
		else
		{
			checkBoolDiff(sameSize2, false) ;
		}
	}
	{
		//-----------------------------
		// check number of Inner points
		//-----------------------------
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

		if (nproc_world == 1)
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
		//---------------------
		// check UpdateI1HALO1
		//---------------------
		print_blank() ;
		string caseName = testCaseName + "UpdateI1HALO1" ;
		printInfo(MASTER, " * Case", caseName) ;

		// fill gridLoc with proc rank
		gridLoc.fill(ALL_POINTS, myid_world) ;
		gridLoc.exchangeHalo(MPI_COMM_MODE_SENDRECV, I1HALO1) ;
		gridLoc.write(caseName+"Loc") ;

		// fill gridRef with expected values
		gridRef.fill(ALL_POINTS, myid_world) ;
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
		gridLoc.fill(ALL_POINTS, myid_world) ;
		gridLoc.exchangeHalo(MPI_COMM_MODE_SENDRECV, I1HALO2) ;
		gridLoc.write(caseName+"Loc") ;

		// fill gridRef with expected values
		gridRef.fill(ALL_POINTS, myid_world) ;
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
		gridLoc.fill(ALL_POINTS, myid_world) ;
		gridLoc.exchangeHalo(MPI_COMM_MODE_SENDRECV, I2HALO1) ;
		gridLoc.write(caseName+"Loc") ;

		// fill gridRef with expected values
		gridRef.fill(ALL_POINTS, myid_world) ;
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
		gridLoc.fill(ALL_POINTS, myid_world) ;
		gridLoc.exchangeHalo(MPI_COMM_MODE_SENDRECV, I2HALO2) ;
		gridLoc.write(caseName+"Loc") ;

		// fill gridRef with expected values
		gridRef.fill(ALL_POINTS, myid_world) ;
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
		gridLoc.fill(ALL_POINTS, myid_world) ;
		gridLoc.exchangeHalo(MPI_COMM_MODE_SENDRECV, I3HALO1) ;
		gridLoc.write(caseName+"Loc") ;

		// fill gridRef with expected values
		gridRef.fill(ALL_POINTS, myid_world) ;
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
		gridLoc.fill(ALL_POINTS, myid_world) ;
		gridLoc.exchangeHalo(MPI_COMM_MODE_SENDRECV, I3HALO2) ;
		gridLoc.write(caseName+"Loc") ;

		// fill gridRef with expected values
		gridRef.fill(ALL_POINTS, myid_world) ;
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
		//-----------------------
		// check grid coordinates
		//-----------------------
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
		//---------------------
		// check exchangeHalos
		//---------------------
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

	this->finalize() ;

	printDebug(MID_DEBUG, "Out TestCase_Util::run") ;
	return(RTN_CODE_OK) ;
}

} // namespace hpcscan

