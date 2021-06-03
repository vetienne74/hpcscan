
//-------------------------------------------------------------------------------------------------------
// Parent class of all test cases
// provides basic functions
//-------------------------------------------------------------------------------------------------------

#include "testCase.h"

#include <cmath> // for fabs

#include <mpi.h>
#include <omp.h>

#include "config.h"
#include "constant.h"
#include "global.h"
#include "grid_Factory.h"
#include "output_report.h"

namespace hpcscan {

//-------------------------------------------------------------------------------------------------------

Rtn_code TestCase::initialize(void)
{
	printDebug(MID_DEBUG, "In TestCase::initialize") ;

	// run this testCase
	if ( (Config::Instance()->testCaseName.compare(this->testCaseName) == 0) ||
			(Config::Instance()->testCaseName.compare("All") == 0) )
	{
		print_blank() ;
		print_line2() ;
		printInfo(MASTER, " TestCase name\t", testCaseName) ;
		printInfo(MASTER, " TestCase version", testCaseVersion) ;
		printInfo(MASTER, " TestCase mode\t", Config::Instance()->testMode) ;

		// try to create and initialize a grid to see if everything is all right
		auto gridTest = Grid_Factory::create(Config::Instance()->testMode, GRID_LOCAL) ;
		if (gridTest == nullptr)
		{
			printError("In TestCase::initialize, Not supported or invalid testMode") ;
			return(RTN_CODE_KO) ;
		}
		if (gridTest->initializeGrid() != RTN_CODE_OK)
		{
			printError("In TestCase::initialize, Can not initialize grid with testMode") ;
			return(RTN_CODE_KO) ;
		}

		// open perf log file
		if (myMpiRank == 0)
		{
			string file_name = "hpcscan.perf." + testCaseName + ".log";
			perfLogFile.open(file_name, ios::app) ;

			// all strings first
			perfLogFile << Config::Instance()->hostName << " " ;
			//perfLogFile << Config::Instance()->userName << " " ;
			perfLogFile << Config::Instance()->testCaseName << " " ;
			perfLogFile << Config::Instance()->testMode << " " ;
			perfLogFile << Config::Instance()->propagator << " " ;

			// numeric values follow
			perfLogFile << nMpiProc << " " ; // 1
			perfLogFile << Config::Instance()->nsub1 << " " ; // 2
			perfLogFile << Config::Instance()->nsub2 << " " ; // 3
			perfLogFile << Config::Instance()->nsub3 << " " ; // 4
			perfLogFile << omp_get_max_threads() << " " ; // 5
			perfLogFile << Config::Instance()->n1 << " " ; // 6
			perfLogFile << Config::Instance()->n2 << " " ; // 7
			perfLogFile << Config::Instance()->n3 << " " ; // 8
			perfLogFile << Config::Instance()->fdOrder << " " ; // 9

		}

		// start timer
		testCaseStart = MPI_Wtime() ;

		printDebug(MID_DEBUG, "Out TestCase::initialize") ;
		return(RTN_CODE_OK) ;
	}
	// do not run this testCase
	else
	{
		printDebug(MID_DEBUG, " * SKIP THIS TESTCASE * ") ;
		printDebug(MID_DEBUG, "Out TestCase::initialize") ;
		return(RTN_CODE_KO) ;
	}
}

//-------------------------------------------------------------------------------------------------------

void TestCase::finalize(void)
{
	printDebug(MID_DEBUG, "In TestCase::finalize") ;

	// end timer
	testCaseEnd = MPI_Wtime() ;

	// close perf log file
	if (myMpiRank == 0)
	{
		perfLogFile.close() ;
	}

	print_blank() ;
	printInfo(MASTER, " Case run time (s)", testCaseEnd-testCaseStart) ;
	print_line2() ;

	printDebug(MID_DEBUG, "Out TestCase::finalize") ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code TestCase::checkGridL1Err(Point_type pointType, const Grid& gridIn, const Grid& gridRef, Myfloat maxAllowed)
{

	printDebug(FULL_DEBUG, "In TestCase::checkGridL1Err") ;

	// check L1 error
	Myfloat L1Err = gridIn.L1Err(pointType, gridRef) ;
	if ((!std::isnan(L1Err)) && (L1Err >= 0.0) && (L1Err < maxAllowed))
	{
		printInfo(MASTER, " Test PASSED / L1 err.", L1Err) ;
	}
	else
	{
		printInfo(MASTER, ">> Test FAILED! / L1 err.", L1Err) ;
	}

	printDebug(FULL_DEBUG, "Out TestCase::checkGridL1Err") ;
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code TestCase::checkAllProcGridL1Err(Point_type pointType, const Grid& gridIn, const Grid& gridRef, Myfloat maxAllowed)
{

	printDebug(FULL_DEBUG, "In TestCase::checkAllProcGridL1Err") ;

	// check L1 error
	Myfloat L1ErrLoc = gridIn.L1Err(pointType, gridRef) ;
	Myfloat L1Err ;

	MPI_Reduce(&L1ErrLoc, &L1Err, 1, MPI_MYFLOAT, MPI_MAX, 0, MPI_COMM_WORLD);

	if ((!std::isnan(L1Err)) && (L1Err >= 0.0) && (L1Err < maxAllowed))
	{
		printInfo(MASTER, " Test PASSED / L1 err.", L1Err) ;
	}
	else
	{
		printInfo(MASTER, ">> Test FAILED! / L1 err.", L1Err) ;
	}

	printDebug(FULL_DEBUG, "Out TestCase::checkAllProcGridL1Err") ;
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code TestCase::checkGridMaxErr(Point_type pointType, const Grid& gridIn, const Grid& gridRef, Myfloat maxAllowed)
{

	printDebug(FULL_DEBUG, "In TestCase::checkGridMaxErr") ;

	// check Max error
	Myfloat maxErrLoc = gridIn.maxErr(pointType, gridRef) ;
	Myfloat maxErr ;
	MPI_Reduce(&maxErrLoc, &maxErr, 1, MPI_MYFLOAT, MPI_MAX, 0, MPI_COMM_WORLD);

	if ((!std::isnan(maxErr)) && (maxErr >= 0.0) && (maxErr < maxAllowed))
	{
		printInfo(MASTER, " Test PASSED / Max err.", maxErr) ;
	}
	else
	{
		printInfo(MASTER, ">> Test FAILED! / Max err.", maxErr) ;
	}

	printDebug(FULL_DEBUG, "Out TestCase::checkGridMaxErr") ;
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code TestCase::checkAllProcGridMaxErr(Point_type pointType, const Grid& gridIn, const Grid& gridRef, Myfloat maxAllowed)
{

	printDebug(FULL_DEBUG, "In TestCase::checkAllProcGridMaxErr") ;

	// check Max error
	Myfloat maxErr = gridIn.maxErr(pointType, gridRef) ;
	if ((!std::isnan(maxErr)) && (maxErr >= 0.0) && (maxErr < maxAllowed))
	{
		printInfo(MASTER, " Test PASSED / Max err.", maxErr) ;
	}
	else
	{
		printInfo(MASTER, ">> Test FAILED! / Max err.", maxErr) ;
	}

	printDebug(FULL_DEBUG, "Out TestCase::checkAllProcGridMaxErr") ;
	return(RTN_CODE_OK) ;
}


//-------------------------------------------------------------------------------------------------------

Rtn_code TestCase::checkIntegerDiff(Myint val1, Myint val2)
{

	printDebug(FULL_DEBUG, "In TestCase::checkIntegerDiff") ;

	if (val1 == val2)
	{
		printInfo(MASTER, " Test PASSED") ;
	}
	else
	{
		printInfo(MASTER, ">> Test FAILED!") ;
		printDebug(LIGHT_DEBUG, "val1", val1) ;
		printDebug(LIGHT_DEBUG, "val2", val2) ;
		return(RTN_CODE_KO) ;
	}

	printDebug(FULL_DEBUG, "Out TestCase::checkIntegerDiff") ;
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code TestCase::checkBoolDiff(bool val1, bool val2)
{

	printDebug(FULL_DEBUG, "In TestCase::checkBoolDiff") ;

	if (val1 == val2)
	{
		printInfo(MASTER, " Test PASSED") ;
	}
	else
	{
		printInfo(MASTER, ">> Test FAILED!") ;
		printDebug(LIGHT_DEBUG, "val1", val1) ;
		printDebug(LIGHT_DEBUG, "val2", val2) ;
		return(RTN_CODE_KO) ;
	}

	printDebug(FULL_DEBUG, "Out TestCase::checkBoolDiff") ;
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code TestCase::checkFloatDiff(Myfloat val1, Myfloat val2, Myfloat maxAllowed)
{

	printDebug(FULL_DEBUG, "In TestCase::checkFloatDiff") ;

	Myfloat relErr2 = relErr(val1, val2) ;

	if (relErr2 <= maxAllowed)
	{
		printInfo(MASTER, " Test PASSED / Rel err.", relErr2) ;
	}
	else
	{
		printInfo(MASTER, ">> Test FAILED! / Rel err.", relErr2) ;
		return(RTN_CODE_KO) ;
	}

	printDebug(FULL_DEBUG, "Out TestCase::checkFloatDiff") ;
	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Myfloat TestCase::relErr(Myfloat val1, Myfloat val2)
{
	printDebug(LIGHT_DEBUG, "In TestCase::relErr") ;

	Myfloat retVal ;
	printDebug(LIGHT_DEBUG, "val1", val1) ;
	printDebug(LIGHT_DEBUG, "val2", val2) ;
	if (fabs(val2) > MAX_ERR_FLOAT)
	{
		retVal = fabs(val1 - val2) / fabs(val2) ;
	}
	else
	{
		retVal = fabs(val1 - val2) ;
	}
	printDebug(LIGHT_DEBUG, "retVal", retVal) ;
	printDebug(LIGHT_DEBUG, "Out TestCase::relErr") ;
	return(retVal) ;
}


} // namespace hpcscan

