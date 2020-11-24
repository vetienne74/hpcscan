
//-------------------------------------------------------------------------------------------------------
// test cases to measure the memory bandwidth
// baseline kernel - no optimization
//-------------------------------------------------------------------------------------------------------

#include "testCase_Memory.h"

#include <algorithm> // for min()
#include <cfloat>    // for FLT_MAX

#include <mpi.h>
#include <omp.h>

#include "config.h"
#include "constant.h"
#include "global.h"
#include "grid_Factory.h"
#include "output_report.h"

namespace hpcscan {

TestCase_Memory::TestCase_Memory(void)
{
	testCaseName    = "Memory" ;
	testCaseVersion = "Standard implementation" ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code TestCase_Memory::run(void)
{

	printDebug(MID_DEBUG, "In TestCase_Memory::run") ;

	if (this->initialize() == RTN_CODE_KO) return (RTN_CODE_KO) ;

	const Myfloat a1 = Config::Instance()->param1 ;
	const Myfloat a2 = Config::Instance()->param2 ;

	const string gridMode = Config::Instance()->testMode ;
	auto Ugrid2 = Grid_Factory::create(gridMode, GRID_GLOBAL) ;
	auto Vgrid2 = Grid_Factory::create(gridMode, GRID_GLOBAL) ;
	auto Wgrid2 = Grid_Factory::create(gridMode, GRID_GLOBAL) ;
	auto Rgrid2 = Grid_Factory::create(gridMode, GRID_GLOBAL) ;
	Grid &Ugrid = *Ugrid2 ;
	Grid &Vgrid = *Vgrid2 ;
	Grid &Wgrid = *Wgrid2 ;
	Grid &Rgrid = *Rgrid2 ; // reference grid
	Ugrid.initializeGrid() ;
	Vgrid.initializeGrid() ;
	Wgrid.initializeGrid() ;
	Rgrid.initializeGrid() ;

	if (myid_world == 0) Ugrid.info() ;

	const Myint64 nGridPoint = Ugrid.getNumberOfGridPoint(GRID_GLOBAL, ALL_POINTS) ;

	Myfloat * const u = Ugrid.grid_3d ;
	Myfloat * const v = Vgrid.grid_3d ;
	Myfloat * const w = Wgrid.grid_3d ;

	// one block per case
	{
		//============================================
		// fill grid
		// W = coef
		//============================================

		print_blank() ;
		string caseName = testCaseName + "FillGrid" ;
		printInfo(MASTER, " * Case", caseName) ;

		Wgrid.fill(ALL_POINTS, 0.0) ;

		double testCase_time_best = FLT_MAX ;

		Myint ntry = Config::Instance()->ntry ;

		for (Myint itry = 0; itry < ntry; itry++)
		{
			double t0 = MPI_Wtime() ;
#pragma omp parallel for
			for (Myint64 ii=0; ii<nGridPoint; ii++)
			{
				w[ii] = a1 ;
			}
			double t1 = MPI_Wtime() ;

			double testCase_time = t1-t0 ;
			Myfloat testCase_bw = nGridPoint*sizeof(Myfloat)/testCase_time/1e9 ;
			printDebug(LIGHT_DEBUG, "Time", testCase_time) ;
			printDebug(LIGHT_DEBUG, "Speed", testCase_bw) ;
			testCase_time_best = min(testCase_time_best, testCase_time) ;

			// check testCase results
			if (itry == 0)
			{
				Myfloat minVal = Wgrid.getMin(ALL_POINTS) ;
				checkFloatDiff(minVal, a1, MAX_ERR_FLOAT) ;
				Myfloat maxVal = Wgrid.getMax(ALL_POINTS) ;
				checkFloatDiff(maxVal, a1, MAX_ERR_FLOAT) ;
			}

		}

		printInfo(MASTER, " Best achieved GByte/s",
				Myfloat(nGridPoint*sizeof(Myfloat)/testCase_time_best/1.e9)) ;
		printInfo(MASTER, " Best achieved GPoint/s",
				Myfloat(nGridPoint/testCase_time_best/1.e9)) ;
	}

	{
		//============================================
		// copy grid
		// W = U
		//============================================

		print_blank() ;
		string caseName = testCaseName + "CopyGrid" ;
		printInfo(MASTER, " * Case", caseName) ;

		Ugrid.fill(ALL_POINTS, a1) ;
		Wgrid.fill(ALL_POINTS, 0.0) ;

		double testCase_time_best = FLT_MAX ;

		Myint ntry = Config::Instance()->ntry ;

		for (Myint itry = 0; itry < ntry; itry++)
		{
			double t0 = MPI_Wtime() ;
#pragma omp parallel for
			for (Myint64 ii=0; ii<nGridPoint; ii++)
			{
				w[ii] = u[ii] ;
			}
			double t1 = MPI_Wtime() ;

			double testCase_time = t1-t0 ;
			Myfloat testCase_bw = nGridPoint*2*sizeof(Myfloat)/testCase_time/1e9 ;
			printDebug(LIGHT_DEBUG, "Time", testCase_time) ;
			printDebug(LIGHT_DEBUG, "Speed", testCase_bw) ;
			testCase_time_best = min(testCase_time_best, testCase_time) ;

			// check testCase results
			if (itry == 0)
			{
				checkGridL1Err(ALL_POINTS, Wgrid, Ugrid, MAX_ERR_FLOAT) ;
				checkGridMaxErr(ALL_POINTS, Wgrid, Ugrid, MAX_ERR_FLOAT) ;
			}

		}

		printInfo(MASTER, " Best achieved GByte/s",
				Myfloat(nGridPoint*2*sizeof(Myfloat)/testCase_time_best/1.e9)) ;
		printInfo(MASTER, " Best achieved GPoint/s",
				Myfloat(nGridPoint/testCase_time_best/1.e9)) ;
	}

	{
		//============================================
		// add grid
		// W = U + V
		//============================================

		print_blank() ;
		string caseName = testCaseName + "AddGrid" ;
		printInfo(MASTER, " * Case", caseName) ;

		Ugrid.fill(ALL_POINTS, a1) ;
		Vgrid.fill(ALL_POINTS, a2) ;
		Wgrid.fill(ALL_POINTS, 0.0) ;
		Rgrid.fill(ALL_POINTS, a1+a2) ;

		double testCase_time_best = FLT_MAX ;

		Myint ntry = Config::Instance()->ntry ;

		for (Myint itry = 0; itry < ntry; itry++)
		{
			double t0 = MPI_Wtime() ;
#pragma omp parallel for
			for (Myint64 ii=0; ii<nGridPoint; ii++)
			{
				w[ii] = u[ii] + v[ii];
			}
			double t1 = MPI_Wtime() ;

			double testCase_time = t1-t0 ;
			Myfloat testCase_bw = nGridPoint*3*sizeof(Myfloat)/testCase_time/1e9 ;
			printDebug(LIGHT_DEBUG, "Time", testCase_time) ;
			printDebug(LIGHT_DEBUG, "Speed", testCase_bw) ;
			testCase_time_best = min(testCase_time_best, testCase_time) ;

			// check testCase results
			if (itry == 0)
			{
				checkGridL1Err(ALL_POINTS, Wgrid, Rgrid, MAX_ERR_FLOAT) ;
				checkGridMaxErr(ALL_POINTS, Wgrid, Rgrid, MAX_ERR_FLOAT) ;
			}

		}

		printInfo(MASTER, " Best achieved GByte/s",
				Myfloat(nGridPoint*3*sizeof(Myfloat)/testCase_time_best/1.e9)) ;
		printInfo(MASTER, " Best achieved GPoint/s",
				Myfloat(nGridPoint/testCase_time_best/1.e9)) ;
	}

	{
		//============================================
		// multiply grid
		// W = U * V
		//============================================

		print_blank() ;
		string caseName = testCaseName + "MultiplyGrid" ;
		printInfo(MASTER, " * Case", caseName) ;

		Ugrid.fill(ALL_POINTS, a1) ;
		Vgrid.fill(ALL_POINTS, a2) ;
		Wgrid.fill(ALL_POINTS, 0.0) ;
		Rgrid.fill(ALL_POINTS, a1*a2) ;

		double testCase_time_best = FLT_MAX ;

		Myint ntry = Config::Instance()->ntry ;

		for (Myint itry = 0; itry < ntry; itry++)
		{
			double t0 = MPI_Wtime() ;
#pragma omp parallel for
			for (Myint64 ii=0; ii<nGridPoint; ii++)
			{
				w[ii] = u[ii] * v[ii];
			}
			double t1 = MPI_Wtime() ;

			double testCase_time = t1-t0 ;
			Myfloat testCase_bw = nGridPoint*3*sizeof(Myfloat)/testCase_time/1e9 ;
			printDebug(LIGHT_DEBUG, "Time", testCase_time) ;
			printDebug(LIGHT_DEBUG, "Speed", testCase_bw) ;
			testCase_time_best = min(testCase_time_best, testCase_time) ;

			// check testCase results
			if (itry == 0)
			{
				checkGridL1Err(ALL_POINTS, Wgrid, Rgrid, MAX_ERR_FLOAT) ;
				checkGridMaxErr(ALL_POINTS, Wgrid, Rgrid, MAX_ERR_FLOAT) ;
			}

		}

		printInfo(MASTER, " Best achieved GByte/s",
				Myfloat(nGridPoint*3*sizeof(Myfloat)/testCase_time_best/1.e9)) ;
		printInfo(MASTER, " Best achieved GPoint/s",
				Myfloat(nGridPoint/testCase_time_best/1.e9)) ;
	}

	{
		//============================================
		// add and update grid
		// W = W + U
		//============================================

		print_blank() ;
		string caseName = testCaseName + "AddUpdateGrid" ;
		printInfo(MASTER, " * Case", caseName) ;

		Ugrid.fill(ALL_POINTS, a1) ;
		Wgrid.fill(ALL_POINTS, a2) ;
		Rgrid.fill(ALL_POINTS, a1+a2) ;

		double testCase_time_best = FLT_MAX ;

		Myint ntry = Config::Instance()->ntry ;

		for (Myint itry = 0; itry < ntry; itry++)
		{
			double t0 = MPI_Wtime() ;
#pragma omp parallel for
			for (Myint64 ii=0; ii<nGridPoint; ii++)
			{
				w[ii] = w[ii] + u[ii];
			}
			double t1 = MPI_Wtime() ;

			double testCase_time = t1-t0 ;
			Myfloat testCase_bw = nGridPoint*3*sizeof(Myfloat)/testCase_time/1e9 ;
			printDebug(LIGHT_DEBUG, "Time", testCase_time) ;
			printDebug(LIGHT_DEBUG, "Speed", testCase_bw) ;
			testCase_time_best = min(testCase_time_best, testCase_time) ;

			// check testCase results
			if (itry == 0)
			{
				checkGridL1Err(ALL_POINTS, Wgrid, Rgrid, MAX_ERR_FLOAT) ;
				checkGridMaxErr(ALL_POINTS, Wgrid, Rgrid, MAX_ERR_FLOAT) ;
			}

		}

		printInfo(MASTER, " Best achieved GByte/s",
				Myfloat(nGridPoint*3*sizeof(Myfloat)/testCase_time_best/1.e9)) ;
		printInfo(MASTER, " Best achieved GPoint/s",
				Myfloat(nGridPoint/testCase_time_best/1.e9)) ;
	}

	this->finalize() ;

	printDebug(MID_DEBUG, "Out TestCase_Memory::run") ;
	return(RTN_CODE_OK) ;
}

} // namespace hpcscan

