
//-------------------------------------------------------------------------------------------------------
// test cases to benchmark basic operations on grids
// baseline kernel - no optimization
//-------------------------------------------------------------------------------------------------------

#include "testCase_Grid.h"

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

TestCase_Grid::TestCase_Grid(void)
								{
	testCaseName    = "Grid" ;
	testCaseVersion = "Standard implementation" ;
								}

//-------------------------------------------------------------------------------------------------------

Rtn_code TestCase_Grid::run(void)
{

	printDebug(MID_DEBUG, "In TestCase_Grid::run") ;

	if (this->initialize() == RTN_CODE_KO) return (RTN_CODE_KO) ;

	// create a block for each test

	// create grids for the following test cases
	const Myfloat64 a1 = Config::Instance()->param1 ;
	const Myfloat64 a2 = Config::Instance()->param2 ;
	const Myfloat64 a3 = Config::Instance()->param3 ;
	const Myfloat64 a4 = Config::Instance()->param4 ;

	const string gridMode = Config::Instance()->testMode ;
	auto W2grid = Grid_Factory::create(gridMode, GRID_LOCAL) ;
	auto R2grid = Grid_Factory::create(gridMode, GRID_LOCAL) ;
	Grid &Wgrid = *W2grid ;
	Grid &Rgrid = *R2grid ;
	Wgrid.initializeGrid() ;
	Rgrid.initializeGrid() ;
	Wgrid.info() ;

	Myfloat * const w = Wgrid.grid_3d ;

	Myint64 nGridPointLoc  = Wgrid.getNumberOfGridPoint(GRID_LOCAL, INNER_POINTS) ;
	Myint64 nGridPointGlob = Wgrid.getNumberOfGridPoint(GRID_GLOBAL, INNER_POINTS) ;
	printInfo(MASTER, " Inner pts (Loc.)", nGridPointLoc) ;
	printInfo(MASTER, " Inner pts (Glob.)", nGridPointGlob) ;

	Myint ntry = Config::Instance()->ntry ;

	// for perf log
	Myfloat FillGB=0, FillGpoint=0, MaxErrGB=0, MaxErrGpoint=0, L1ErrGB=0, L1ErrGpoint=0 ;
	Myfloat GetSumAbsGB=0, GetSumAbsGpoint=0, GetSumAbsDiffGB=0, GetSumAbsDiffGpoint=0 ;
	Myfloat GetMaxGB=0, GetMaxGpoint=0, GetMinGB=0, GetMinGpoint=0 ;
	Myfloat UpdatePressureGB=0, UpdatePressureGpoint=0, ApplyBoundaryConditionGB=0, ApplyBoundaryConditionGpoint=0 ;

	{
		//============================================
		// Fill grid
		// W = coef
		//============================================

		print_blank() ;
		string caseName = testCaseName + "Fill" ;
		printInfo(MASTER, " * Case", caseName) ;

		Rgrid.fill(INNER_POINTS, a1) ;

		double testCase_time_best = FLT_MAX ;

		for (Myint itry = 0; itry < ntry; itry++)
		{
			double t0 = MPI_Wtime() ;
			Wgrid.fill(INNER_POINTS, a1) ;
			// wait all process completed before ending time
			MPI_Barrier(MPI_COMM_WORLD) ;
			double t1 = MPI_Wtime() ;

			double testCase_time = t1-t0 ;
			Myfloat testCase_bw = nGridPointGlob*sizeof(Myfloat)/testCase_time/1e9 ;
			printDebug(LIGHT_DEBUG, "Time", testCase_time) ;
			printDebug(LIGHT_DEBUG, "Speed", testCase_bw) ;
			testCase_time_best = min(testCase_time_best, testCase_time) ;

			// check testCase results
			if (itry == 0)
			{
				checkAllProcGridL1Err(INNER_POINTS, Wgrid, Rgrid, MAX_ERR_FLOAT) ;
				// checkAllProcGridMaxErr(INNER_POINTS, Wgrid, Rgrid, MAX_ERR_FLOAT) ;
			}

		}

		FillGB     = Myfloat(nGridPointGlob*sizeof(Myfloat)/testCase_time_best/1.e9) ;
		FillGpoint = Myfloat(nGridPointGlob/testCase_time_best/1.e9) ;
		printInfo(MASTER, " Best achieved GByte/s", FillGB) ;
		printInfo(MASTER, " Best achieved GPoint/s", FillGpoint) ;
	}

	// initialize grids
	Rgrid.fill(INNER_POINTS, a1) ;
	Wgrid.fill(INNER_POINTS, a2) ;

	// change value of only one point
	if (myid_world == 0) Wgrid.fill(MIDDLE_POINT, 2*a2) ;

	if (false)
	{
		//============================================
		// maxErr = max of abs(W-R)/R
		//============================================

		print_blank() ;
		string caseName = testCaseName + "MaxErr" ;
		printInfo(MASTER, " * Case", caseName) ;

		double testCase_time_best = FLT_MAX ;

		for (Myint itry = 0; itry < ntry; itry++)
		{
			double t0 = MPI_Wtime() ;
			Myfloat maxErrLoc = Wgrid.maxErr(INNER_POINTS, Rgrid) ;
			// wait all process completed before ending time
			MPI_Barrier(MPI_COMM_WORLD) ;
			double t1 = MPI_Wtime() ;

			double testCase_time = t1-t0 ;
			printDebug(LIGHT_DEBUG, "Time", testCase_time) ;
			testCase_time_best = min(testCase_time_best, testCase_time) ;

			// check testCase results
			if (itry == 0)
			{
				Myfloat maxErr ;
				MPI_Reduce(&maxErrLoc, &maxErr, 1, MPI_MYFLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
				Myfloat maxErrRef = abs(2*a2-a1)/a1 ;
				checkFloatDiff(maxErr, maxErrRef, MAX_ERR_FLOAT) ;
			}

		}

		MaxErrGB = Myfloat(2*nGridPointGlob*sizeof(Myfloat)/testCase_time_best/1.e9) ;
		MaxErrGpoint = Myfloat(nGridPointGlob/testCase_time_best/1.e9) ;
		printInfo(MASTER, " Best achieved GByte/s", MaxErrGB) ;
		printInfo(MASTER, " Best achieved GPoint/s", MaxErrGpoint) ;
	}

	{
		//============================================
		// L1Err = sum of abs(W-R)/abs(R)
		//============================================

		print_blank() ;
		string caseName = testCaseName + "L1Err" ;
		printInfo(MASTER, " * Case", caseName) ;

		double testCase_time_best = FLT_MAX ;

		for (Myint itry = 0; itry < ntry; itry++)
		{
			double t0 = MPI_Wtime() ;
			// Myfloat L1Err = Wgrid.allProcL1Err(INNER_POINTS, Rgrid) ;
			Myfloat L1Err = Wgrid.L1Err(INNER_POINTS, Rgrid) ;
			// wait all process completed before ending time
			MPI_Barrier(MPI_COMM_WORLD) ;
			double t1 = MPI_Wtime() ;

			double testCase_time = t1-t0 ;
			printDebug(LIGHT_DEBUG, "Time", testCase_time) ;
			testCase_time_best = min(testCase_time_best, testCase_time) ;

			// check testCase results
			if (itry == 0)
			{
				Myfloat64 sum1 = (nGridPointGlob-1) * abs(a2-a1) + abs(2*a2-a1) ;
				Myfloat64 sum2 = nGridPointGlob * a1 ;
				Myfloat L1ErrRef = sum1 / sum2 ;
				checkFloatDiff(L1Err, L1ErrRef, MAX_ERR_FLOAT) ;
			}

		}

		L1ErrGB = Myfloat(2*nGridPointGlob*sizeof(Myfloat)/testCase_time_best/1.e9) ;
		L1ErrGpoint = Myfloat(nGridPointGlob/testCase_time_best/1.e9) ;
		printInfo(MASTER, " Best achieved GByte/s", L1ErrGB) ;
		printInfo(MASTER, " Best achieved GPoint/s", L1ErrGpoint) ;
	}

	{
		//============================================
		// GetSumAbs = sum of abs(R)
		//============================================

		print_blank() ;
		string caseName = testCaseName + "GetSumAbs" ;
		printInfo(MASTER, " * Case", caseName) ;

		double testCase_time_best = FLT_MAX ;

		for (Myint itry = 0; itry < ntry; itry++)
		{
			double t0 = MPI_Wtime() ;
			Myfloat sumAbs = Rgrid.getSumAbs(INNER_POINTS) ;
			// wait all process completed before ending time
			MPI_Barrier(MPI_COMM_WORLD) ;
			double t1 = MPI_Wtime() ;

			double testCase_time = t1-t0 ;
			printDebug(LIGHT_DEBUG, "Time", testCase_time) ;
			testCase_time_best = min(testCase_time_best, testCase_time) ;

			// check testCase results
			if (itry == 0)
			{
				Myfloat64 sum2 = nGridPointGlob * a1 ;
				checkFloatDiff(sumAbs, sum2, MAX_ERR_FLOAT) ;
			}
		}

		GetSumAbsGB = Myfloat(1*nGridPointGlob*sizeof(Myfloat)/testCase_time_best/1.e9) ;
		GetSumAbsGpoint = Myfloat(nGridPointGlob/testCase_time_best/1.e9) ;
		printInfo(MASTER, " Best achieved GByte/s", GetSumAbsGB) ;
		printInfo(MASTER, " Best achieved GPoint/s", GetSumAbsGpoint) ;
	}

	{
		//============================================
		// GetSumAbs = sum of abs (W-R)
		//============================================

		print_blank() ;
		string caseName = testCaseName + "GetSumAbsDiff" ;
		printInfo(MASTER, " * Case", caseName) ;

		double testCase_time_best = FLT_MAX ;

		for (Myint itry = 0; itry < ntry; itry++)
		{
			double t0 = MPI_Wtime() ;
			Myfloat sumAbsDiff = Wgrid.getSumAbsDiff(INNER_POINTS, Rgrid) ;
			// wait all process completed before ending time
			MPI_Barrier(MPI_COMM_WORLD) ;
			double t1 = MPI_Wtime() ;

			double testCase_time = t1-t0 ;
			printDebug(LIGHT_DEBUG, "Time", testCase_time) ;
			testCase_time_best = min(testCase_time_best, testCase_time) ;

			// check testCase results
			if (itry == 0)
			{
				Myfloat64 sum1 = (nGridPointGlob-1) * abs(a2-a1) + abs(2*a2-a1) ;
				checkFloatDiff(sumAbsDiff, sum1, MAX_ERR_FLOAT) ;
			}
		}

		GetSumAbsDiffGB = Myfloat(2*nGridPointGlob*sizeof(Myfloat)/testCase_time_best/1.e9) ;
		GetSumAbsDiffGpoint = Myfloat(nGridPointGlob/testCase_time_best/1.e9) ;
		printInfo(MASTER, " Best achieved GByte/s", GetSumAbsDiffGB) ;
		printInfo(MASTER, " Best achieved GPoint/s", GetSumAbsDiffGpoint) ;
	}

	{
		//============================================
		// getMax = max of R
		//============================================

		print_blank() ;
		string caseName = testCaseName + "GetMax" ;
		printInfo(MASTER, " * Case", caseName) ;

		double testCase_time_best = FLT_MAX ;

		for (Myint itry = 0; itry < ntry; itry++)
		{
			double t0 = MPI_Wtime() ;
			Myfloat maxVal = Wgrid.getMax(INNER_POINTS) ;
			// wait all process completed before ending time
			MPI_Barrier(MPI_COMM_WORLD) ;
			double t1 = MPI_Wtime() ;

			double testCase_time = t1-t0 ;
			printDebug(LIGHT_DEBUG, "Time", testCase_time) ;
			testCase_time_best = min(testCase_time_best, testCase_time) ;

			// check testCase results
			if (itry == 0)
			{
				Myfloat maxValRef = 2*a2 ;
				checkFloatDiff(maxVal, maxValRef, MAX_ERR_FLOAT) ;
			}

		}

		GetMaxGB = Myfloat(nGridPointGlob*sizeof(Myfloat)/testCase_time_best/1.e9) ;
		GetMaxGpoint = Myfloat(nGridPointGlob/testCase_time_best/1.e9) ;
		printInfo(MASTER, " Best achieved GByte/s", GetMaxGB) ;
		printInfo(MASTER, " Best achieved GPoint/s", GetMaxGpoint) ;
	}

	{
		//============================================
		// getMin = min of R
		//============================================

		print_blank() ;
		string caseName = testCaseName + "GetMin" ;
		printInfo(MASTER, " * Case", caseName) ;

		double testCase_time_best = FLT_MAX ;

		for (Myint itry = 0; itry < ntry; itry++)
		{
			double t0 = MPI_Wtime() ;
			Myfloat minVal = Wgrid.getMin(INNER_POINTS) ;
			// wait all process completed before ending time
			MPI_Barrier(MPI_COMM_WORLD) ;
			double t1 = MPI_Wtime() ;

			double testCase_time = t1-t0 ;
			printDebug(LIGHT_DEBUG, "Time", testCase_time) ;
			testCase_time_best = min(testCase_time_best, testCase_time) ;

			// check testCase results
			if (itry == 0)
			{
				Myfloat minValRef = a2 ;
				checkFloatDiff(minVal, minValRef, MAX_ERR_FLOAT) ;
			}

		}

		GetMinGB = Myfloat(nGridPointGlob*sizeof(Myfloat)/testCase_time_best/1.e9) ;
		GetMinGpoint = Myfloat(nGridPointGlob/testCase_time_best/1.e9) ;
		printInfo(MASTER, " Best achieved GByte/s", GetMinGB) ;
		printInfo(MASTER, " Best achieved GPoint/s", GetMinGpoint) ;
	}

	{
		//============================================
		// updatePressure
		// W = coef
		//============================================

		print_blank() ;
		string caseName = testCaseName + "UpdatePressure" ;
		printInfo(MASTER, " * Case", caseName) ;

		Rgrid.fill(INNER_POINTS, a1) ;

		double testCase_time_best = FLT_MAX ;

		auto prcGrid2   = Grid_Factory::create(gridMode, GRID_LOCAL) ;
		auto coefGrid2  = Grid_Factory::create(gridMode, GRID_LOCAL) ;
		auto laplaGrid2 = Grid_Factory::create(gridMode, GRID_LOCAL) ;
		Grid &prcGrid   = *prcGrid2 ;
		Grid &coefGrid  = *coefGrid2 ;
		Grid &laplaGrid = *laplaGrid2 ;
		prcGrid.initializeGrid() ;
		coefGrid.initializeGrid() ;
		laplaGrid.initializeGrid() ;

		Wgrid.fill(INNER_POINTS, a1) ;
		prcGrid.fill(INNER_POINTS, a2) ;
		coefGrid.fill(INNER_POINTS, a3) ;
		laplaGrid.fill(INNER_POINTS, a4) ;

		Myfloat valRef = 2.0 * a2 - a1 + a3 * a4 ;
		Rgrid.fill(INNER_POINTS, valRef) ;

		for (Myint itry = 0; itry < ntry; itry++)
		{
			double t0 = MPI_Wtime() ;
			Wgrid.updatePressure(INNER_POINTS, prcGrid, coefGrid, laplaGrid) ;
			// wait all process completed before ending time
			MPI_Barrier(MPI_COMM_WORLD) ;
			double t1 = MPI_Wtime() ;

			double testCase_time = t1-t0 ;
			Myfloat testCase_bw = nGridPointGlob*sizeof(Myfloat)/testCase_time/1e9 ;
			printDebug(LIGHT_DEBUG, "Time", testCase_time) ;
			printDebug(LIGHT_DEBUG, "Speed", testCase_bw) ;
			testCase_time_best = min(testCase_time_best, testCase_time) ;

			// check testCase results
			if (itry == 0)
			{
				checkAllProcGridL1Err(INNER_POINTS, Wgrid, Rgrid, MAX_ERR_FLOAT) ;
			}

		}

		UpdatePressureGB = Myfloat(5 * nGridPointGlob*sizeof(Myfloat)/testCase_time_best/1.e9) ;
		UpdatePressureGpoint = Myfloat(nGridPointGlob/testCase_time_best/1.e9) ;
		printInfo(MASTER, " Best achieved GByte/s", UpdatePressureGB) ;
		printInfo(MASTER, " Best achieved GPoint/s", UpdatePressureGpoint) ;
	}

	{
		//-----------------------------
		// applyBoundaryCondition
		//-----------------------------
		print_blank() ;
		string caseName = testCaseName + "ApplyBoundaryCondition" ;
		printInfo(MASTER, " * Case", caseName) ;

		// set unit grid
		Wgrid.defineUnitGrid() ;
		Rgrid.defineUnitGrid() ;

		// fill Wgrid with sine function
		// note: only sine function can be used since it is compatible with the
		// boundary condition of type BOUND_COND_ANTI_MIRROR
		Wgrid.fill(ALL_POINTS, FUNC_SINE, FUNC_SINE, FUNC_SINE, PI, PI, PI, 1.0) ;

		// set all edge halos to zeros
		if (Wgrid.dim >= DIM1)
		{
			if (Wgrid.getNeighbourProc(I1HALO1) == MPI_PROC_NULL) Wgrid.fill(I1HALO1, -1) ;
			if (Wgrid.getNeighbourProc(I1HALO2) == MPI_PROC_NULL) Wgrid.fill(I1HALO2, -1) ;
		}
		if (Wgrid.dim >= DIM2)
		{
			if (Wgrid.getNeighbourProc(I2HALO1) == MPI_PROC_NULL) Wgrid.fill(I2HALO1, -1) ;
			if (Wgrid.getNeighbourProc(I2HALO2) == MPI_PROC_NULL) Wgrid.fill(I2HALO2, -1) ;
		}
		if (Wgrid.dim >= DIM3)
		{
			if (Wgrid.getNeighbourProc(I3HALO1) == MPI_PROC_NULL) Wgrid.fill(I3HALO1, -1) ;
			if (Wgrid.getNeighbourProc(I3HALO2) == MPI_PROC_NULL) Wgrid.fill(I3HALO2, -1) ;
		}

		// update all halos
		Wgrid.exchangeHalos(MPI_COMM_MODE_SENDRECV) ;

		// fill gridRef with expected values
		Rgrid.fill(ALL_POINTS, FUNC_SINE, FUNC_SINE, FUNC_SINE, PI, PI, PI, 1.0) ;
		Rgrid.write(caseName+"Ref") ;

		double testCase_time_best = FLT_MAX ;

		Myint64 nGridPointBoundLoc = Wgrid.getNumberOfGridPointBoundaryCondition(BOUND_COND_ANTI_MIRROR) ;
		Myint64 nGridPointBound = 0 ;
		MPI_Reduce(&nGridPointBoundLoc, &nGridPointBound, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
		printInfo(MASTER, " Boundary pts (Loc.)", nGridPointBoundLoc) ;
		printInfo(MASTER, " Boundary pts (Glob.)", nGridPointBound) ;
		printInfo(MASTER, " % Boundary vs Inner", (Myfloat) nGridPointBound/nGridPointGlob*100) ;

		for (Myint itry = 0; itry < ntry; itry++)
		{
			double t0 = MPI_Wtime() ;
			Wgrid.applyBoundaryCondition(BOUND_COND_ANTI_MIRROR) ;
			// wait all process completed before ending time
			MPI_Barrier(MPI_COMM_WORLD) ;
			double t1 = MPI_Wtime() ;

			double testCase_time = t1-t0 ;
			Myfloat testCase_bw = 2 * nGridPointGlob*sizeof(Myfloat)/testCase_time/1e9 ;
			printDebug(LIGHT_DEBUG, "Time", testCase_time) ;
			printDebug(LIGHT_DEBUG, "Speed", testCase_bw) ;
			testCase_time_best = min(testCase_time_best, testCase_time) ;

			// check testCase results
			if (itry == 0)
			{
				checkAllProcGridL1Err(INNER_POINTS, Wgrid, Rgrid, MAX_ERR_FLOAT) ;
				Wgrid.write(caseName+"W") ;
			}

		}

		// we multiply by 2 the number of point to get an approximate of the bytes
		// this number is not exact
		ApplyBoundaryConditionGB = Myfloat(2 * nGridPointGlob*sizeof(Myfloat)/testCase_time_best/1.e9) ;
		printInfo(MASTER, " Best achieved GByte/s", ApplyBoundaryConditionGB) ;

		// Number of point is exact
		ApplyBoundaryConditionGpoint = Myfloat(nGridPointGlob/testCase_time_best/1.e9) ;
		printInfo(MASTER, " Best achieved GPoint/s", ApplyBoundaryConditionGpoint) ;
	}

	// log perf
	if (myid_world == 0)
	{
		perfLogFile
		// 10, 11, 12, 13
		<< FillGB << " " << FillGpoint << " " << MaxErrGB << " " << MaxErrGpoint << " "
		// 14, 15, 16, 17
		<< L1ErrGB << " " << L1ErrGpoint << " " << GetSumAbsGB << " " << GetSumAbsGpoint << " "
		// 18, 19, 20, 21
		<< GetSumAbsDiffGB << " " << GetSumAbsDiffGpoint << " " << GetMaxGB << " " << GetMaxGpoint << " "
		// 22, 23, 24, 25
		<< GetMinGB << " " << GetMinGpoint << " " << UpdatePressureGB << " " << UpdatePressureGpoint << " "
		// 26, 27
		<< ApplyBoundaryConditionGB << " " << ApplyBoundaryConditionGpoint << " "
		<< "\n" ;
	}

	this->finalize() ;

	printDebug(MID_DEBUG, "Out TestCase_Grid::run") ;
	return(RTN_CODE_OK) ;
}

} // namespace hpcscan

