
//-------------------------------------------------------------------------------------------------------
// test cases to measure the finite-difference computation bandwidth
// baseline kernel - no optimization
//-------------------------------------------------------------------------------------------------------

#include "testCase_FD_D2.h"

#include <algorithm>
#include <cfloat>

#include <mpi.h>
#include <omp.h>

#include "config.h"
#include "constant.h"
#include "fdm.h"
#include "global.h"
#include "grid_Factory.h"
#include "output_report.h"

namespace hpcscan {

TestCase_FD_D2::TestCase_FD_D2(void)
{
	testCaseName    = "FD_D2" ;
	testCaseVersion = "Standard implementation" ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code TestCase_FD_D2::run(void)
{

	printDebug(MID_DEBUG, "In TestCase_FD_D2::run") ;

	if (this->initialize() == RTN_CODE_KO) return (RTN_CODE_KO) ;

	const string gridMode = Config::Instance()->testMode ;
	auto Ugrid2 = Grid_Factory::create(gridMode, GRID_LOCAL) ;
	auto Wgrid2 = Grid_Factory::create(gridMode, GRID_LOCAL) ;
	auto Rgrid2 = Grid_Factory::create(gridMode, GRID_LOCAL) ;
	Grid &Ugrid = *Ugrid2 ; // input grid
	Grid &Wgrid = *Wgrid2 ; // output of FD computations
	Grid &Rgrid = *Rgrid2 ; // reference grid
	Ugrid.initializeGrid() ;
	Wgrid.initializeGrid() ;
	Rgrid.initializeGrid() ;

	Ugrid.info() ;

	// Number of points exchanged in Halos
	Myint64 nGridPointHaloLoc = Ugrid.getNumberOfGridPointCommHalo(GRID_LOCAL) ;
	Myint64 nGridPointHaloGlob = Ugrid.getNumberOfGridPointCommHalo(GRID_GLOBAL) ;
	print_blank() ;
	printInfo(MASTER, " Halos Com Local (Pts)", nGridPointHaloLoc) ;
	printInfo(MASTER, " Halos Com Global (Pts)", nGridPointHaloGlob) ;
	printInfo(MASTER, " Halos Com Global (MB)", Myfloat(nGridPointHaloGlob*sizeof(Myfloat)/1.e6)) ;

	const Myfloat64 a1 = Config::Instance()->param1 ;
	const Myfloat64 a2 = Config::Instance()->param2 ;
	const Myfloat64 a3 = Config::Instance()->param3 ;
	const Myfloat64 a4 = Config::Instance()->param4 ;

	//............................................
	// fill U grid with
	// sin(a1*x1) * sin(a2*x2) * sin(a3*x3) * a4
	//............................................
	Ugrid.fill(ALL_POINTS, FUNC_SINE, FUNC_SINE, FUNC_SINE, a1, a2, a3, a4) ;
	Ugrid.write(testCaseName+"U") ;

	Myfloat * const u = Ugrid.grid_3d ;
	Myfloat * const w = Wgrid.grid_3d ;

	double testCase_time_best, testCase_time_com ;

	Myint ntry = Config::Instance()->ntry ;

	// define loop
	Myint64 i1Start, i1End, i2Start, i2End, i3Start, i3End ;
	Ugrid.getGridIndex(INNER_POINTS, &i1Start, &i1End, &i2Start, &i2End, &i3Start, &i3End) ;

	Myint64 nGridPointLoc  = Ugrid.getNumberOfGridPoint(GRID_LOCAL, INNER_POINTS) ;
	Myint64 nGridPointGlob = Ugrid.getNumberOfGridPoint(GRID_GLOBAL, INNER_POINTS) ;
	printInfo(MASTER, " Grid Inner Loc. (Pts)", nGridPointLoc) ;
	printInfo(MASTER, " Grid Inner Glob. (Pts)", nGridPointGlob) ;
	printInfo(MASTER, " Grid Inner Glob. (MB)",  Myfloat(nGridPointGlob*sizeof(Myfloat)/1.e6)) ;
	print_blank() ;

	const Myint fdOrder = 2 * Ugrid.haloWidth ;
	printInfo(MASTER, " FD space order\t", fdOrder) ;

	const Myfloat maxErr = 1.e-3 ;
	printInfo(MASTER, " Max allowed error", maxErr) ;

	// for perf log
	Myfloat D2Axis1Gflop=0, D2Axis2Gflop=0  , D2Axis3Gflop=0  , D2LaplaGflop=0 ;
	Myfloat D2Axis1GpointEff=0, D2Axis2GpointEff=0, D2Axis3GpointEff=0, D2LaplaGpointEff=0 ;
	Myfloat D2Axis1GpointFD=0, D2Axis2GpointFD=0, D2Axis3GpointFD=0, D2LaplaGpointFD=0 ;
	Myfloat D2Axis1GB=0, D2Axis2GB=0, D2Axis3GB=0, D2LaplaGB=0 ;

	if (Config::Instance()->dim >= DIM1)
	{
		//====================================================
		// W = d2(U)/dx1x1
		// reference solution:
		// -a1^2 * sin(a1*x1) * sin(a2*x2) * sin(a3*x3) * a4
		//====================================================

		print_blank() ;
		string caseName = testCaseName + "Axis1" ;
		printInfo(MASTER, " * Case", caseName) ;

		Myfloat64 ampRef = - a1 * a1 * a4 ;
		Rgrid.fill(ALL_POINTS, FUNC_SINE, FUNC_SINE, FUNC_SINE, a1, a2, a3, ampRef) ;
		Rgrid.write(testCaseName+"R") ;

		testCase_time_best = FLT_MAX ;
		Wgrid.fill(ALL_POINTS, 0.0) ;

		for (Myint itry = 0; itry < ntry; itry++)
		{
			// start timer
			double t0 = MPI_Wtime() ;

			// exchange halos
			Ugrid.exchangeHalos(MPI_COMM_MODE_SENDRECV) ;
			double t1 = MPI_Wtime() ;

			// compute FD
			Ugrid.FD_D2_N1(INNER_POINTS, Wgrid, fdOrder) ;

			// wait all process completed computations before ending time
			MPI_Barrier(MPI_COMM_WORLD) ;
			double t2 = MPI_Wtime() ;

			double testCase_time = t2-t0 ;
			printDebug(LIGHT_DEBUG, "testCase time", testCase_time) ;
			if (testCase_time < testCase_time_best)
			{
				testCase_time_best = testCase_time ;
				testCase_time_com = t1-t0 ;
			}

			// check testCase results
			if (itry == 0)
			{
				Wgrid.write(caseName+"W") ;
				checkGridL1Err(INNER_POINTS, Wgrid, Rgrid, maxErr) ;
			}
		}

		Myint nOpPerPoint   = Ugrid.getFlopPerPtFD_D2(fdOrder) ;
		Myint nPtPerStencil = Ugrid.getPtPerStencilFD_D2(fdOrder) ;
		Myint nMemOpPerPoint = nPtPerStencil + 1 ; // + 1 store
		printInfo(MASTER, " #Flop per point", nOpPerPoint) ;
		printInfo(MASTER, " #Point in stencil", nPtPerStencil) ;

		double timeInFD = testCase_time_best-testCase_time_com ;
		D2Axis1Gflop     = nGridPointGlob/timeInFD/1.e9 * nOpPerPoint;
		D2Axis1GpointEff = nGridPointGlob/testCase_time_best/1.e9 ;
		D2Axis1GpointFD  = nGridPointGlob/timeInFD/1.e9 ;
		D2Axis1GB        = nGridPointGlob/timeInFD/1.e9 * nMemOpPerPoint * sizeof(Myfloat) ;

		printInfo(MASTER, " Best GFlop/s in FD" ,   D2Axis1Gflop) ;
		printInfo(MASTER, " Best Gpoint/s eff." ,   D2Axis1GpointEff) ;
		printInfo(MASTER, " Best Gpoint/s in FD",   D2Axis1GpointFD) ;
		printInfo(MASTER, " Best Apparent BW GB/s", D2Axis1GB) ;

		if (nGridPointHaloGlob > 0)
		{
			printInfo(MASTER, " % MPI Comm\t", Myfloat(testCase_time_com / testCase_time_best * 100)) ;
			printInfo(MASTER, " MPI Comm. BW GB/s", Myfloat(nGridPointHaloGlob/testCase_time_com
					/ 1.e9 * sizeof(Myfloat))) ;
		}
	}

	if (Config::Instance()->dim >= DIM2)
	{
		//====================================================
		// W = d2(U)/dx2x2
		// reference solution:
		// -a2^2 * sin(a1*x1) * sin(a2*x2) * sin(a3*x3) * a4
		//====================================================

		print_blank() ;
		string caseName = testCaseName + "Axis2" ;
		printInfo(MASTER, " * Case", caseName) ;

		Myfloat64 ampRef = - a2 * a2 * a4 ;
		Rgrid.fill(ALL_POINTS, FUNC_SINE, FUNC_SINE, FUNC_SINE, a1, a2, a3, ampRef) ;
		Rgrid.write(testCaseName+"R") ;

		testCase_time_best = FLT_MAX ;
		Wgrid.fill(ALL_POINTS, 0.0) ;

		for (Myint itry = 0; itry < ntry; itry++)
		{
			// start timer
			double t0 = MPI_Wtime() ;

			// exchange halos
			Ugrid.exchangeHalos(MPI_COMM_MODE_SENDRECV) ;
			double t1 = MPI_Wtime() ;

			// compute FD
			Ugrid.FD_D2_N2(INNER_POINTS, Wgrid, fdOrder) ;

			// wait all process completed computations before ending time
			MPI_Barrier(MPI_COMM_WORLD) ;
			double t2 = MPI_Wtime() ;

			double testCase_time = t2-t0 ;
			printDebug(LIGHT_DEBUG, "testCase time", testCase_time) ;
			if (testCase_time < testCase_time_best)
			{
				testCase_time_best = testCase_time ;
				testCase_time_com = t1-t0 ;
			}

			// check testCase results
			if (itry == 0)
			{
				Wgrid.write(caseName+"W") ;
				checkGridL1Err(INNER_POINTS, Wgrid, Rgrid, maxErr) ;
			}
		}

		Myint nOpPerPoint   = Ugrid.getFlopPerPtFD_D2(fdOrder) ;
		Myint nPtPerStencil = Ugrid.getPtPerStencilFD_D2(fdOrder) ;
		Myint nMemOpPerPoint = nPtPerStencil + 1 ; // + 1 store
		printInfo(MASTER, " #Flop per point", nOpPerPoint) ;
		printInfo(MASTER, " #Point in stencil", nPtPerStencil) ;

		double timeInFD = testCase_time_best-testCase_time_com ;
		D2Axis2Gflop     = nGridPointGlob/timeInFD/1.e9 * nOpPerPoint;
		D2Axis2GpointEff = nGridPointGlob/testCase_time_best/1.e9 ;
		D2Axis2GpointFD  = nGridPointGlob/timeInFD/1.e9 ;
		D2Axis2GB        = nGridPointGlob/timeInFD/1.e9 * nMemOpPerPoint * sizeof(Myfloat) ;

		printInfo(MASTER, " Best GFlop/s in FD" ,   D2Axis2Gflop) ;
		printInfo(MASTER, " Best Gpoint/s eff." ,   D2Axis2GpointEff) ;
		printInfo(MASTER, " Best Gpoint/s in FD",   D2Axis2GpointFD) ;
		printInfo(MASTER, " Best Apparent BW GB/s", D2Axis2GB) ;

		if (nGridPointHaloGlob > 0)
		{
			printInfo(MASTER, " % MPI Comm\t", Myfloat(testCase_time_com / testCase_time_best * 100)) ;
			printInfo(MASTER, " MPI Comm. BW GB/s", Myfloat(nGridPointHaloGlob/testCase_time_com
					/ 1.e9 * sizeof(Myfloat))) ;
		}
	}

	if (Config::Instance()->dim >= DIM3)
	{
		//====================================================
		// W = d2(U)/dx3x3
		// reference solution:
		// -a3^2 * sin(a1*x1) * sin(a2*x2) * sin(a3*x3) * a4
		//====================================================

		print_blank() ;
		string caseName = testCaseName + "Axis3" ;
		printInfo(MASTER, " * Case", caseName) ;

		Myfloat64 ampRef = - a3 * a3 * a4 ;
		Rgrid.fill(ALL_POINTS, FUNC_SINE, FUNC_SINE, FUNC_SINE, a1, a2, a3, ampRef) ;
		Rgrid.write(testCaseName+"R") ;

		testCase_time_best = FLT_MAX ;
		Wgrid.fill(ALL_POINTS, 0.0) ;

		for (Myint itry = 0; itry < ntry; itry++)
		{
			// start timer
			double t0 = MPI_Wtime() ;

			// exchange halos
			Ugrid.exchangeHalos(MPI_COMM_MODE_SENDRECV) ;
			double t1 = MPI_Wtime() ;

			// compute FD
			Ugrid.FD_D2_N3(INNER_POINTS, Wgrid, fdOrder) ;

			// wait all process completed computations before ending time
			MPI_Barrier(MPI_COMM_WORLD) ;
			double t2 = MPI_Wtime() ;

			double testCase_time = t2-t0 ;
			printDebug(LIGHT_DEBUG, "testCase time", testCase_time) ;
			if (testCase_time < testCase_time_best)
			{
				testCase_time_best = testCase_time ;
				testCase_time_com = t1-t0 ;
			}

			// check testCase results
			if (itry == 0)
			{
				Wgrid.write(caseName+"W") ;
				checkGridL1Err(INNER_POINTS, Wgrid, Rgrid, maxErr) ;
			}
		}

		Myint nOpPerPoint    = Ugrid.getFlopPerPtFD_D2(fdOrder) ;
		Myint nPtPerStencil  = Ugrid.getPtPerStencilFD_D2(fdOrder) ;
		Myint nMemOpPerPoint = nPtPerStencil + 1 ; // + 1 store
		printInfo(MASTER, " #Flop per point", nOpPerPoint) ;
		printInfo(MASTER, " #Point in stencil", nPtPerStencil) ;

		double timeInFD = testCase_time_best-testCase_time_com ;
		D2Axis3Gflop     = nGridPointGlob/timeInFD/1.e9 * nOpPerPoint;
		D2Axis3GpointEff = nGridPointGlob/testCase_time_best/1.e9 ;
		D2Axis3GpointFD  = nGridPointGlob/timeInFD/1.e9 ;
		D2Axis3GB        = nGridPointGlob/timeInFD/1.e9 * nMemOpPerPoint * sizeof(Myfloat) ;

		printInfo(MASTER, " Best GFlop/s in FD" ,   D2Axis3Gflop) ;
		printInfo(MASTER, " Best Gpoint/s eff." ,   D2Axis3GpointEff) ;
		printInfo(MASTER, " Best Gpoint/s in FD",   D2Axis3GpointFD) ;
		printInfo(MASTER, " Best Apparent BW GB/s", D2Axis3GB) ;

		if (nGridPointHaloGlob > 0)
		{
			printInfo(MASTER, " % MPI Comm\t", Myfloat(testCase_time_com / testCase_time_best * 100)) ;
			printInfo(MASTER, " MPI Comm. BW GB/s", Myfloat(nGridPointHaloGlob/testCase_time_com
					/ 1.e9 * sizeof(Myfloat))) ;
		}
	}

	{
		//===========================================================
		// W = Laplacian(U)
		//===========================================================

		string caseName ;
		Myfloat64 ampRef= 0.0 ;

		if (Config::Instance()->dim == DIM1)
		{
			//===========================================================
			// W = Laplacian(U) for 1D
			// For 1D same as same W = d2(U)/dx1x1
			//===========================================================

			caseName = testCaseName + "Laplacian1D" ;
			ampRef = - a1 * a1 * a4 ;
		}
		else if (Config::Instance()->dim == DIM2)
		{
			//===========================================================
			// W = Laplacian(U) for 2D
			// reference solution:
			// -(a1^2+a2^2) * sin(a1*x1) * sin(a2*x2) * sin(a3*x3) * a4
			//===========================================================

			caseName = testCaseName + "Laplacian2D" ;
			ampRef = - (a1 * a1 + a2 * a2) * a4 ;
		}
		else if (Config::Instance()->dim == DIM3)
		{
			//================================================================
			// W = Laplacian(U) for 3D
			// reference solution:
			// -(a1^2+a2^2+a3^2) * sin(a1*x1) * sin(a2*x2) * sin(a3*x3) * a4
			//================================================================

			caseName = testCaseName + "Laplacian3D" ;
			ampRef = - (a1 * a1 + a2 * a2 + a3 * a3) * a4 ;
		}

		print_blank() ;
		printInfo(MASTER, " * Case", caseName) ;
		Rgrid.fill(ALL_POINTS, FUNC_SINE, FUNC_SINE, FUNC_SINE, a1, a2, a3, ampRef) ;
		Rgrid.write(testCaseName+"R") ;

		testCase_time_best = FLT_MAX ;
		Wgrid.fill(ALL_POINTS, 0.0) ;

		for (Myint itry = 0; itry < ntry; itry++)
		{
			// start timer
			double t0 = MPI_Wtime() ;

			// exchange halos
			Ugrid.exchangeHalos(MPI_COMM_MODE_SENDRECV) ;
			double t1 = MPI_Wtime() ;

			// compute FD
			Ugrid.FD_LAPLACIAN(INNER_POINTS, Wgrid, fdOrder) ;

			// wait all process completed computations before ending time
			MPI_Barrier(MPI_COMM_WORLD) ;
			double t2 = MPI_Wtime() ;

			double testCase_time = t2-t0 ;
			printDebug(LIGHT_DEBUG, "testCase time", testCase_time) ;
			if (testCase_time < testCase_time_best)
			{
				testCase_time_best = testCase_time ;
				testCase_time_com = t1-t0 ;
			}

			// check testCase results
			if (itry == 0)
			{
				Wgrid.write(caseName+"W") ;
				checkGridL1Err(INNER_POINTS, Wgrid, Rgrid, maxErr) ;
			}
		}

		Myint nPtPerStencil = Ugrid.getPtPerStencilFD_LAPLACIAN(fdOrder) ;
		Myint nMemOpPerPoint = nPtPerStencil + 1 ; // + 1 store
		Myint nOpPerPoint = Ugrid.getFlopPerPtFD_LAPLACIAN(fdOrder) ;

		printInfo(MASTER, " #Flop per point", nOpPerPoint) ;
		printInfo(MASTER, " #Point in stencil", nPtPerStencil) ;

		double timeInFD = testCase_time_best-testCase_time_com ;
		D2LaplaGflop     = nGridPointGlob/timeInFD/1.e9 * nOpPerPoint;
		D2LaplaGpointEff = nGridPointGlob/testCase_time_best/1.e9 ;
		D2LaplaGpointFD  = nGridPointGlob/timeInFD/1.e9 ;
		D2LaplaGB        = nGridPointGlob/timeInFD/1.e9 * nMemOpPerPoint * sizeof(Myfloat) ;

		printInfo(MASTER, " Best GFlop/s in FD" ,   D2LaplaGflop) ;
		printInfo(MASTER, " Best Gpoint/s eff." ,   D2LaplaGpointEff) ;
		printInfo(MASTER, " Best Gpoint/s in FD",   D2LaplaGpointFD) ;
		printInfo(MASTER, " Best Apparent BW GB/s", D2LaplaGB) ;

		if (nGridPointHaloGlob > 0)
		{
			printInfo(MASTER, " % MPI Comm\t", Myfloat(testCase_time_com / testCase_time_best * 100)) ;
			printInfo(MASTER, " MPI Comm. BW GB/s", Myfloat(nGridPointHaloGlob/testCase_time_com
					/ 1.e9 * sizeof(Myfloat))) ;
		}
	}

	// log perf
	if (myid_world == 0)
	{
		perfLogFile
		// 10, 11, 12, 13
		<< D2Axis1Gflop << " " << D2Axis1GpointFD << " " << D2Axis1GpointEff << " " << D2Axis1GB << " "
		// 14, 15, 16, 17
		<< D2Axis2Gflop << " " << D2Axis2GpointFD << " " << D2Axis2GpointEff << " " << D2Axis2GB << " "
		// 18, 19, 20, 21
		<< D2Axis3Gflop << " " << D2Axis3GpointFD << " " << D2Axis3GpointEff << " " << D2Axis3GB << " "
		// 22, 23, 24, 25
		<< D2LaplaGflop << " " << D2LaplaGpointFD << " " << D2LaplaGpointEff << " " << D2LaplaGB << " "
		<< "\n" ;
	}

	this->finalize() ;

	printDebug(MID_DEBUG, "Out TestCase_FD_D2::run") ;
	return(RTN_CODE_OK) ;
}

} // namespace hpcscan

