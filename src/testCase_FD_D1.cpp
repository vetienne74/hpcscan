
//-------------------------------------------------------------------------------------------------------
// test cases to measure the finite-difference computation bandwidth
// first derivatice operator D1 (i.e. d/dx)
// baseline kernel - no optimization
//-------------------------------------------------------------------------------------------------------

#include "testCase_FD_D1.h"

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

TestCase_FD_D1::TestCase_FD_D1(void)
{
	testCaseName    = "FD_D1" ;
	testCaseVersion = "Standard implementation" ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code TestCase_FD_D1::run(void)
{

	printDebug(MID_DEBUG, "In TestCase_FD_D1::run") ;

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
	Myfloat D1Axis1Gflop=0, D1Axis2Gflop=0  , D1Axis3Gflop=0  , D1LaplaGflop=0 ;
	Myfloat D1Axis1GpointEff=0, D1Axis2GpointEff=0, D1Axis3GpointEff=0, D1LaplaGpointEff=0 ;
	Myfloat D1Axis1GpointFD=0, D1Axis2GpointFD=0, D1Axis3GpointFD=0, D1LaplaGpointFD=0 ;
	Myfloat D1Axis1GB=0, D1Axis2GB=0, D1Axis3GB=0, D1LaplaGB=0 ;

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
			//TODO Ugrid.FD_D1_N1(INNER_POINTS, Wgrid, fdOrder) ;
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

		//TODO Myint nOpPerPoint   = Ugrid.getFlopPerPtFD_D1(fdOrder) ;
		//TODO Myint nPtPerStencil = Ugrid.getPtPerStencilFD_D1(fdOrder) ;
		Myint nOpPerPoint   = Ugrid.getFlopPerPtFD_D2(fdOrder) ;
		Myint nPtPerStencil = Ugrid.getPtPerStencilFD_D2(fdOrder) ;

        Myint nMemOpPerPoint = nPtPerStencil + 1 ; // + 1 store
		printInfo(MASTER, " #Flop per point", nOpPerPoint) ;
		printInfo(MASTER, " #Point in stencil", nPtPerStencil) ;

		double timeInFD = testCase_time_best-testCase_time_com ;
		D1Axis1Gflop     = nGridPointGlob/timeInFD/1.e9 * nOpPerPoint;
		D1Axis1GpointEff = nGridPointGlob/testCase_time_best/1.e9 ;
		D1Axis1GpointFD  = nGridPointGlob/timeInFD/1.e9 ;
		D1Axis1GB        = nGridPointGlob/timeInFD/1.e9 * nMemOpPerPoint * sizeof(Myfloat) ;

		printInfo(MASTER, " Best GFlop/s in FD" ,   D1Axis1Gflop) ;
		printInfo(MASTER, " Best Gpoint/s eff." ,   D1Axis1GpointEff) ;
		printInfo(MASTER, " Best Gpoint/s in FD",   D1Axis1GpointFD) ;
		printInfo(MASTER, " Best Apparent BW GB/s", D1Axis1GB) ;

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
			//TODO Ugrid.FD_D1_N2(INNER_POINTS, Wgrid, fdOrder) ;
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

		//TODO Myint nOpPerPoint   = Ugrid.getFlopPerPtFD_D1(fdOrder) ;
		//TODO Myint nPtPerStencil = Ugrid.getPtPerStencilFD_D1(fdOrder) ;
		Myint nOpPerPoint   = Ugrid.getFlopPerPtFD_D2(fdOrder) ;
		Myint nPtPerStencil = Ugrid.getPtPerStencilFD_D2(fdOrder) ;
        Myint nMemOpPerPoint = nPtPerStencil + 1 ; // + 1 store
		printInfo(MASTER, " #Flop per point", nOpPerPoint) ;
		printInfo(MASTER, " #Point in stencil", nPtPerStencil) ;

		double timeInFD = testCase_time_best-testCase_time_com ;
		D1Axis2Gflop     = nGridPointGlob/timeInFD/1.e9 * nOpPerPoint;
		D1Axis2GpointEff = nGridPointGlob/testCase_time_best/1.e9 ;
		D1Axis2GpointFD  = nGridPointGlob/timeInFD/1.e9 ;
		D1Axis2GB        = nGridPointGlob/timeInFD/1.e9 * nMemOpPerPoint * sizeof(Myfloat) ;

		printInfo(MASTER, " Best GFlop/s in FD" ,   D1Axis2Gflop) ;
		printInfo(MASTER, " Best Gpoint/s eff." ,   D1Axis2GpointEff) ;
		printInfo(MASTER, " Best Gpoint/s in FD",   D1Axis2GpointFD) ;
		printInfo(MASTER, " Best Apparent BW GB/s", D1Axis2GB) ;

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
			//TODO Ugrid.FD_D1_N3(INNER_POINTS, Wgrid, fdOrder) ;
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

		//TODO Myint nOpPerPoint    = Ugrid.getFlopPerPtFD_D1(fdOrder) ;
		//TODO Myint nPtPerStencil  = Ugrid.getPtPerStencilFD_D1(fdOrder) ;
        Myint nOpPerPoint    = Ugrid.getFlopPerPtFD_D2(fdOrder) ;
		Myint nPtPerStencil  = Ugrid.getPtPerStencilFD_D2(fdOrder) ;
		Myint nMemOpPerPoint = nPtPerStencil + 1 ; // + 1 store
		printInfo(MASTER, " #Flop per point", nOpPerPoint) ;
		printInfo(MASTER, " #Point in stencil", nPtPerStencil) ;

		double timeInFD = testCase_time_best-testCase_time_com ;
		D1Axis3Gflop     = nGridPointGlob/timeInFD/1.e9 * nOpPerPoint;
		D1Axis3GpointEff = nGridPointGlob/testCase_time_best/1.e9 ;
		D1Axis3GpointFD  = nGridPointGlob/timeInFD/1.e9 ;
		D1Axis3GB        = nGridPointGlob/timeInFD/1.e9 * nMemOpPerPoint * sizeof(Myfloat) ;

		printInfo(MASTER, " Best GFlop/s in FD" ,   D1Axis3Gflop) ;
		printInfo(MASTER, " Best Gpoint/s eff." ,   D1Axis3GpointEff) ;
		printInfo(MASTER, " Best Gpoint/s in FD",   D1Axis3GpointFD) ;
		printInfo(MASTER, " Best Apparent BW GB/s", D1Axis3GB) ;

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
		D1LaplaGflop     = nGridPointGlob/timeInFD/1.e9 * nOpPerPoint;
		D1LaplaGpointEff = nGridPointGlob/testCase_time_best/1.e9 ;
		D1LaplaGpointFD  = nGridPointGlob/timeInFD/1.e9 ;
		D1LaplaGB        = nGridPointGlob/timeInFD/1.e9 * nMemOpPerPoint * sizeof(Myfloat) ;

		printInfo(MASTER, " Best GFlop/s in FD" ,   D1LaplaGflop) ;
		printInfo(MASTER, " Best Gpoint/s eff." ,   D1LaplaGpointEff) ;
		printInfo(MASTER, " Best Gpoint/s in FD",   D1LaplaGpointFD) ;
		printInfo(MASTER, " Best Apparent BW GB/s", D1LaplaGB) ;

		if (nGridPointHaloGlob > 0)
		{
			printInfo(MASTER, " % MPI Comm\t", Myfloat(testCase_time_com / testCase_time_best * 100)) ;
			printInfo(MASTER, " MPI Comm. BW GB/s", Myfloat(nGridPointHaloGlob/testCase_time_com
					/ 1.e9 * sizeof(Myfloat))) ;
		}
	}

	// log perf
	if (myMpiRank == 0)
	{
		perfLogFile
		// 10, 11, 12, 13
		<< D1Axis1Gflop << " " << D1Axis1GpointFD << " " << D1Axis1GpointEff << " " << D1Axis1GB << " "
		// 14, 15, 16, 17
		<< D1Axis2Gflop << " " << D1Axis2GpointFD << " " << D1Axis2GpointEff << " " << D1Axis2GB << " "
		// 18, 19, 20, 21
		<< D1Axis3Gflop << " " << D1Axis3GpointFD << " " << D1Axis3GpointEff << " " << D1Axis3GB << " "
		// 22, 23, 24, 25
		<< D1LaplaGflop << " " << D1LaplaGpointFD << " " << D1LaplaGpointEff << " " << D1LaplaGB << " "

		// cache block sizes
		<< Config::Instance()->cb1 << " " << Config::Instance()->cb2 << " " << Config::Instance()->cb3
		<< "\n" ;
	}

	this->finalize() ;

	printDebug(MID_DEBUG, "Out TestCase_FD_D1::run") ;
	return(RTN_CODE_OK) ;
}

} // namespace hpcscan

