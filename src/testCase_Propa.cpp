
//-------------------------------------------------------------------------------------------------------
// Test cases to measure the acoustic propagator bandwidth
//-------------------------------------------------------------------------------------------------------

#include "testCase_Propa.h"

#include <algorithm>
#include <cfloat>

#include <mpi.h>
#include <omp.h>

#include "config.h"
#include "fdm.h"
#include "global.h"
#include "grid_Factory.h"
#include "output_report.h"
#include "propagator_Factory.h"

namespace hpcscan {

TestCase_Propa::TestCase_Propa(void)
{
	testCaseName    = "Propa" ;
	testCaseVersion = "Standard implementation" ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code TestCase_Propa::run(void)
{

	printDebug(MID_DEBUG, "In TestCase_Propa::run") ;

	if (this->initialize() == RTN_CODE_KO) return (RTN_CODE_KO) ;

	// instantiate propagator
	const string propagator = Config::Instance()->propagator ;
	auto propa = Propagator_Factory::create(propagator) ;
	if (propa == nullptr) return(RTN_CODE_KO) ;

	// initialize propagator
	Rtn_code rtnCode = propa->initialize(EIGEN_MODE) ;
	if (rtnCode != RTN_CODE_OK)
	{
		printError("In TestCase_Propa::run, propa.initialize() not Ok") ;
		return(RTN_CODE_OK) ;
	}

	// print propagator info
	propa->info() ;

	// get pointers on propagators grids
	auto prnGrid  = propa->prnGrid ;
	auto prcGrid  = propa->prcGrid ;

	// allocate and initialize reference grid
	const string gridMode = Config::Instance()->testMode ;
	auto refGrid = Grid_Factory::create(gridMode, GRID_LOCAL) ;
	refGrid->initializeGrid() ;
	refGrid->defineUnitGrid() ;

	if (myMpiRank == 0) prnGrid->info() ;

	// Number of points exchanged in Halos
	Myint64 nGridPointHaloLoc  = prnGrid->getNumberOfGridPointCommHalo(GRID_LOCAL) ;
	Myint64 nGridPointHaloGlob = prnGrid->getNumberOfGridPointCommHalo(GRID_GLOBAL) ;
	print_blank() ;
	printInfo(MASTER, " Halos Com Local (Pts)", nGridPointHaloLoc) ;
	printInfo(MASTER, " Halos Com Global (Pts)", nGridPointHaloGlob) ;
	printInfo(MASTER, " Halos Com Global (MB)", Myfloat(nGridPointHaloGlob*sizeof(Myfloat)/1.e6)) ;

	double testCase_time_best = 0.0, testCase_time_com = 0.0 ;

	Myint ntry = Config::Instance()->ntry ;

	Myint64 nGridPointLoc  = prnGrid->getNumberOfGridPoint(GRID_LOCAL, INNER_POINTS) ;
	Myint64 nGridPointGlob = prnGrid->getNumberOfGridPoint(GRID_GLOBAL, INNER_POINTS) ;
	printInfo(MASTER, " Grid Inner Loc. (Pts)", nGridPointLoc) ;
	printInfo(MASTER, " Grid Inner Glob. (Pts)", nGridPointGlob) ;
	printInfo(MASTER, " Grid Inner Glob. (MB)",  Myfloat(nGridPointGlob*sizeof(Myfloat)/1.e6)) ;
	print_blank() ;

	// get parameters from propagator
	const Myint fdOrder = propa->fdOrder ;
	const Myint nt      = propa->nt ;
	const Myfloat dt    = propa->dt ;
	const Myint snapInc = propa->snapInc ;

	const Myfloat maxErr = 2.e-2 ;
	printInfo(MASTER, " Max allowed error", maxErr) ;

	Myint nPtPerStencil  = prcGrid->getPtPerStencilFD_LAPLACIAN(fdOrder) ;
	Myint nMemOpPerPoint = nPtPerStencil + 4 ; // + 1 store + 3 load
	Myint nOpPerPoint    = prcGrid->getFlopPerPtFD_LAPLACIAN(fdOrder) + 4 ; // + 1 ADD + 1 SUB + 2 MUL

	// for perf log
	Myfloat propaGflop=0, propaGpointEff=0, propaGpointFD=0, propaGB=0 ;

	print_blank() ;
	string caseName = testCaseName + "EigenMode" ;
	printInfo(MASTER, " * Case", caseName) ;

	testCase_time_best = FLT_MAX ;

	Myfloat errTestCase, sum1, sum2 ;
	Myint   ntCheck ;

	for (Myint itry = 0; itry < ntry; itry++)
	{
		ntCheck = 0 ;
		sum1    = 0.0 ;
		sum2    = 0.0 ;

		// start timer
		double t0 = MPI_Wtime() ;

		for (Myint it = 0; it < nt; it++)
		{
			Myfloat * prn = prnGrid->grid_3d ;
			Myfloat * prc = prcGrid->grid_3d ;

			// compute wavefield at next time step
			propa->computeWavefieldNextTimeStep(*prnGrid, *prcGrid) ;

			// check testCase results
			if (itry == 0)
			{
				if (it%snapInc == 0)
				{
					// initialize R at t = it * dt
					Myfloat timeSec = it * dt ;
					printInfo(MASTER, " Check at time", timeSec) ;
					propa->initializeGrid(*refGrid, EIGEN_MODE, timeSec) ;
					refGrid->write(caseName+"Ref") ;
					prnGrid->write(caseName+"Prn") ;

					ntCheck++ ;
					sum1 += prnGrid->getSumAbsDiff(INNER_POINTS, *refGrid) ;
					sum2 += refGrid->getSumAbs(INNER_POINTS) ;
				}
			}

			// swap grid
			auto tempGrid = prnGrid ;
			prnGrid = prcGrid ;
			prcGrid = tempGrid ;

			// update hardware counter (at regular time interval)
			hw->watchTimeAndUpdateHwCounter() ;

		} // for (Myint it = 0; it < nt; it++)

		// wait all process completed computations before ending time
		MPI_Barrier(MPI_COMM_WORLD) ;
		double t2 = MPI_Wtime() ;

		double testCase_time = t2-t0 ;
		printDebug(LIGHT_DEBUG, "testCase time", testCase_time) ;
		if (testCase_time < testCase_time_best)
		{
			testCase_time_best = testCase_time ;
		}

		// check error
		if (itry == 0)
		{

			printInfo(MASTER, " ntCheck", ntCheck) ;
			errTestCase = sum1 / sum2 ;
			printInfo(MASTER, " errTestCase2", errTestCase) ;
			checkBoolDiff((errTestCase<=maxErr), true) ;

		}

	} // for (Myint itry = 0; itry < ntry; itry++)

	// display perf
	printInfo(MASTER, " #Flop per point", nOpPerPoint) ;
	printInfo(MASTER, " #Point in stencil", nPtPerStencil) ;

	double timeInFD = testCase_time_best-testCase_time_com ;
	propaGflop     = nt * nGridPointGlob/timeInFD/1.e9 * nOpPerPoint;
	propaGpointEff = nt * nGridPointGlob/testCase_time_best/1.e9 ;
	propaGpointFD  = nt * nGridPointGlob/timeInFD/1.e9 ;
	propaGB        = nt * nGridPointGlob/timeInFD/1.e9 * nMemOpPerPoint * sizeof(Myfloat) ;

	printInfo(MASTER, " Best GFlop/s in FD" ,   propaGflop) ;
	printInfo(MASTER, " Best Gpoint/s eff." ,   propaGpointEff) ;
	printInfo(MASTER, " Best Gpoint/s in FD",   propaGpointFD) ;
	printInfo(MASTER, " Best Apparent BW GB/s", propaGB) ;

	if (nGridPointHaloGlob > 0)
	{
		printInfo(MASTER, " % MPI Comm\t", Myfloat(testCase_time_com / testCase_time_best * 100)) ;
		printInfo(MASTER, " MPI Comm. BW GB/s", Myfloat(nGridPointHaloGlob/testCase_time_com
				/ 1.e9 * sizeof(Myfloat))) ;
	}

	// log perf
	if (myMpiRank == 0)
	{
		perfLogFile
		// 10, 11, 12, 13
		<< propaGflop << " " << propaGpointFD << " " << propaGpointEff << " " << propaGB << " "

		// 14, 15
		<< testCase_time_best << " " << errTestCase << " "

		// 16, 17, 18
		<< propa->nt << " " << propa->dt << " " << propa->stableDt

		<< "\n" ;
	}

	this->finalize() ;

	printDebug(MID_DEBUG, "Out TestCase_Propa::run") ;
	return(RTN_CODE_OK) ;
}

} // namespace hpcscan

