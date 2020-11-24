
//-------------------------------------------------------------------------------------------------------
// test cases to measure the MPI communication bandwidth
// baseline kernel - no optimization
//-------------------------------------------------------------------------------------------------------

#include "testCase_Comm.h"

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

TestCase_Comm::TestCase_Comm(void)
		{
	testCaseName    = "Comm" ;
	testCaseVersion = "Standard implementation" ;
		}

//-------------------------------------------------------------------------------------------------------

Rtn_code TestCase_Comm::run(void)
{

	printDebug(MID_DEBUG, "In TestCase_Comm::run") ;

	if (this->initialize() == RTN_CODE_KO) return (RTN_CODE_KO) ;

	// test case requires at least 2 MPI proc
	if (nproc_world <= 1)
	{
		printWarning(" This test case requires at least 2 MPI procs") ;
		return(RTN_CODE_OK) ;
	}

	// for perf log file
	Myfloat SendGB=0     , SendRecvGB=0     , HaloExchGB=0 ;
	Myfloat SendGPoint=0 , SendRecvGPoint=0 , HaloExchGPoint=0 ;

	const string gridMode = Config::Instance()->testMode ;

	{
		//============================================
		// uni-directionnal MPI comm with MPI_Send
		// grid (proc x) -> grid (proc 0)
		//============================================

		// *** NOTE ***
		// For this test case, we overwrite config parameters
		// to keep grid size constant whatever the config is

		Dim_type dim = DIM3 ;
		Myint64 n1 = 1000 ;
		Myint64 n2 = 100 ;
		Myint64 n3 = 100 ;

		auto gridSrc2  = Grid_Factory::create(gridMode, GRID_GLOBAL, dim, n1, n2, n3) ;
		auto gridDest2 = Grid_Factory::create(gridMode, GRID_GLOBAL, dim, n1, n2, n3) ;
		Grid &gridSrc  = *gridSrc2 ;
		Grid &gridDest = *gridDest2 ;
		gridSrc.initializeGrid() ;
		gridDest.initializeGrid() ;

		Myint64 nGridPoint = gridSrc.getNumberOfGridPoint(GRID_GLOBAL, ALL_POINTS) ;

		// display grid info
		gridSrc.info() ;

		print_blank() ;
		string caseName = testCaseName + "MPI_Send" ;
		printInfo(MASTER, " * Case", caseName) ;

		Myint procSrcId, procDestId ;
		Myint ntry = Config::Instance()->ntry ;

		// loop on the number of MPI procs
		for (Myint iproc = 1; iproc < nproc_world; iproc++)
		{
			// sender is iproc and receiver is 0
			// other procs are not involved in communication
			if (myid_world == 0)
			{
				procSrcId = iproc ;
			}
			else
			{
				procSrcId = MPI_PROC_NULL ;
			}

			if (myid_world == iproc)
			{
				procDestId = 0 ;
			}
			else
			{
				procDestId = MPI_PROC_NULL ;
			}

			// loop on the number of tries
			double testCase_time_best = FLT_MAX ;
			for (Myint itry = 0; itry < ntry; itry++)
			{
				gridSrc.fill(ALL_POINTS, iproc) ;
				gridDest.fill(ALL_POINTS, -1) ;

				// start MPI_Send & MPI_Recv
				double t0 = MPI_Wtime() ;

				MPI_Status status ;

				// MPI_Recv
				MPI_Recv(gridDest.grid_3d, nGridPoint, MPI_MYFLOAT, procSrcId, 0, MPI_COMM_WORLD, &status);
				if (status.MPI_ERROR != MPI_SUCCESS)
				{
					//printError("MPI ERROR", status.MPI_ERROR) ;
				}

				// MPI_Send
				MPI_Send(gridSrc.grid_3d, nGridPoint, MPI_MYFLOAT, procDestId, 0, MPI_COMM_WORLD);

				double t1 = MPI_Wtime() ;
				// end MPI comm

				double testCase_time = t1-t0 ;
				Myfloat testCase_bw = nGridPoint*sizeof(Myfloat)/testCase_time/1e9 ;
				printDebug(LIGHT_DEBUG, "Time", testCase_time) ;
				printDebug(LIGHT_DEBUG, "Speed", testCase_bw) ;
				testCase_time_best = min(testCase_time_best, testCase_time) ;

				// check grid
				if (itry == 0)
				{
					printInfo(MASTER, " * MPI_Send ===>", procSrcId, 0) ;
					if (myid_world == 0)
					{
						gridDest.write(caseName+to_string(procSrcId)+to_string(0)) ;
					}

					Myint nbFailed = 0 ;
					if (myid_world == 0)
					{
						Myfloat minVal = gridDest.getMin(ALL_POINTS) ;
						if (relErr(minVal, iproc) > MAX_ERR_FLOAT) nbFailed++ ;
					}

					if (myid_world == 0)
					{
						Myfloat maxVal = gridDest.getMax(ALL_POINTS) ;
						if (relErr(maxVal, iproc) > MAX_ERR_FLOAT) nbFailed++ ;
					}

					Myint nbFailedTot ;
					MPI_Reduce(&nbFailed, &nbFailedTot, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
					checkIntegerDiff(nbFailedTot, 0) ;
				}
			} // for (Myint itry = 0; itry < ntry; itry++)

			SendGPoint = nGridPoint/testCase_time_best/1.e9 ;
			SendGB = SendGPoint *sizeof(Myfloat) ;
			printInfo(MASTER, " Best achieved GByte/s", SendGB) ;
			printInfo(MASTER, " Best achieved GPoint/s", SendGPoint) ;

			MPI_Barrier(MPI_COMM_WORLD) ;

		} // for (Myint iproc = 0; iproc < nproc_world; iproc++)
	}

	{
		//============================================
		// bi-directionnal MPI comm with MPI_Sendrecv
		// grid (proc x) <--> grid (proc 0)
		//============================================

		// *** NOTE ***
		// For this test case, we overwrite config parameters
		// to keep grid size constant whatever the config is

		Dim_type dim = DIM3 ;
		Myint64 n1 = 1000 ;
		Myint64 n2 = 100 ;
		Myint64 n3 = 100 ;

		auto gridSrc2  = Grid_Factory::create(gridMode, GRID_GLOBAL, dim, n1, n2, n3) ;
		auto gridDest2 = Grid_Factory::create(gridMode, GRID_GLOBAL, dim, n1, n2, n3) ;
		Grid &gridSrc  = *gridSrc2 ;
		Grid &gridDest = *gridDest2 ;
		gridSrc.initializeGrid() ;
		gridDest.initializeGrid() ;

		Myint64 nGridPoint = gridSrc.getNumberOfGridPoint(GRID_GLOBAL, ALL_POINTS) ;

		// display grid info
		gridSrc.info() ;

		print_blank() ;
		string caseName = testCaseName + "MPI_Sendrecv" ;
		printInfo(MASTER, " * Case", caseName) ;

		Myint ntry = Config::Instance()->ntry ;
		Myint procDestId = 0 ;

		// loop on the number of MPI procs
		for (Myint iproc = 1; iproc < nproc_world; iproc++)
		{
			Myint procSrcId = iproc ;

			// loop on the number of tries
			double testCase_time_best = FLT_MAX ;
			for (Myint itry = 0; itry < ntry; itry++)
			{

				Myfloat *bufSend, *bufRecv ;
				Myint idSend, idRecv ;

				gridSrc.fill(ALL_POINTS, myid_world) ;
				gridDest.fill(ALL_POINTS, -1) ;
				bufSend = gridSrc.grid_3d ;
				bufRecv = gridDest.grid_3d ;

				if (myid_world == procSrcId)
				{
					idSend  = procDestId ;
					idRecv  = procDestId ;
				}
				else if (myid_world == procDestId)
				{
					idSend  = procSrcId ;
					idRecv  = procSrcId ;
				}
				else
				{
					idSend  = MPI_PROC_NULL ;
					idRecv  = MPI_PROC_NULL ;
				}

				// start MPI_Sendrecv
				double t0 = MPI_Wtime() ;
				MPI_Status status ;
				printDebug(LIGHT_DEBUG, " idSend ", idSend) ;
				printDebug(LIGHT_DEBUG, " idRecv ", idRecv) ;
				MPI_Sendrecv(bufSend, nGridPoint, MPI_MYFLOAT, idSend, 0,
						bufRecv, nGridPoint, MPI_MYFLOAT, idRecv, 0,
						MPI_COMM_WORLD, &status);

				double t1 = MPI_Wtime() ;
				if (status.MPI_ERROR != MPI_SUCCESS)
				{
					//printError("MPI ERROR", status.MPI_ERROR) ;
				}
				// end MPI comm

				double testCase_time = t1-t0 ;
				Myfloat testCase_bw = 2*nGridPoint*sizeof(Myfloat)/testCase_time/1e9 ;
				printDebug(LIGHT_DEBUG, "Time", testCase_time) ;
				printDebug(LIGHT_DEBUG, "Speed", testCase_bw) ;
				testCase_time_best = min(testCase_time_best, testCase_time) ;

				// check grid
				if (itry == 0)
				{
					printInfo(MASTER, " * MPI_Sendrecv <==>", procSrcId, procDestId) ;

					if ((myid_world == procSrcId) || (myid_world == procDestId))
					{
						gridDest.write(caseName+to_string(procSrcId)+to_string(procDestId)) ;
					}

					Myint nbFailed = 0 ;
					if (myid_world == procSrcId)
					{
						Myfloat minVal = gridDest.getMin(ALL_POINTS) ;
						if (relErr(minVal, procDestId) > MAX_ERR_FLOAT) nbFailed++ ;
						Myfloat maxVal = gridDest.getMax(ALL_POINTS) ;
						if (relErr(maxVal, procDestId) > MAX_ERR_FLOAT) nbFailed++ ;
					}

					else if (myid_world == procDestId)
					{
						Myfloat minVal = gridDest.getMin(ALL_POINTS) ;
						if (relErr(minVal, procSrcId) > MAX_ERR_FLOAT) nbFailed++ ;
						Myfloat maxVal = gridDest.getMax(ALL_POINTS) ;
						if (relErr(maxVal, procSrcId) > MAX_ERR_FLOAT) nbFailed++ ;
					}

					Myint nbFailedTot ;
					MPI_Reduce(&nbFailed, &nbFailedTot, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
					checkIntegerDiff(nbFailedTot, 0) ;
				}

			} // for (Myint itry = 0; itry < ntry; itry++)

			SendRecvGPoint = 2*nGridPoint/testCase_time_best/1.e9 ;
			SendRecvGB = SendRecvGPoint * sizeof(Myfloat) ;

			printInfo(MASTER, " Best achieved GByte/s", SendRecvGB) ;
			printInfo(MASTER, " Best achieved GPoint/s", SendRecvGPoint) ;

		} // for (Myint iproc = 0; iproc < nproc_world; iproc++)
	}

	{
		//============================================
		// Exchange Halos
		//============================================

		print_blank() ;
		string caseName = testCaseName + "ExchangeHalos" ;
		printInfo(MASTER, " * Case", caseName) ;

		auto Ugrid2  = Grid_Factory::create(gridMode, GRID_LOCAL) ;
		auto Rgrid2 = Grid_Factory::create(gridMode, GRID_LOCAL) ;
		Grid &Ugrid = *Ugrid2 ;
		Grid &Rgrid = *Rgrid2 ;
		Ugrid.initializeGrid() ;
		Rgrid.initializeGrid() ;

		Ugrid.info() ;

		const Myfloat64 a1 = Config::Instance()->param1 ;
		const Myfloat64 a2 = Config::Instance()->param2 ;
		const Myfloat64 a3 = Config::Instance()->param3 ;
		const Myfloat64 a4 = Config::Instance()->param4 ;

		Ugrid.fill(ALL_POINTS, FUNC_LINEAR, FUNC_LINEAR, FUNC_LINEAR, a1, a2, a3, a4) ;
		Rgrid.fill(ALL_POINTS, FUNC_LINEAR, FUNC_LINEAR, FUNC_LINEAR, a1, a2, a3, a4) ;

		// Number of inner points
		Myint64 nGridPointInnerLoc  = Ugrid.getNumberOfGridPoint(GRID_LOCAL, INNER_POINTS) ;
		Myint64 nGridPointInnerGlob = Ugrid.getNumberOfGridPoint(GRID_GLOBAL, INNER_POINTS) ;
		printInfo(MASTER, " Inner pts (Loc.)", nGridPointInnerLoc) ;
		printInfo(MASTER, " Inner pts (Glob.)", nGridPointInnerGlob) ;

		// Number of points exchanged in Halos
		Myint64 nGridPointHaloLoc = Ugrid.getNumberOfGridPointCommHalo(GRID_LOCAL) ;
		Myint64 nGridPointHaloGlob = Ugrid.getNumberOfGridPointCommHalo(GRID_GLOBAL) ;
		print_blank() ;
		printInfo(MASTER, " Halos Local (Points)", nGridPointHaloLoc) ;
		printInfo(MASTER, " Halos Global (Points)", nGridPointHaloGlob) ;
		printInfo(MASTER, " Halos Global (MB)", Myfloat(nGridPointHaloGlob*sizeof(Myfloat)/1.e6)) ;

		printInfo(MASTER, " % Halo vs Inner", (Myfloat) nGridPointHaloGlob/nGridPointInnerGlob*100) ;

		double testCase_time_best = FLT_MAX ;

		Myint ntry = Config::Instance()->ntry ;

		for (Myint itry = 0; itry < ntry; itry++)
		{
			double t0 = MPI_Wtime() ;
			Ugrid.exchangeHalos(MPI_COMM_MODE_SENDRECV) ;
			double t1 = MPI_Wtime() ;

			double testCase_time = t1-t0 ;
			printDebug(LIGHT_DEBUG, "Time", testCase_time) ;
			testCase_time_best = min(testCase_time_best, testCase_time) ;

			// check testCase results
			if (itry == 0)
			{
				checkAllProcGridL1Err(ALL_POINTS, Ugrid, Rgrid, MAX_ERR_FLOAT) ;
			}

		}

		HaloExchGPoint = nGridPointHaloGlob/testCase_time_best/1.e9 ;
		HaloExchGB = HaloExchGPoint * sizeof(Myfloat) ;

		printInfo(MASTER, " Best achieved GByte/s", HaloExchGB) ;
		printInfo(MASTER, " Best achieved GPoint/s", HaloExchGPoint) ;
	}

	// log perf
	if (myid_world == 0)
	{
		perfLogFile
		<< SendGB << " " << SendGPoint << " "
		<< SendRecvGB << " " << SendRecvGPoint << " "
		<< HaloExchGB << " " << HaloExchGPoint << " "
		<< "\n" ;
	}

	this->finalize() ;

	printDebug(MID_DEBUG, "Out TestCase_Comm::run") ;
	return(RTN_CODE_OK) ;
}

} // namespace hpcscan

