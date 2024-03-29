
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
#include "table_Results.h"

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
	if (nMpiProc <= 1)
	{
		printWarning(" This test case requires at least 2 MPI procs") ;
		return(RTN_CODE_OK) ;
	}

	// for perf logging
	Table_Results *tabSend, *tabSendRecv ;
	Myfloat HaloExchGB = 0, HaloExchGPoint = 0, HaloExchBestTime = 0, HaloExchSize = 0 ;

	const string gridMode = Config::Instance()->testMode ;

	{
		//============================================
		// uni-directionnal MPI comm with MPI_Send
		// grid (proc x) -> grid (proc y)
		//============================================

		print_blank() ;
		string caseName = testCaseName + "MPI_Send" ;
		printInfo(MASTER, " * Case", caseName) ;
		tabSend = new Table_Results(caseName, nMpiProc, nMpiProc) ;

		// create the grid used in the MPI comm.
		Dim_type dim = DIM1 ;
		auto n1      = Config::Instance()->n1 ;
		auto n2      = Config::Instance()->n2 ;
		auto n3      = Config::Instance()->n3 ;
		auto fdOrder = Config::Instance()->fdOrder ;

		// the size of MPI message is set to the number of points in one halo in n1 direction

		// number of points in halo
		Myint64 nHaloGridPoint = (n2 * n3) *  Config::Instance()->fdOrder / 2 ;

		// minimum is set to 1000 points
		if (nHaloGridPoint < 1000) nHaloGridPoint = 1000 ;

		// halo size
		Myint64 haloSize = nHaloGridPoint * sizeof(Myfloat) ;

		// total number of values in MPI comm
		Myint64 nCommPoint = nHaloGridPoint ;

		// MPI comm size
		Myint64 commSize = nCommPoint * sizeof(Myfloat) ;
		printInfo(MASTER, " Total comm size (MB)", commSize / 1e6) ;

		Myint64 nPointForCreateGrid = nHaloGridPoint - Config::Instance()->fdOrder ;
		auto gridSrc2  = Grid_Factory::create(gridMode, GRID_GLOBAL, dim, nPointForCreateGrid, 0, 0) ;
		auto gridDest2 = Grid_Factory::create(gridMode, GRID_GLOBAL, dim, nPointForCreateGrid, 0, 0) ;
		Grid &gridSrc  = *gridSrc2 ;
		Grid &gridDest = *gridDest2 ;
		gridSrc.initializeGrid() ;
		gridDest.initializeGrid() ;

		Myint64 nGridPoint = gridSrc.getNumberOfGridPoint(GRID_GLOBAL, ALL_POINTS) ;
		if (nGridPoint != nHaloGridPoint)
		{
			printError("In TestCase_Comm::run, nGridPoint != nHaloGridPoint") ;
			printInfo(MASTER, "nGridPoint", nGridPoint) ;
			printInfo(MASTER, "nHaloGridPoint", nHaloGridPoint) ;
			return(RTN_CODE_KO) ;
		}

		Myint procSrcId, procDestId ;
		Myint ntry = Config::Instance()->ntry ;

		print_blank() ;
		// loop on the number of MPI procs senders
		for (Myint iprocSend = 0; iprocSend < nMpiProc; iprocSend++)
		{
			// loop on the number of MPI procs receivers
			for (Myint iprocRecv = 0; iprocRecv < nMpiProc; iprocRecv++)
			{

				Myfloat perfGB     = UNSPECIFIED ;
				Myfloat perfGPoint = UNSPECIFIED ;
				procSrcId  = iprocSend ;
				procDestId = iprocRecv ;

				// MPI_Send is blocking, src and dest can not be the same
				if (procSrcId != procDestId)
				{
					printInfo(MASTER, " * MPI_Send ===>", iprocSend, iprocRecv) ;

					// loop on the number of tries
					double testCase_time_best = FLT_MAX ;
					for (Myint itry = 0; itry < ntry; itry++)
					{
						gridSrc.fill(ALL_POINTS, iprocSend) ;
						gridDest.fill(ALL_POINTS, -1) ;

						// start MPI_Send & MPI_Recv
						MPI_Barrier(MPI_COMM_WORLD) ;
						double t0 = MPI_Wtime() ;

						// send grid with MPI_Send
						if (myMpiRank == procSrcId) gridSrc.sendWithMPI(nGridPoint, procDestId) ;

						// receive grid with MPI_Recv
						if (myMpiRank == procDestId) gridDest.recvWithMPI(nGridPoint, procSrcId) ;

						MPI_Barrier(MPI_COMM_WORLD) ;
						double t1 = MPI_Wtime() ;
						// end MPI comm

						double testCase_time = t1-t0 ;
						Myfloat testCase_bw = nCommPoint*sizeof(Myfloat)/testCase_time/1e9 ;
						printDebug(LIGHT_DEBUG, "Time", testCase_time) ;
						printDebug(LIGHT_DEBUG, "Speed", testCase_bw) ;
						testCase_time_best = min(testCase_time_best, testCase_time) ;

						// check grid
						if (itry == 0)
						{
							gridDest.write(caseName+to_string(procSrcId)+to_string(0)) ;

							Myint nbFailed = 0 ;
							if (myMpiRank == procDestId)
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

					if (myMpiRank == procDestId)
					{
						perfGPoint = nCommPoint/testCase_time_best/1.e9 ;
					}

				} // if ((procSrcId != MPI_PROC_NULL) && (procDestId != MPI_PROC_NULL))

				// get perf on master node
				Myfloat perfGPointglob ;
				MPI_Reduce(&perfGPoint, &perfGPointglob, 1, MPI_MYFLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
				if (perfGPointglob != UNSPECIFIED)
				{
					perfGB = perfGPointglob *sizeof(Myfloat) ;
					printInfo(MASTER, " Best achieved GByte/s", perfGB) ;
					printInfo(MASTER, " Best achieved GPoint/s", perfGPointglob) ;
					// store perf in table
					tabSend->seOneValue(procSrcId, procDestId, perfGB) ;
				}

			} // for (Myint iprocRecv = 0; iprocRecv < nMpiProc; iprocRecv++)
		} // for (Myint iprocSend = 0; iprocSend < nMpiProc; iprocSend++)
	}

	{
		//============================================
		// bi-directionnal MPI comm with MPI_Sendrecv
		// grid (proc x) <--> grid (proc y)
		//============================================

		print_blank() ;
		string caseName = testCaseName + "MPI_Sendrecv" ;
		printInfo(MASTER, " * Case", caseName) ;
		tabSendRecv = new Table_Results(caseName, nMpiProc, nMpiProc) ;

		// size of the MPI message corresponds to the number of point
		// in one halo in n1 direction (n2 * n3 * fdOrder / 2)
		// the total size is x2 (send and receive)
		// min is 1000 points

		// create the grid used in the MPI comm.
		Dim_type dim = DIM1 ;
		auto n1      = Config::Instance()->n1 ;
		auto n2      = Config::Instance()->n2 ;
		auto n3      = Config::Instance()->n3 ;
		auto fdOrder = Config::Instance()->fdOrder ;

		// the size of MPI message is set to the number of points in one halo in n1 direction

		// number of points in halo
		Myint64 nHaloGridPoint = (n2 * n3) *  Config::Instance()->fdOrder / 2 ;

		// minimum is set to 1000 points
		if (nHaloGridPoint < 1000) nHaloGridPoint = 1000 ;

		// halo size
		Myint64 haloSize = nHaloGridPoint * sizeof(Myfloat) ;

		// total number of values in MPI comm
		// x2 halo because send and receive
		Myint64 nCommPoint = nHaloGridPoint * 2 ;

		// MPI comm size
		Myint64 commSize = nCommPoint * sizeof(Myfloat) ;
		printInfo(MASTER, " Total comm size (MB)", commSize / 1e6) ;

		Myint64 nPointForCreateGrid = nHaloGridPoint - Config::Instance()->fdOrder ;
		auto gridSrc2  = Grid_Factory::create(gridMode, GRID_GLOBAL, dim, nPointForCreateGrid, 0, 0) ;
		auto gridDest2 = Grid_Factory::create(gridMode, GRID_GLOBAL, dim, nPointForCreateGrid, 0, 0) ;
		Grid &gridSrc  = *gridSrc2 ;
		Grid &gridDest = *gridDest2 ;
		gridSrc.initializeGrid() ;
		gridDest.initializeGrid() ;

		Myint64 nGridPoint = gridSrc.getNumberOfGridPoint(GRID_GLOBAL, ALL_POINTS) ;
		if (nGridPoint != nHaloGridPoint)
		{
			printError("In TestCase_Comm::run, nGridPoint != nHaloGridPoint") ;
			printInfo(MASTER, "nGridPoint", nGridPoint) ;
			printInfo(MASTER, "nHaloGridPoint", nHaloGridPoint) ;
			return(RTN_CODE_KO) ;
		}

		Myint procSrcId, procDestId ;
		Myint ntry = Config::Instance()->ntry ;

		print_blank() ;
		// loop on the number of MPI procs senders
		for (Myint iprocSend = 0; iprocSend < nMpiProc; iprocSend++)
		{
			// loop on the number of MPI procs receivers
			for (Myint iprocRecv = 0; iprocRecv < nMpiProc; iprocRecv++)
			{

				Myfloat perfGB     = UNSPECIFIED ;
				Myfloat perfGPoint = UNSPECIFIED ;
				procSrcId  = iprocSend ;
				procDestId = iprocRecv ;

				// MPI_Sendrecv is blocking, src and dest can not be the same
				if (procSrcId != procDestId)
				{
					printInfo(MASTER, " * MPI_Sendrecv <==>", procSrcId, procDestId) ;

					// loop on the number of tries
					double testCase_time_best = FLT_MAX ;
					for (Myint itry = 0; itry < ntry; itry++)
					{

						Myfloat *bufSend, *bufRecv ;
						Myint idSend, idRecv ;

						gridSrc.fill(ALL_POINTS, myMpiRank) ;
						gridDest.fill(ALL_POINTS, -1) ;
						bufSend = gridSrc.grid_3d ;
						bufRecv = gridDest.grid_3d ;

						if (myMpiRank == procSrcId)
						{
							idSend  = procDestId ;
							idRecv  = procDestId ;
						}
						else if (myMpiRank == procDestId)
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
						MPI_Barrier(MPI_COMM_WORLD) ;
						double t0 = MPI_Wtime() ;

						printDebug(LIGHT_DEBUG, " idSend ", idSend) ;
						printDebug(LIGHT_DEBUG, " idRecv ", idRecv) ;

						gridSrc.sendRecvWithMPI(gridDest, idSend, idRecv, nGridPoint) ;

						MPI_Barrier(MPI_COMM_WORLD) ;
						double t1 = MPI_Wtime() ;
						// end MPI comm

						double testCase_time = t1-t0 ;
						Myfloat testCase_bw = nCommPoint*sizeof(Myfloat)/testCase_time/1e9 ;
						printDebug(LIGHT_DEBUG, "Time", testCase_time) ;
						printDebug(LIGHT_DEBUG, "Speed", testCase_bw) ;
						testCase_time_best = min(testCase_time_best, testCase_time) ;

						// check grid
						if (itry == 0)
						{
							if ((myMpiRank == procSrcId) || (myMpiRank == procDestId))
							{
								gridDest.write(caseName+to_string(procSrcId)+to_string(procDestId)) ;
							}

							Myint nbFailed = 0 ;
							if (myMpiRank == procSrcId)
							{
								Myfloat minVal = gridDest.getMin(ALL_POINTS) ;
								if (relErr(minVal, procDestId) > MAX_ERR_FLOAT) nbFailed++ ;
								Myfloat maxVal = gridDest.getMax(ALL_POINTS) ;
								if (relErr(maxVal, procDestId) > MAX_ERR_FLOAT) nbFailed++ ;
							}

							else if (myMpiRank == procDestId)
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

					if (myMpiRank == procDestId)
					{
						perfGPoint = nCommPoint/testCase_time_best/1.e9 ;
					}

				} // if (procSrcId != procDestId)

				// get perf on master node
				Myfloat perfGPointglob ;
				MPI_Reduce(&perfGPoint, &perfGPointglob, 1, MPI_MYFLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
				if (perfGPointglob != UNSPECIFIED)
				{
					perfGB = perfGPointglob * sizeof(Myfloat) ;
					printInfo(MASTER, " Best achieved GByte/s", perfGB) ;
					printInfo(MASTER, " Best achieved GPoint/s", perfGPointglob) ;
					// store perf in table
					tabSendRecv->seOneValue(procSrcId, procDestId, perfGB) ;
				}

			} // for (Myint iprocRecv = 0; iprocRecv < nMpiProc; iprocRecv++)
		} // for (Myint iprocSend = 0; iprocSend < nMpiProc; iprocSend++)
	}

	{
		//============================================
		// Exchange Halos with Sendrecv
		// All procs to all procs
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
			MPI_Barrier(MPI_COMM_WORLD) ;
			double t0 = MPI_Wtime() ;
			Ugrid.exchangeHalos(MPI_COMM_MODE_SENDRECV) ;
			MPI_Barrier(MPI_COMM_WORLD) ;
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

		HaloExchBestTime = testCase_time_best ;
		HaloExchSize = nGridPointHaloGlob * sizeof(Myfloat) ;
		HaloExchGPoint = nGridPointHaloGlob / HaloExchBestTime / 1.e9 ;
		HaloExchGB = HaloExchGPoint * sizeof(Myfloat) ;

		printInfo(MASTER, " Best achieved GByte/s", HaloExchGB) ;
		printInfo(MASTER, " Best achieved GPoint/s", HaloExchGPoint) ;
	}

	// display performance summary
	print_blank() ;
	print_line2() ;
	printInfo(MASTER, " Performance summary GByte/s") ;
	tabSend->display() ;
	tabSendRecv->display() ;
	print_blank() ;
	printInfo(MASTER, " *** Exchange halos", HaloExchGB) ;

	delete(tabSend) ;
	delete(tabSendRecv) ;

	// log perf
	if (myMpiRank == 0)
	{
		// first number is at position 10 in log file for numeric values
		for (Myint ii = 0; ii < tabSend->nVal ; ii++)
		{
			perfLogFile << tabSend->val[ii] << " " ;
		}
		for (Myint ii = 0; ii < tabSendRecv->nVal ; ii++)
		{
			perfLogFile << tabSendRecv->val[ii] << " " ;
		}
		perfLogFile << HaloExchSize/1e6 << " " << HaloExchBestTime << " " << HaloExchGB << "\n" ;
	}

	this->finalize() ;

	printDebug(MID_DEBUG, "Out TestCase_Comm::run") ;
	return(RTN_CODE_OK) ;
}

} // namespace hpcscan

