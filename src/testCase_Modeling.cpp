
//-------------------------------------------------------------------------------------------------------
// test case to perform modeling
//-------------------------------------------------------------------------------------------------------

#include "testCase_Modeling.h"

#include <algorithm>
#include <cfloat>

#include <mpi.h>
#include <omp.h>

#include "config.h"
#include "data_Acquisition_Factory.h"
#include "fdm.h"
#include "global.h"
#include "grid_Factory.h"
#include "output_report.h"
#include "propagator_Factory.h"

namespace hpcscan
{

	TestCase_Modeling::TestCase_Modeling(void)
	{
		testCaseName = "Modeling";
		testCaseVersion = "Standard implementation";
	}

	//-------------------------------------------------------------------------------------------------------

	Rtn_code TestCase_Modeling::run(void)
	{

		printDebug(MID_DEBUG, "In TestCase_Modeling::run");

		if (this->initialize() == RTN_CODE_KO)
			return (RTN_CODE_KO);

		// instantiate propagator
		const string propagator = Config::Instance()->propagator;
		auto propa = Propagator_Factory::create(propagator);
		if (propa == nullptr)
			return (RTN_CODE_KO);

		// initialize propagator
		Rtn_code rtnCode = propa->initialize(MODELING);
		if (rtnCode != RTN_CODE_OK)
		{
			printError("In TestCase_Modeling::run, propa->initialize() not Ok");
			return (RTN_CODE_OK);
		}

		// print propagator info
		propa->info();

		// get pointers on propagators grids
		auto prnGrid = propa->prnGrid;
		auto prcGrid = propa->prcGrid;

		if (myMpiRank == 0)
			prnGrid->info();

		// Number of points exchanged in Halos
		Myint64 nGridPointHaloLoc = prnGrid->getNumberOfGridPointCommHalo(GRID_LOCAL);
		Myint64 nGridPointHaloGlob = prnGrid->getNumberOfGridPointCommHalo(GRID_GLOBAL);
		print_blank();
		printInfo(MASTER, " Halos Com Local (Pts)", nGridPointHaloLoc);
		printInfo(MASTER, " Halos Com Global (Pts)", nGridPointHaloGlob);
		printInfo(MASTER, " Halos Com Global (MB)", Myfloat(nGridPointHaloGlob * sizeof(Myfloat) / 1.e6));

		double testCase_time_best = 0.0, testCase_time_com = 0.0;

		Myint ntry = Config::Instance()->ntry;

		Myint64 nGridPointLoc = prnGrid->getNumberOfGridPoint(GRID_LOCAL, INNER_POINTS);
		Myint64 nGridPointGlob = prnGrid->getNumberOfGridPoint(GRID_GLOBAL, INNER_POINTS);
		printInfo(MASTER, " Grid Inner Loc. (Pts)", nGridPointLoc);
		printInfo(MASTER, " Grid Inner Glob. (Pts)", nGridPointGlob);
		printInfo(MASTER, " Grid Inner Glob. (MB)", Myfloat(nGridPointGlob * sizeof(Myfloat) / 1.e6));
		print_blank();

		// get parameters from propagator
		const Myint fdOrder = propa->fdOrder;
		const Myint nt = propa->nt;
		const Myfloat dt = propa->dt;
		const Myint snapInc = propa->snapInc;

		const Myfloat maxErr = 1.0;
		printInfo(MASTER, " Max allowed error", maxErr);

		Myint nPtPerStencil = prcGrid->getPtPerStencilFD_LAPLACIAN(fdOrder);
		Myint nMemOpPerPoint = nPtPerStencil + 4;							// + 1 store + 3 load
		Myint nOpPerPoint = prcGrid->getFlopPerPtFD_LAPLACIAN(fdOrder) + 4; // + 1 ADD + 1 SUB + 2 MUL

		// for perf log
		Myfloat propaGflop = 0, propaGpointEff = 0, propaGpointFD = 0, propaGB = 0;

		// initialize data acquisition		
		auto acqui2 = DataAcquisition_Factory::create(Config::Instance()->testMode) ;		
		DataAcquisition &acqui = *acqui2 ;
		acqui.initialize(ACQUI_BUILTIN, *prnGrid, nt) ;

		print_blank();
		string caseName = testCaseName + "Forward";
		printInfo(MASTER, " * Case", caseName);

		testCase_time_best = FLT_MAX;

		Myfloat errTestCase;
		Myint ntCheck;
		Myint n_iter_snap=0;
		double total_compute_time = 0;

		for (Myint itry = 0; itry < ntry; itry++)
		{
			ntCheck = 0;

			// start timer
			double t0 = MPI_Wtime();
			total_compute_time = 0;

			for (Myint it = 0; it < nt; it++)
			{

				// compute wavefield at next time step
				double t0_compute = MPI_Wtime();
				
				propa->computeWavefieldNextTimeStep(*prnGrid, *prcGrid);

				double t1_compute = MPI_Wtime();

				total_compute_time += t1_compute - t0_compute;

				// apply source term
				acqui.appliSourceTerm(*prnGrid, FORWARD, it, dt) ;

				// output snapshots
				if (it % snapInc == 0)
				{
					Myfloat timeSec = it * dt;
					printInfo(MASTER, " Snapshot at time", timeSec);
					//prnGrid->write(caseName + "Prn");
					prnGrid->writeGlobal(INNER_POINTS,caseName + "Prn",n_iter_snap);
					n_iter_snap++;
				}	

				// record traces
				acqui.recordTrace(*prnGrid, FORWARD, it, dt) ;			

				// swap grid
				auto tempGrid = prnGrid;
				prnGrid = prcGrid;
				prcGrid = tempGrid;

				// update hardware counter (at regular time interval)
				hw->watchTimeAndUpdateHwCounter();

			} // for (Myint it = 0; it < nt; it++)

			// write traces
			acqui.writeTrace(caseName + "Prn") ;

			// wait all process completed computations before ending time
			MPI_Barrier(MPI_COMM_WORLD);
			double t2 = MPI_Wtime();

			double testCase_time = t2 - t0;
			printDebug(LIGHT_DEBUG, "testCase time", testCase_time);
			if (testCase_time < testCase_time_best)
			{
				testCase_time_best = testCase_time;
			}

			// check error
			if (itry == 0)
			{

				printInfo(MASTER, " ntCheck", ntCheck);
				// errTestCase = sum1 / sum2;
				// printInfo(MASTER, " errTestCase2", errTestCase);
				// checkBoolDiff((errTestCase <= maxErr), true);
			}

		} // for (Myint itry = 0; itry < ntry; itry++)

		// display perf
		printInfo(MASTER, " #Flop per point", nOpPerPoint);
		printInfo(MASTER, " #Point in stencil", nPtPerStencil);

		// double timeInFD = testCase_time_best - testCase_time_com;
		double timeInFD = propa->compute_pressure_time;
		propaGflop		= nt * nGridPointGlob / timeInFD / 1.e9 * nOpPerPoint;
		propaGpointEff	= nt * nGridPointGlob / testCase_time_best / 1.e9;
		propaGpointFD	= nt * nGridPointGlob / timeInFD / 1.e9;
		propaGB			= nt * nGridPointGlob / timeInFD / 1.e9 * nMemOpPerPoint * sizeof(Myfloat);

		printInfo(MASTER, " Total compute time ", total_compute_time);
		printInfo(MASTER, " Total compute speed ", nt * nGridPointGlob / total_compute_time / 1.e9);
		
		// i/o is = write grid (if asked) + writetrace (after gather) + read model (plane_by_plane)
		double i_o_time = propa->i_o_read_time + acqui.trace_write_time 
						+ prnGrid->grid_writeGlobal_time + propa->i_o_write_coefgrid_time;
		printInfo(MASTER, "\n Total i/o time ", i_o_time);
		
		double MPI_COMM_time = propa->halo_comm_time + acqui.trace_gather_time;
		// communication is = exchange_halo (sendrecv) + writeTrace (gather)
		printInfo(MASTER, " Total communication time ", MPI_COMM_time);

		// summarizing each steps
		printInfo(MASTER, "\n total Read model time ", propa->coef_read_total_time);
		printInfo(MASTER, " i_o Read model time ", propa->i_o_read_time);
		printInfo(MASTER, " i_o Write trace time ", acqui.trace_write_time);

		printInfo(MASTER, " Halo communication time ", propa->halo_comm_time);
		printInfo(MASTER, " Apply boundary condition time ", propa->boundary_condition_time);
		printInfo(MASTER, " Compute pressure time ", propa->compute_pressure_time);
		printInfo(MASTER, " Record trace time ", acqui.record_trace_time);
		if(Config::Instance()->writeGrid)
		{
			printInfo(MASTER, " Total write snapshotgrid time ", prnGrid->grid_writeGlobal_time);
			printInfo(MASTER, " Total write coefGrid time ", propa->i_o_write_coefgrid_time);
		}
		printInfo(MASTER, " Total write trace time ", acqui.total_timeWriteTrace);

		printInfo(MASTER, " \nBest GFlop/s in FD", propaGflop);
		printInfo(MASTER, " Best Gpoint/s eff.", propaGpointEff);
		printInfo(MASTER, " Best Gpoint/s in FD", propaGpointFD);
		printInfo(MASTER, " Best Apparent BW GB/s", propaGB);

		if (nGridPointHaloGlob > 0)
		{
			printInfo(MASTER, " % MPI Comm\t", Myfloat(testCase_time_com / testCase_time_best * 100));
			printInfo(MASTER, " MPI Comm. BW GB/s", Myfloat(nGridPointHaloGlob / testCase_time_com / 1.e9 * sizeof(Myfloat)));
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
			<< propa->nt << " " << propa->dt << " " << propa->stableDt << " "

			// 19, 20, 21 (main loop)
			<< propa->halo_comm_time << " " << propa->boundary_condition_time << " " << propa->compute_pressure_time << " "

			// 22, 23 (compute total time, + speed)
			<< total_compute_time << " " << nt * nGridPointGlob / total_compute_time / 1.e9 << " "

			// 24, 25 (2 main i/o)
			<< propa->i_o_read_time << " " << acqui.trace_write_time << " " 

			// 26, 27, 28 (total time spend for each part)
			<< propa->coef_read_total_time << " " << acqui.record_trace_time << " " << acqui.total_timeWriteTrace << " "

			// 29 (trace gather communication time)
			<< acqui.trace_gather_time << " "

			// 30, 31 (total i/o & comm)
			<< i_o_time << " " << MPI_COMM_time << " ";

			// 32, 33
			if(Config::Instance()->writeGrid)
			{
				perfLogFile << prnGrid->grid_writeGlobal_time << " " << propa->i_o_write_coefgrid_time << " ";
			}

			perfLogFile << "\n";
		}

		this->finalize();

		printDebug(MID_DEBUG, "Out TestCase_Modeling::run");
		return (RTN_CODE_OK);
	}

} // namespace hpcscan
