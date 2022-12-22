
//-------------------------------------------------------------------------------------------------------
// test cases to measure the memory bandwidth
// baseline kernel - no optimization
//-------------------------------------------------------------------------------------------------------

#include "testCase_Memory.h"

#include <algorithm> // for min()
#include <cfloat>	 // for FLT_MAX

#include <mpi.h>
#include <omp.h>

#include "config.h"
#include "constant.h"
#include "global.h"
#include "grid_Factory.h"
#include "output_report.h"

namespace hpcscan
{

	TestCase_Memory::TestCase_Memory(void)
	{
		testCaseName = "Memory";
		testCaseVersion = "Standard implementation";
	}

	//-------------------------------------------------------------------------------------------------------

	Rtn_code TestCase_Memory::run(void)
	{

		printDebug(MID_DEBUG, "In TestCase_Memory::run");

		if (this->initialize() == RTN_CODE_KO)
			return (RTN_CODE_KO);

		const Myfloat a1 = Config::Instance()->param1;
		const Myfloat a2 = Config::Instance()->param2;

		const string gridMode = Config::Instance()->testMode;
		auto Ugrid2 = Grid_Factory::create(gridMode, GRID_GLOBAL);
		auto Vgrid2 = Grid_Factory::create(gridMode, GRID_GLOBAL);
		auto Wgrid2 = Grid_Factory::create(gridMode, GRID_GLOBAL);
		auto Rgrid2 = Grid_Factory::create(gridMode, GRID_GLOBAL);
		Grid &Ugrid = *Ugrid2;
		Grid &Vgrid = *Vgrid2;
		Grid &Wgrid = *Wgrid2;
		Grid &Rgrid = *Rgrid2; // reference grid
		Ugrid.initializeGrid();
		Vgrid.initializeGrid();
		Wgrid.initializeGrid();
		Rgrid.initializeGrid();

		if (myMpiRank == 0)
			Ugrid.info();

		const Myint64 nGridPoint = Ugrid.getNumberOfGridPoint(GRID_GLOBAL, ALL_POINTS);

		// for perf log
		Myfloat FillGridGB = 0, FillGridGpoint = 0, CopyGridGB = 0, CopyGridGpoint = 0;
		Myfloat AddGridGB = 0, AddGridGpoint = 0, MultiplyGridGB = 0, MultiplyGridGpoint = 0;
		Myfloat AddUpdateGridGB = 0, AddUpdateGridGpoint = 0;

		// one block per case
		{
			//============================================
			// fill array
			// W = coef
			//============================================

			print_blank();
			string caseName = testCaseName + "FillArray";
			printInfo(MASTER, " * Case", caseName);

			Wgrid.fill(ALL_POINTS, 0.0);

			double testCase_time_best = FLT_MAX;

			Myint ntry = Config::Instance()->ntry;

			for (Myint itry = 0; itry < ntry; itry++)
			{
				double t0 = MPI_Wtime();
				Wgrid.fillArray(a1);
				double t1 = MPI_Wtime();

				double testCase_time = t1 - t0;
				Myfloat testCase_bw = nGridPoint * sizeof(Myfloat) / testCase_time / 1e9;
				printDebug(LIGHT_DEBUG, "Time", testCase_time);
				printDebug(LIGHT_DEBUG, "Speed", testCase_bw);
				testCase_time_best = min(testCase_time_best, testCase_time);

				// check testCase results
				if (itry == 0)
				{
					Myfloat minVal = Wgrid.getMin(ALL_POINTS);
					checkFloatDiff(minVal, a1, MAX_ERR_FLOAT);
					Myfloat maxVal = Wgrid.getMax(ALL_POINTS);
					checkFloatDiff(maxVal, a1, MAX_ERR_FLOAT);
				}
			}

			FillGridGB = Myfloat(nGridPoint * sizeof(Myfloat) / testCase_time_best / 1.e9);
			FillGridGpoint = Myfloat(nGridPoint / testCase_time_best / 1.e9);
			printInfo(MASTER, " Best achieved GByte/s", FillGridGB);
			printInfo(MASTER, " Best achieved GPoint/s", FillGridGpoint);
		}

		{
			//============================================
			// copy array
			// W = U
			//============================================

			print_blank();
			string caseName = testCaseName + "CopyArray";
			printInfo(MASTER, " * Case", caseName);

			Ugrid.fill(ALL_POINTS, a1);
			Wgrid.fill(ALL_POINTS, 0.0);

			double testCase_time_best = FLT_MAX;

			Myint ntry = Config::Instance()->ntry;

			for (Myint itry = 0; itry < ntry; itry++)
			{
				double t0 = MPI_Wtime();
				Wgrid.copyArray(Ugrid);
				double t1 = MPI_Wtime();

				double testCase_time = t1 - t0;
				Myfloat testCase_bw = nGridPoint * 2 * sizeof(Myfloat) / testCase_time / 1e9;
				printDebug(LIGHT_DEBUG, "Time", testCase_time);
				printDebug(LIGHT_DEBUG, "Speed", testCase_bw);
				testCase_time_best = min(testCase_time_best, testCase_time);

				// check testCase results
				if (itry == 0)
				{
					checkGridL1Err(ALL_POINTS, Wgrid, Ugrid, MAX_ERR_FLOAT);
					checkGridMaxErr(ALL_POINTS, Wgrid, Ugrid, MAX_ERR_FLOAT);
				}
			}

			CopyGridGB = Myfloat(nGridPoint * 2 * sizeof(Myfloat) / testCase_time_best / 1.e9);
			CopyGridGpoint = Myfloat(nGridPoint / testCase_time_best / 1.e9);
			printInfo(MASTER, " Best achieved GByte/s", CopyGridGB);
			printInfo(MASTER, " Best achieved GPoint/s", CopyGridGpoint);
		}

		{
			//============================================
			// add array
			// W = U + V
			//============================================

			print_blank();
			string caseName = testCaseName + "AddArray";
			printInfo(MASTER, " * Case", caseName);

			Ugrid.fill(ALL_POINTS, a1);
			Vgrid.fill(ALL_POINTS, a2);
			Wgrid.fill(ALL_POINTS, 0.0);
			Rgrid.fill(ALL_POINTS, a1 + a2);

			double testCase_time_best = FLT_MAX;

			Myint ntry = Config::Instance()->ntry;

			for (Myint itry = 0; itry < ntry; itry++)
			{
				double t0 = MPI_Wtime();
				Wgrid.addArray(Ugrid, Vgrid);
				double t1 = MPI_Wtime();

				double testCase_time = t1 - t0;
				Myfloat testCase_bw = nGridPoint * 3 * sizeof(Myfloat) / testCase_time / 1e9;
				printDebug(LIGHT_DEBUG, "Time", testCase_time);
				printDebug(LIGHT_DEBUG, "Speed", testCase_bw);
				testCase_time_best = min(testCase_time_best, testCase_time);

				// check testCase results
				if (itry == 0)
				{
					checkGridL1Err(ALL_POINTS, Wgrid, Rgrid, MAX_ERR_FLOAT);
					checkGridMaxErr(ALL_POINTS, Wgrid, Rgrid, MAX_ERR_FLOAT);
				}
			}

			AddGridGB = Myfloat(nGridPoint * 3 * sizeof(Myfloat) / testCase_time_best / 1.e9);
			AddGridGpoint = Myfloat(nGridPoint / testCase_time_best / 1.e9);
			printInfo(MASTER, " Best achieved GByte/s", AddGridGB);
			printInfo(MASTER, " Best achieved GPoint/s", AddGridGpoint);
		}

		{
			//============================================
			// multiply array
			// W = U * V
			//============================================

			print_blank();
			string caseName = testCaseName + "MultiplyArray";
			printInfo(MASTER, " * Case", caseName);

			Ugrid.fill(ALL_POINTS, a1);
			Vgrid.fill(ALL_POINTS, a2);
			Wgrid.fill(ALL_POINTS, 0.0);
			Rgrid.fill(ALL_POINTS, a1 * a2);

			double testCase_time_best = FLT_MAX;

			Myint ntry = Config::Instance()->ntry;

			for (Myint itry = 0; itry < ntry; itry++)
			{
				double t0 = MPI_Wtime();
				Wgrid.multiplyArray(Ugrid, Vgrid);
				double t1 = MPI_Wtime();

				double testCase_time = t1 - t0;
				Myfloat testCase_bw = nGridPoint * 3 * sizeof(Myfloat) / testCase_time / 1e9;
				printDebug(LIGHT_DEBUG, "Time", testCase_time);
				printDebug(LIGHT_DEBUG, "Speed", testCase_bw);
				testCase_time_best = min(testCase_time_best, testCase_time);

				// check testCase results
				if (itry == 0)
				{
					checkGridL1Err(ALL_POINTS, Wgrid, Rgrid, MAX_ERR_FLOAT);
					checkGridMaxErr(ALL_POINTS, Wgrid, Rgrid, MAX_ERR_FLOAT);
				}
			}

			MultiplyGridGB = Myfloat(nGridPoint * 3 * sizeof(Myfloat) / testCase_time_best / 1.e9);
			MultiplyGridGpoint = Myfloat(nGridPoint / testCase_time_best / 1.e9);
			printInfo(MASTER, " Best achieved GByte/s", MultiplyGridGB);
			printInfo(MASTER, " Best achieved GPoint/s", MultiplyGridGpoint);
		}

		{
			//============================================
			// add and update array
			// W = W + U
			//============================================

			print_blank();
			string caseName = testCaseName + "AddUpdateArray";
			printInfo(MASTER, " * Case", caseName);

			Ugrid.fill(ALL_POINTS, a1);
			Wgrid.fill(ALL_POINTS, a2);
			Rgrid.fill(ALL_POINTS, a1 + a2);

			double testCase_time_best = FLT_MAX;

			Myint ntry = Config::Instance()->ntry;

			for (Myint itry = 0; itry < ntry; itry++)
			{
				double t0 = MPI_Wtime();
				Wgrid.addUpdateArray(Ugrid);
				double t1 = MPI_Wtime();

				double testCase_time = t1 - t0;
				Myfloat testCase_bw = nGridPoint * 3 * sizeof(Myfloat) / testCase_time / 1e9;
				printDebug(LIGHT_DEBUG, "Time", testCase_time);
				printDebug(LIGHT_DEBUG, "Speed", testCase_bw);
				testCase_time_best = min(testCase_time_best, testCase_time);

				// check testCase results
				if (itry == 0)
				{
					checkGridL1Err(ALL_POINTS, Wgrid, Rgrid, MAX_ERR_FLOAT);
					checkGridMaxErr(ALL_POINTS, Wgrid, Rgrid, MAX_ERR_FLOAT);
				}
			}

			AddUpdateGridGB = Myfloat(nGridPoint * 3 * sizeof(Myfloat) / testCase_time_best / 1.e9);
			AddUpdateGridGpoint = Myfloat(nGridPoint / testCase_time_best / 1.e9);
			printInfo(MASTER, " Best achieved GByte/s", AddUpdateGridGB);
			printInfo(MASTER, " Best achieved GPoint/s", AddUpdateGridGpoint);
		}

		// log perf
		if (myMpiRank == 0)
		{
			perfLogFile
				// 10, 11, 12, 13
				<< FillGridGB << " " << FillGridGpoint << " " << CopyGridGB << " " << CopyGridGpoint << " "
				// 14, 15, 16, 17
				<< AddGridGB << " " << AddGridGpoint << " " << MultiplyGridGB << " " << MultiplyGridGpoint << " "
				// 18, 19
				<< AddUpdateGridGB << " " << AddUpdateGridGpoint << " "
				<< "\n";
		}

		this->finalize();

		printDebug(MID_DEBUG, "Out TestCase_Memory::run");
		return (RTN_CODE_OK);
	}

} // namespace hpcscan
