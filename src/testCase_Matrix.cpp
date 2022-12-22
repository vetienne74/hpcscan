
//-------------------------------------------------------------------------------------------------------
// test cases to measure the matrix operations bandwidth
// baseline kernel - no optimization
//-------------------------------------------------------------------------------------------------------

#include "testCase_Matrix.h"

#include <algorithm> // for min()
#include <cfloat>	 // for FLT_MAX
#include <iostream>

#include <mpi.h>
#include <omp.h>

#include "config.h"
#include "constant.h"
#include "global.h"
#include "grid_Factory.h"
#include "output_report.h"

#ifdef _FTRACE
#include <ftrace.h>
#define FTRACE_BEGIN(a) ftrace_region_begin(a);
#define FTRACE_END(a) ftrace_region_end(a);
#else
#define FTRACE_BEGIN(a)
#define FTRACE_END(a)
#endif

namespace hpcscan
{

	TestCase_Matrix::TestCase_Matrix(void)
	{
		testCaseName = "Matrix";
		testCaseVersion = "Standard implementation";
	}

	//-------------------------------------------------------------------------------------------------------

	Rtn_code TestCase_Matrix::run(void)
	{

		printDebug(MID_DEBUG, "In TestCase_Matrix::run");

		if (this->initialize() == RTN_CODE_KO)
			return (RTN_CODE_KO);

		const Myint nMatrix = 5;
		const Myint nElem = 100000;

		// one block per case
		{
			//============================================
			// Copy Vector
			//============================================

			print_blank();
			string caseName = testCaseName + "CopyVector";
			printInfo(MASTER, " * Case", caseName);

			printInfo(MASTER, " Number of elements", nElem);

			double testCase_time_best = FLT_MAX;

			Myint ntry = Config::Instance()->ntry;

			// const Myint nNodePerElem = nMatrix * nMatrix * nMatrix;
			Myint nNodePerElem = 512;
			printInfo(MASTER, " Nnode per element", nNodePerElem);
			const Myint nNode = nElem * nNodePerElem;
			Myfloat *globData = new Myfloat[nNode];
			Myfloat *B = new Myfloat[nNode];

			// Myint *globId = new Myint[nNode];

			// initialize globIb & globData
			for (Myint iNode = 0; iNode < nNode; iNode++)
			{
				// globId[iNode] = iNode;
				globData[iNode] = 1.0;
			}

			for (Myint itry = 0; itry < ntry; itry++)
			{
				double t0 = MPI_Wtime();
				FTRACE_BEGIN(caseName.c_str())

				// loop over all elements
#pragma omp parallel for
				for (Myint iElem = 0; iElem < nElem; iElem++)
				{

					// copy global data to local data
#pragma _NEC packed_stencil
					for (Myint ijk = 0; ijk < nNodePerElem; ijk++)
					{
						Myint iGlob = iElem * nNodePerElem + ijk;
						B[iGlob] = globData[iGlob];
					}

				} // for (Myint iElem = 0; iElem < nElem; iElem++)

				FTRACE_END(caseName.c_str())
				double t1 = MPI_Wtime();

				double testCase_time = t1 - t0;

				Myfloat testCase_bw = nNode * 2 * sizeof(Myfloat) / testCase_time / 1e9;
				printDebug(LIGHT_DEBUG, "Time", testCase_time);
				printInfo(MASTER, " Time", testCase_time);
				printDebug(LIGHT_DEBUG, "Speed", testCase_bw);
				testCase_time_best = min(testCase_time_best, testCase_time);

				// check testCase results
				if (itry == 0)
				{
					// check globData
					Myfloat minVal = FLT_MAX;
					Myfloat maxVal = -FLT_MAX;
					for (Myint iNode = 0; iNode < nNode; iNode++)
					{
						if (B[iNode] < minVal)
							minVal = globData[iNode];
						if (B[iNode] > maxVal)
							maxVal = globData[iNode];
					}

					checkFloatDiff(minVal, 1.0, MAX_ERR_FLOAT);
					checkFloatDiff(maxVal, 1.0, MAX_ERR_FLOAT);
				}
			}

			Myfloat GB = Myfloat(nNode * 2 * sizeof(Myfloat) / testCase_time_best / 1.e9);
			Myfloat GNode = Myfloat(nNode / testCase_time_best / 1.e9);
			// Myfloat GFlops = GNode * 9;
			printInfo(MASTER, " Best achieved GByte/s", GB);
			printInfo(MASTER, " Best achieved GNode/s", GNode);
			// printInfo(MASTER, " Best achieved GFlop/s", GFlops);

			// delete[] A;
			delete[] B;
			// delete[] C;
			delete[] globData;
			// delete[] globId;
		}

		{
			//============================================
			// Copy Vector + Indirection
			//============================================

			print_blank();
			string caseName = testCaseName + "CopyVectorWithIndirection";
			printInfo(MASTER, " * Case", caseName);

			printInfo(MASTER, " Number of elements", nElem);

			double testCase_time_best = FLT_MAX;

			Myint ntry = Config::Instance()->ntry;

			// const Myint nNodePerElem = nMatrix * nMatrix * nMatrix;
			const Myint nNodePerElem = 512;
			printInfo(MASTER, " Nnode per element", nNodePerElem);
			const Myint nNode = nElem * nNodePerElem;
			Myfloat *globData = new Myfloat[nNode];
			Myfloat *B = new Myfloat[nNode];
			Myint *globId = new Myint[nNode];

			// initialize globIb & globData
			for (Myint iNode = 0; iNode < nNode; iNode++)
			{
				globId[iNode] = iNode;
				globData[iNode] = 1.0;
			}

			for (Myint itry = 0; itry < ntry; itry++)
			{
				double t0 = MPI_Wtime();
				FTRACE_BEGIN(caseName.c_str())

				// loop over all elements
#pragma _NEC novector
#pragma omp parallel for
				for (Myint iElem = 0; iElem < nElem; iElem++)
				{
#pragma _NEC ivdep
// copy global data to local data
#pragma _NEC packed_stencil
					for (Myint ijk = 0; ijk < nNodePerElem; ijk++)
					{
						// Myint iGlob = globId[iElem * nNodePerElem + ijk];
						Myint iGlob = iElem * nNodePerElem + globId[ijk];
						// Myint iGlob = iElem * nNodePerElem + ijk;
						B[iElem * nNodePerElem + ijk] = globData[iGlob];
					}

				} // for (Myint iElem = 0; iElem < nElem; iElem++)

				FTRACE_END(caseName.c_str())
				double t1 = MPI_Wtime();

				double testCase_time = t1 - t0;

				Myfloat testCase_bw = nNode * 2 * sizeof(Myfloat) / testCase_time / 1e9;
				printDebug(LIGHT_DEBUG, "Time", testCase_time);
				printInfo(MASTER, " Time", testCase_time);
				printDebug(LIGHT_DEBUG, "Speed", testCase_bw);
				testCase_time_best = min(testCase_time_best, testCase_time);

				// check testCase results
				if (itry == 0)
				{
					// check globData
					Myfloat minVal = FLT_MAX;
					Myfloat maxVal = -FLT_MAX;
					for (Myint iNode = 0; iNode < nNode; iNode++)
					{
						if (B[iNode] < minVal)
							minVal = globData[iNode];
						if (B[iNode] > maxVal)
							maxVal = globData[iNode];
					}

					checkFloatDiff(minVal, 1.0, MAX_ERR_FLOAT);
					checkFloatDiff(maxVal, 1.0, MAX_ERR_FLOAT);
				}
			}

			Myfloat GB = Myfloat(nNode * 2 * sizeof(Myfloat) / testCase_time_best / 1.e9);
			Myfloat GNode = Myfloat(nNode / testCase_time_best / 1.e9);
			// Myfloat GFlops = GNode * 9;
			printInfo(MASTER, " Best achieved GByte/s", GB);
			printInfo(MASTER, " Best achieved GNode/s", GNode);
			// printInfo(MASTER, " Best achieved GFlop/s", GFlops);

			// delete[] A;
			delete[] B;
			// delete[] C;
			delete[] globData;
			// delete[] globId;
		}

		{
			//============================================
			// C = A * B
			// Algo1 (initial version)
			//============================================

			print_blank();
			string caseName = testCaseName + "C=AxB_algo1";
			printInfo(MASTER, " * Case", caseName);

			string sizeA = to_string(nMatrix) + " x " + to_string(nMatrix);
			printInfo(MASTER, " A matrix size\t", sizeA);
			string sizeB = to_string(nMatrix * nMatrix) + " x " + to_string(nMatrix);
			printInfo(MASTER, " B matrix size\t", sizeB);
			string sizeC = to_string(nMatrix * nMatrix) + " x " + to_string(nMatrix);
			printInfo(MASTER, " C matrix size\t", sizeC);
			printInfo(MASTER, " Number of elements", nElem);

			double testCase_time_best = FLT_MAX;

			Myint ntry = Config::Instance()->ntry;

			const Myint nNodePerElem = nMatrix * nMatrix * nMatrix;
			printInfo(MASTER, " Nnode per element", nNodePerElem);
			Myfloat *A = new Myfloat[nMatrix * nMatrix];
			Myfloat *B = new Myfloat[nNodePerElem];
			Myfloat *C = new Myfloat[nNodePerElem];
			const Myint nNode = nElem * nNodePerElem;
			Myfloat *globData = new Myfloat[nNode];
			Myint *globId = new Myint[nNode];

			// initialize globIb & globData
			for (Myint iNode = 0; iNode < nNode; iNode++)
			{
				globId[iNode] = iNode;
				globData[iNode] = 1.0;
			}

			// initialize matrix A
			for (Myint ii = 0; ii < nMatrix * nMatrix; ii++)
			{
				A[ii] = 1.0;
			}

			for (Myint itry = 0; itry < ntry; itry++)
			{
				double t0 = MPI_Wtime();
				FTRACE_BEGIN(caseName.c_str())

				// loop over all elements
				for (Myint iElem = 0; iElem < nElem; iElem++)
				{

					// copy global data to local data
					for (Myint ijk = 0; ijk < nNodePerElem; ijk++)
					{
						Myint iGlob = globId[iElem * nNodePerElem + ijk];
						B[ijk] = globData[iGlob];
					}

					// matrix operations
					for (Myint jj = 0; jj < nMatrix * nMatrix; jj++)
					{
						for (Myint ii = 0; ii < nMatrix; ii++)
						{
							Myint iiA = ii * nMatrix;
							Myint iiC = ii * nMatrix * nMatrix;

							C[iiC + jj] =
								A[iiA] * B[jj] + A[iiA + 1] * B[nMatrix * nMatrix + jj] + A[iiA + 2] * B[nMatrix * nMatrix * 2 + jj] + A[iiA + 3] * B[nMatrix * nMatrix * 3 + jj] + A[iiA + 4] * B[nMatrix * nMatrix * 4 + jj];
						}
					}

					// copy local data to global data
					for (Myint ijk = 0; ijk < nNodePerElem; ijk++)
					{
						Myint iGlob = globId[iElem * nNodePerElem + ijk];
						globData[iGlob] = C[ijk];
					}

				} // for (Myint iElem = 0; iElem < nElem; iElem++)

				FTRACE_END(caseName.c_str())
				double t1 = MPI_Wtime();

				double testCase_time = t1 - t0;

				Myfloat testCase_bw = nNode * sizeof(Myfloat) / testCase_time / 1e9;
				printDebug(LIGHT_DEBUG, "Time", testCase_time);
				printInfo(MASTER, " Time", testCase_time);
				printDebug(LIGHT_DEBUG, "Speed", testCase_bw);
				testCase_time_best = min(testCase_time_best, testCase_time);

				// check testCase results
				if (itry == 0)
				{
					// check globData
					Myfloat minVal = FLT_MAX;
					Myfloat maxVal = -FLT_MAX;
					for (Myint iNode = 0; iNode < nNode; iNode++)
					{
						if (globData[iNode] < minVal)
							minVal = globData[iNode];
						if (globData[iNode] > maxVal)
							maxVal = globData[iNode];
					}

					checkFloatDiff(minVal, 5.0, MAX_ERR_FLOAT);
					checkFloatDiff(maxVal, 5.0, MAX_ERR_FLOAT);
				}
			}

			Myfloat GB = Myfloat(nNode * 4 * sizeof(Myfloat) / testCase_time_best / 1.e9);
			Myfloat GNode = Myfloat(nNode / testCase_time_best / 1.e9);
			Myfloat GFlops = GNode * 9;
			printInfo(MASTER, " Best achieved GByte/s", GB);
			printInfo(MASTER, " Best achieved GNode/s", GNode);
			printInfo(MASTER, " Best achieved GFlop/s", GFlops);

			delete[] A;
			delete[] B;
			delete[] C;
			delete[] globData;
			delete[] globId;
		}

		{
			//============================================
			// C = A * B
			// Algo2 (optimzied version)
			//============================================

			print_blank();
			string caseName = testCaseName + "C=AxB_algo2";
			printInfo(MASTER, " * Case", caseName);

			string sizeA = to_string(nMatrix) + " x " + to_string(nMatrix);
			printInfo(MASTER, " A matrix size\t", sizeA);
			string sizeB = to_string(nMatrix * nMatrix) + " x " + to_string(nMatrix);
			printInfo(MASTER, " B matrix size\t", sizeB);
			string sizeC = to_string(nMatrix * nMatrix) + " x " + to_string(nMatrix);
			printInfo(MASTER, " C matrix size\t", sizeC);
			printInfo(MASTER, " Number of elements", nElem);

			double testCase_time_best = FLT_MAX;

			Myint ntry = Config::Instance()->ntry;

			const Myint nNodePerElem = nMatrix * nMatrix * nMatrix;
			Myfloat *A = new Myfloat[nMatrix * nMatrix];
			const Myint nNode = nElem * nNodePerElem;
			Myfloat *B = new Myfloat[nNode];
			Myfloat *C = new Myfloat[nNode];
			Myfloat *globData = new Myfloat[nNode];
			Myint *globId = new Myint[nNode];

			// initialize globIb & globData
			for (Myint iNode = 0; iNode < nNode; iNode++)
			{
				globId[iNode] = iNode;
				globData[iNode] = 1.0;
			}

			// initialize matrix A
			for (Myint ii = 0; ii < nMatrix * nMatrix; ii++)
			{
				A[ii] = 1.0;
			}

			for (Myint itry = 0; itry < ntry; itry++)
			{
				double t0 = MPI_Wtime();
				FTRACE_BEGIN(caseName.c_str())

				// loop over all elements
				for (Myint iElem = 0; iElem < nElem; iElem++)
				{

					// copy global data to local data
					for (Myint ijk = 0; ijk < nNodePerElem; ijk++)
					{
						Myint iGlob = globId[iElem * nNodePerElem + ijk];
						B[iElem * nNodePerElem + ijk] = globData[iGlob];
					}
				}

				// loop over jj
				for (Myint jj = 0; jj < nMatrix * nMatrix; jj++)
				{

					// loop over all elements
					for (Myint iElem = 0; iElem < nElem; iElem++)
					{
// matrix operations
#pragma _NEC outerloop_unroll(nMatrix)
						for (Myint ii = 0; ii < nMatrix; ii++)
						{
							Myint iiA = ii * nMatrix;
							Myint iiC = ii * nMatrix * nMatrix;

							C[iElem * nNodePerElem + iiC + jj] =
								A[iiA] * B[iElem * nNodePerElem + jj] + A[iiA + 1] * B[iElem * nNodePerElem + nMatrix * nMatrix + jj] + A[iiA + 2] * B[iElem * nNodePerElem + nMatrix * nMatrix * 2 + jj] + A[iiA + 3] * B[iElem * nNodePerElem + nMatrix * nMatrix * 3 + jj] + A[iiA + 4] * B[iElem * nNodePerElem + nMatrix * nMatrix * 4 + jj];
						}
					}
				} // for (Myint jj = 0; jj < nMatrix * nMatrix; jj++)

				// loop over all elements
				for (Myint iElem = 0; iElem < nElem; iElem++)
				{
					// copy local data to global data
					for (Myint ijk = 0; ijk < nNodePerElem; ijk++)
					{
						Myint iGlob = globId[iElem * nNodePerElem + ijk];
						globData[iGlob] = C[ijk];
					}

				} // for (Myint iElem = 0; iElem < nElem; iElem++)

				FTRACE_END(caseName.c_str())
				double t1 = MPI_Wtime();

				double testCase_time = t1 - t0;

				Myfloat testCase_bw = nNode * sizeof(Myfloat) / testCase_time / 1e9;
				printDebug(LIGHT_DEBUG, "Time", testCase_time);
				printInfo(MASTER, " Time", testCase_time);
				printDebug(LIGHT_DEBUG, "Speed", testCase_bw);
				testCase_time_best = min(testCase_time_best, testCase_time);

				// check testCase results
				if (itry == 0)
				{
					// check globData
					Myfloat minVal = FLT_MAX;
					Myfloat maxVal = -FLT_MAX;
					for (Myint iNode = 0; iNode < nNode; iNode++)
					{
						if (globData[iNode] < minVal)
							minVal = globData[iNode];
						if (globData[iNode] > maxVal)
							maxVal = globData[iNode];
					}

					checkFloatDiff(minVal, 5.0, MAX_ERR_FLOAT);
					checkFloatDiff(maxVal, 5.0, MAX_ERR_FLOAT);
				}
			}

			Myfloat GB = Myfloat(nNode * 4 * sizeof(Myfloat) / testCase_time_best / 1.e9);
			Myfloat GNode = Myfloat(nNode / testCase_time_best / 1.e9);
			Myfloat GFlops = GNode * 9;
			printInfo(MASTER, " Best achieved GByte/s", GB);
			printInfo(MASTER, " Best achieved GNode/s", GNode);
			printInfo(MASTER, " Best achieved GFlop/s", GFlops);

			delete[] A;
			delete[] B;
			delete[] C;
			delete[] globData;
			delete[] globId;
		}

		// log perf
		if (myMpiRank == 0)
		{
			perfLogFile
				// 10, 11, 12, 13
				//<< FillGridGB << " " << FillGridGpoint << " " << CopyGridGB << " " << CopyGridGpoint << " "
				// 14, 15, 16, 17
				//<< AddGridGB << " " << AddGridGpoint << " " << MultiplyGridGB << " " << MultiplyGridGpoint << " "
				// 18, 19
				//<< AddUpdateGridGB << " " << AddUpdateGridGpoint << " "
				<< "\n";
		}

		this->finalize();

		printDebug(MID_DEBUG, "Out TestCase_Matrix::run");
		return (RTN_CODE_OK);
	}

} // namespace hpcscan
