
#include <mpi.h>

#include "config.h"
#include "global.h"
#include "output_report.h"
#include "testCase_Comm.h"
#include "testCase_FD_D2.h"
#include "testCase_Grid.h"
#include "testCase_Matrix.h"
#include "testCase_Memory.h"
#include "testCase_Modeling.h"
#include "testCase_Propa.h"
#include "testCase_Template.h"
#include "testCase_Util.h"
#include "type_def.h"

//-------------------------------------------------------------------------------------------------------
// global variables
//-------------------------------------------------------------------------------------------------------

// debug level switch
hpcscan::Debug_level hpcscan::debug = NO_DEBUG;

// global number of MPI process
int hpcscan::nMpiProc;

// global rank of MPI process
int hpcscan::myMpiRank;

int main(int argc, char *argv[])
{
	// start MPI environment
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &hpcscan::nMpiProc);
	MPI_Comm_rank(MPI_COMM_WORLD, &hpcscan::myMpiRank);

	// start timer
	double tStart = MPI_Wtime();

	//-------------------------------------------------------------------------------------------------------
	// parse arguments in command line
	// and initialize Config singleton with all parameters
	//-------------------------------------------------------------------------------------------------------
	hpcscan::Rtn_code rtn_code;
	rtn_code = hpcscan::Config::Instance()->parse_argument(argc, argv);
	if (rtn_code != hpcscan::RTN_CODE_OK)
	{
		hpcscan::Config::Instance()->finalize(rtn_code);
		return (hpcscan::RTN_CODE_OK);
	}

	// other initializations
	hpcscan::Config::Instance()->initialize();

	//-------------------------------------------------------------------------------------------------------
	// print header of output report
	//-------------------------------------------------------------------------------------------------------
	rtn_code = hpcscan::print_header_of_output_report();
	if (rtn_code != hpcscan::RTN_CODE_OK)
	{
		hpcscan::Config::Instance()->finalize(rtn_code);
		return (hpcscan::RTN_CODE_OK);
	}

	//-------------------------------------------------------------------------------------------------------
	// print configuration parameters
	//-------------------------------------------------------------------------------------------------------
	rtn_code = hpcscan::Config::Instance()->info();
	if (rtn_code != hpcscan::RTN_CODE_OK)
	{
		hpcscan::Config::Instance()->finalize(rtn_code);
		return (hpcscan::RTN_CODE_OK);
	}

	//-------------------------------------------------------------------------------------------------------
	// run testCases
	//-------------------------------------------------------------------------------------------------------
	hpcscan::Myint numTestCase = 0;
	{
		hpcscan::TestCase_Comm testCase;
		rtn_code = testCase.run();
		if (rtn_code == hpcscan::RTN_CODE_OK)
			numTestCase++;
	}
	{
		hpcscan::TestCase_FD_D2 testCase;
		rtn_code = testCase.run();
		if (rtn_code == hpcscan::RTN_CODE_OK)
			numTestCase++;
	}
	{
		hpcscan::TestCase_Grid testCase;
		rtn_code = testCase.run();
		if (rtn_code == hpcscan::RTN_CODE_OK)
			numTestCase++;
	}
	{
		hpcscan::TestCase_Memory testCase;
		rtn_code = testCase.run();
		if (rtn_code == hpcscan::RTN_CODE_OK)
			numTestCase++;
	}
	{
		hpcscan::TestCase_Matrix testCase;
		rtn_code = testCase.run();
		if (rtn_code == hpcscan::RTN_CODE_OK)
			numTestCase++;
	}
	{
		hpcscan::TestCase_Modeling testCase;
		rtn_code = testCase.run();
		if (rtn_code == hpcscan::RTN_CODE_OK)
			numTestCase++;
	}
	{
		hpcscan::TestCase_Propa testCase;
		rtn_code = testCase.run();
		if (rtn_code == hpcscan::RTN_CODE_OK)
			numTestCase++;
	}
	{
		hpcscan::TestCase_Template testCase;
		rtn_code = testCase.run();
		if (rtn_code == hpcscan::RTN_CODE_OK)
			numTestCase++;
	}
	{
		hpcscan::TestCase_Util testCase;
		rtn_code = testCase.run();
		if (rtn_code == hpcscan::RTN_CODE_OK)
			numTestCase++;
	}

	if (numTestCase == 0)
	{
		hpcscan::printWarning(" No test case was executed");
	}
	else
	{
		hpcscan::printInfo(hpcscan::MASTER, " Number of test cases executed", numTestCase);
	}

	// end timer
	double tEnd = MPI_Wtime();
	hpcscan::printInfo(hpcscan::MASTER, " Total run time (s)", tEnd - tStart);

	// terminate program
	hpcscan::Config::Instance()->finalize(hpcscan::RTN_CODE_OK);

	return (0);
}
