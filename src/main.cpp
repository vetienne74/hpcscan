
#include <mpi.h>

#include "config.h"
#include "global.h"
#include "output_report.h"
#include "testCase_Comm.h"
#include "testCase_FD_D2.h"
#include "testCase_Grid.h"
#include "testCase_Memory.h"
#include "testCase_Propa.h"
#include "testCase_Template.h"
#include "testCase_Util.h"
#include "type_def.h"

//-------------------------------------------------------------------------------------------------------
// global variables
//-------------------------------------------------------------------------------------------------------

// debug level switch
hpcscan::Debug_level hpcscan::debug = NO_DEBUG ;

// global number of MPI process
int hpcscan::nMpiProc ;

// global rank of MPI process
int hpcscan::myMpiRank ;

int main(int argc, char* argv[])
{
	// start MPI environment
	MPI_Init(&argc, &argv) ;
	MPI_Comm_size(MPI_COMM_WORLD, &hpcscan::nMpiProc) ;
	MPI_Comm_rank(MPI_COMM_WORLD, &hpcscan::myMpiRank) ;

	// start timer
	double tStart = MPI_Wtime() ;

	//-------------------------------------------------------------------------------------------------------
	// parse arguments in command line
	// and initialize Config singleton with all parameters
	//-------------------------------------------------------------------------------------------------------
	hpcscan::Rtn_code rtn_code ;
	rtn_code = hpcscan::Config::Instance()->parse_argument(argc, argv) ;
	if (rtn_code != hpcscan::RTN_CODE_OK) {
		hpcscan::Config::Instance()->finalize(rtn_code) ;
		return(hpcscan::RTN_CODE_OK) ;
	}

	// other initializations
	hpcscan::Config::Instance()->initialize() ;

	//-------------------------------------------------------------------------------------------------------
	// print header of output report
	//-------------------------------------------------------------------------------------------------------
	rtn_code = hpcscan::print_header_of_output_report() ;
	if (rtn_code != hpcscan::RTN_CODE_OK) {
		hpcscan::Config::Instance()->finalize(rtn_code) ;
		return(hpcscan::RTN_CODE_OK) ;
	}

	//-------------------------------------------------------------------------------------------------------
	// print configuration parameters
	//-------------------------------------------------------------------------------------------------------
	rtn_code = hpcscan::Config::Instance()->info() ;
	if (rtn_code != hpcscan::RTN_CODE_OK) {
		hpcscan::Config::Instance()->finalize(rtn_code) ;
		return(hpcscan::RTN_CODE_OK) ;
	}

	//-------------------------------------------------------------------------------------------------------
	// print hardware info
	//-------------------------------------------------------------------------------------------------------
	if (hpcscan::Config::Instance()->hw != nullptr)
		hpcscan::Config::Instance()->hw->info() ;

	//-------------------------------------------------------------------------------------------------------
	// run testCases
	//-------------------------------------------------------------------------------------------------------
	{
		hpcscan::TestCase_Comm testCase ;
		testCase.run() ;
	}
	{
		hpcscan::TestCase_FD_D2 testCase ;
		testCase.run() ;
	}
	{
		hpcscan::TestCase_Grid testCase ;
		testCase.run() ;
	}
	{
		hpcscan::TestCase_Memory testCase ;
		testCase.run() ;
	}
	{
		hpcscan::TestCase_Propa testCase ;
		testCase.run() ;
	}
	{
		hpcscan::TestCase_Template testCase ;
		testCase.run() ;
	}
	{
		hpcscan::TestCase_Util testCase ;
		testCase.run() ;
	}

	// end timer
	double tEnd = MPI_Wtime() ;
	hpcscan::printInfo(hpcscan::MASTER, " Total run time (s)", tEnd-tStart) ;

	// terminate program
	hpcscan::Config::Instance()->finalize(hpcscan::RTN_CODE_OK) ;

	return(0);
}
