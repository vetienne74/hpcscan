
//-------------------------------------------------------------------------------------------------------
// use this template to create a new test case
//-------------------------------------------------------------------------------------------------------

#include "testCase_Template.h"

#include <algorithm>
#include <cfloat>

#include <mpi.h>
#include <omp.h>

#include "config.h"
#include "constant.h"
#include "grid_Factory.h"
#include "output_report.h"

namespace hpcscan {

TestCase_Template::TestCase_Template(void)
{
	testCaseName    = "Template" ;
	testCaseVersion = "Standard implementation" ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code TestCase_Template::run(void)
{

	printDebug(MID_DEBUG, "In TestCase_Template::run") ;

	if (this->initialize() == RTN_CODE_KO) return (RTN_CODE_KO) ;

	// create a block for each test
	{
		// example check grids have same size
		print_blank() ;
		string caseName = testCaseName + "SameSize" ;
		printInfo(MASTER, " * Case", caseName) ;

		const string gridMode = Config::Instance()->testMode ;
		auto gridIn2  = Grid_Factory::create(gridMode, GRID_GLOBAL) ;
		auto gridOut2 = Grid_Factory::create(gridMode, GRID_GLOBAL) ;
		Grid &gridIn  = *gridIn2 ;
		Grid &gridOut = *gridOut2 ;
		gridIn.initializeGrid() ;
		gridOut.initializeGrid() ;

		bool sameSize = gridIn.sameSize(gridOut) ;
		checkBoolDiff(sameSize, true) ;
	}

	this->finalize() ;

	printDebug(MID_DEBUG, "Out TestCase_Template::run") ;
	return(RTN_CODE_OK) ;
}

} // namespace hpcscan

