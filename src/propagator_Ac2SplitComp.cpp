
//-------------------------------------------------------------------------------------------------------
// * 2nd order acoustic wave equation with constant density = 1
//
// ### This class derives from propagator_Ac2 ###
//
// Computations are split:
// First computation of the Laplacian with operator Grid::FD_LAPLACIAN
// Second update pressure with with Grid::updatePressure
//-------------------------------------------------------------------------------------------------------

#include "propagator_Ac2SplitComp.h"

#include "config.h"
#include "constant.h"
#include "fdm.h"
#include "grid_Factory.h"
#include "output_report.h"

using namespace std;

namespace hpcscan {

//-------------------------------------------------------------------------------------------------------

Propagator_Ac2SplitComp::Propagator_Ac2SplitComp(void) : Propagator_Ac2()
																																														{
	printDebug(FULL_DEBUG, "In Propagator_Ac2SplitComp::Propagator_Ac2SplitComp") ;

	propaKernelType = "With Laplacian operator" ;

	laplaGrid   = nullptr ;

	printDebug(FULL_DEBUG, "Out Propagator_Ac2SplitComp::Propagator_Ac2SplitComp") ;
																																														}

//-------------------------------------------------------------------------------------------------------

Rtn_code Propagator_Ac2SplitComp::initialize(PropaInit_type propaInitType)
{

	printDebug(FULL_DEBUG, "In Propagator_Ac2SplitComp::initialize") ;

	// initialize parent
	Propagator_Ac2::initialize(propaInitType) ;

	// allocate Laplacian grid
	const string gridMode = Config::Instance()->testMode ;
	laplaGrid = Grid_Factory::create(gridMode, GRID_LOCAL) ;
	laplaGrid->initializeGrid();

	if (propaInitType == EIGEN_MODE)
	{
		laplaGrid->defineUnitGrid() ;
	}
	else
	{
		printError("In Propagator_Ac2SplitComp::initialize, propaInitType not supported", propaInitType) ;
		return(RTN_CODE_KO) ;
	}

	printDebug(FULL_DEBUG, "Out Propagator_Ac2SplitComp::initialize") ;
	return(RTN_CODE_OK) ;
}


//-------------------------------------------------------------------------------------------------------

Rtn_code Propagator_Ac2SplitComp::computePressureWithFD(Grid& prnGridIn, Grid& prcGridIn)
{

	printDebug(FULL_DEBUG, "In Propagator_Ac2SplitComp::computePressureWithFD") ;

	// compute Laplacian
	prcGridIn.FD_LAPLACIAN(INNER_POINTS, *laplaGrid, fdOrder) ;

	// update pressure
	prnGridIn.updatePressure(INNER_POINTS, prcGridIn, *coefGrid, *laplaGrid) ;

	printDebug(FULL_DEBUG, "Out Propagator_Ac2SplitComp::computePressureWithFD") ;
	return(RTN_CODE_OK) ;
}

} // namespace hpcscan
