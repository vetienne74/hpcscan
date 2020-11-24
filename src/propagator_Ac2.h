
#ifndef HPCSCAN_PROPAGATOR_AC2_H_
#define HPCSCAN_PROPAGATOR_AC2_H_

#include <memory>

#include "grid.h"
#include "type_def.h"

namespace hpcscan {

class Propagator_Ac2
{
public:

	// constructor
	Propagator_Ac2(void) ;

	// compute wavefield at next time step
	Rtn_code computeWavefieldNextTimeStep(Grid& prnGridIn, Grid& prcGridIn) ;

	// specific computation of this class
	virtual Rtn_code computePressureWithFD(Grid& prnGridIn, Grid& prcGridIn) ;

	// determine stable dt
	Myfloat getCFL(void) ;

	// initialize
	virtual Rtn_code initialize(PropaInit_type propaInitType) ;

	// initialize grid
	virtual Rtn_code initializeGrid(Grid& gridIn, PropaInit_type propaInitType, Myfloat timeSec) ;

	// print info
	void info(void) ;

	// the grids of the propagator
	shared_ptr<Grid> coefGrid ;   // physical coefficient
	shared_ptr<Grid> prnGrid ;    // pressure wavefield new time step
	shared_ptr<Grid> prcGrid ;    // pressure wavefield previous time step

	// FD scheme order in space
	Myint fdOrder ;

	// Number of time steps
	Myint nt ;

	// Time steps
	Myfloat dt, stableDt ;

	// Maximum time
	Myfloat tmax ;

	// Snaphot increment (no. of time steps)
	Myint snapInc ;

	// Min and max velocity in model
	Myfloat minVelocity, maxVelocity ;

	// Max frequency
	Myfloat maxFreq ;

	// Wave equation and kernel type
	string propaWaveEq ;
	string propaKernelType ;

protected:

} ;

} // namespace hpcscan

#endif
