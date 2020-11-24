
#ifndef HPCSCAN_PROPAGATOR_AC2_SPLIT_COMP_H_
#define HPCSCAN_PROPAGATOR_AC2_SPLIT_COMP_H_

#include <memory>

#include "grid.h"
#include "propagator_Ac2.h"
#include "type_def.h"

namespace hpcscan {

class Propagator_Ac2SplitComp : public Propagator_Ac2
{
public:

	// constructor
	Propagator_Ac2SplitComp(void) ;

	// specific computation of this class
	virtual Rtn_code computePressureWithFD(Grid& prnGridIn, Grid& prcGridIn) ;

	// initialize
	virtual Rtn_code initialize(PropaInit_type propaInitType) ;

private:

	// the specific grids of this propagator
	shared_ptr<Grid> laplaGrid ;  // laplacian of pressure

} ;

} // namespace hpcscan

#endif
