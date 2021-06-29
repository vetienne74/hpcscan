
//-------------------------------------------------------------------------------------------------------
// This class handles all queries related to the hardware
// Generic class that targets CPUs
// Associated to the test modes Baseline & CacheBlk
//-------------------------------------------------------------------------------------------------------

#ifndef HPCSCAN_HARDWARE_H_
#define HPCSCAN_HARDWARE_H_

#include <string>

#include "mpi.h"

#include "type_def.h"

namespace hpcscan {

class Hardware
{
public:

	// constructor
	Hardware(string gridMode) ;

	// destructor
	~Hardware(void) ;

	// print info
	virtual void info(void) ;

protected:

} ;

} // namespace hpcscan

#endif
