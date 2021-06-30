
//-------------------------------------------------------------------------------------------------------
// This class handles all queries related to the hardware
// Generic class that targets CPUs
// Associated to the test modes Baseline & CacheBlk
//-------------------------------------------------------------------------------------------------------

#ifndef HPCSCAN_HARDWARE_H_
#define HPCSCAN_HARDWARE_H_

#include <string>
#include <vector>

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

	// supports electric power reading
	virtual bool supportGetPowerUsage(void) ;

	// update hardware counters
	// add one entry in hwCounterVec
	void updateHwCounter(void) ;

	// measure current power usage
	Myfloat measureCurrentPower(void) ;

	// display hardware counters statistics
	void displayCounterStat(void) ;

protected:

	// vector of hwCounter_struct_type
	vector<hwCounter_struct_type> hwCounterVec ;

} ;

} // namespace hpcscan

#endif
