
//-------------------------------------------------------------------------------------------------------
// This class handles all queries related to the hardware
// Generic class that targets CPUs
// Associated to the test modes Baseline & CacheBlk
//-------------------------------------------------------------------------------------------------------

#ifndef HPCSCAN_HARDWARE_H_
#define HPCSCAN_HARDWARE_H_

#include <string>
#include <vector>

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

	// update hardware counters
	// add one entry in hwCounterVec
	void updateHwCounter(void) ;

	// update hardware counters at regular time interval
	void watchTimeAndUpdateHwCounter(void) ;

	// measure current power usage
	virtual Myfloat measureCurrentPower(void) ;

	// display hardware counters statistics
	void displayCounterStat(void) ;

	// supports electric power reading
	bool supportGetPowerUsage ;

protected:

	// check if hw supports electric power reading
	virtual bool checkSupportGetPowerUsage(void) ;


	// vector of hwCounter_struct_type
	vector<hwCounter_struct_type> hwCounterVec ;

} ;

} // namespace hpcscan

#endif
