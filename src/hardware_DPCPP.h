
//-------------------------------------------------------------------------------------------------------
// This class handles all queries related to the hardware
// Specialized class that targets INTEL CPUs, GPUs and FPGAs
// Derived class from Hardware
// Associated to the test modes DPCPP
//-------------------------------------------------------------------------------------------------------

#ifndef HPCSCAN_HARDWARE_DPCPP_H_
#define HPCSCAN_HARDWARE_DPCPP_H_

#include <string>

#include "mpi.h"

#include "hardware.h"
#include "type_def.h"

namespace hpcscan {

class Hardware_DPCPP : public Hardware
{
public:

	// constructor
	Hardware_DPCPP(string gridMode) ;

	// destructor
	~Hardware_DPCPP(void) ;

	// print info
	virtual void info(void) ;

	// measure current power usage
	virtual Myfloat measureCurrentPower(void) ;

protected:

	// check if hw supports electric power reading
	virtual bool checkSupportGetPowerUsage(void) ;

} ;

} // namespace hpcscan

#endif
