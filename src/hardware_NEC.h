
//-------------------------------------------------------------------------------------------------------
// This class handles all queries related to the hardware
// Specialized class that targets NEC SX-Aurora
// Derived class from Hardware
// Associated to all test modes when hpcscan is compiled with NEC compiler
//-------------------------------------------------------------------------------------------------------

#ifndef HPCSCAN_HARDWARE_NEC_H_
#define HPCSCAN_HARDWARE_NEC_H_

#include <string>

#include "mpi.h"

#include "hardware.h"
#include "type_def.h"

namespace hpcscan {

class Hardware_NEC : public Hardware
{
public:

	// constructor
	Hardware_NEC(string gridMode) ;

	// destructor
	~Hardware_NEC(void) ;

	// print info
	virtual void info(void) ;

	// measure current power usage
	virtual Myfloat measureCurrentPower(void) ;

protected:

	// check if hw supports electric power reading
	virtual bool checkSupportGetPowerUsage(void) ;

	// device id
	// TODO adapt to multiple devices
	int deviceId;

} ;

} // namespace hpcscan

#endif
