
//-------------------------------------------------------------------------------------------------------
// This class handles all queries related to the hardware
// Specialized class that targets NVIDIA GPUs
// Derived class from Hardware
// Associated to the test modes CUDA, CUDA_Ref & CUDA_Opt
//-------------------------------------------------------------------------------------------------------

#ifndef HPCSCAN_HARDWARE_CUDA_H_
#define HPCSCAN_HARDWARE_CUDA_H_

#include <string>

#include "mpi.h"
#include <nvml.h>

#include "hardware.h"
#include "type_def.h"

namespace hpcscan {

class Hardware_Cuda : public Hardware
{
public:

	// constructor
	Hardware_Cuda(string gridMode) ;

	// destructor
	~Hardware_Cuda(void) ;

	// print info
	virtual void info(void) ;

	// measure current power usage
	virtual Myfloat measureCurrentPower(void) ;

protected:

	// check if hw supports electric power reading
	virtual bool checkSupportGetPowerUsage(void) ;

	// device token
	// TODO adapt to multiple devices
	nvmlDevice_t myDevice;

} ;

} // namespace hpcscan

#endif
