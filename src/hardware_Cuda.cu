
//-------------------------------------------------------------------------------------------------------
// This class handles all queries related to the hardware
// Specialized class that targets NVIDIA GPUs
// Derived class from Hardware
// Associated to the test modes CUDA, CUDA_Ref & CUDA_Opt
//-------------------------------------------------------------------------------------------------------

#include "hardware_Cuda.h"

#include "config.h"
#include "constant.h"
#include "global.h"
#include "output_report.h"

using namespace std;

namespace hpcscan {

//-------------------------------------------------------------------------------------------------------

Hardware_Cuda::Hardware_Cuda(string gridMode) : Hardware(gridMode)
{
	printDebug(MID_DEBUG, "IN Hardware_Cuda::Hardware_Cuda");

	// TODO

	printDebug(MID_DEBUG, "OUT Hardware_Cuda::Hardware_Cuda");
}

//-------------------------------------------------------------------------------------------------------

Hardware_Cuda::~Hardware_Cuda(void)
{
	printDebug(MID_DEBUG, "IN Hardware_Cuda::~Hardware_Cuda");

	// TODO

	printDebug(MID_DEBUG, "OUT Hardware_Cuda::~Hardware_Cuda");
}


} // namespace hpcscan
