
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

//-------------------------------------------------------------------------------------------------------

void Hardware_Cuda::info(void)
{
	printDebug(MID_DEBUG, "IN Hardware_Cuda::info");

	printInfo(MASTER, " Hardware information") ;
	printInfo(MASTER, " NVIDIA GPU") ;

	// display all available GPUs
	int startDevice = 0;
	int endDevice = 0;
	int deviceCount;
	cudaError_t error_id = cudaGetDeviceCount(&deviceCount);
	if (error_id != cudaSuccess) {
		printError(" In Grid_Cuda::info, cudaGetDeviceCount", (int) error_id) ;
	}

	if (deviceCount == 0) {
		printError(" No GPU found") ;
	}
	else
	{
		printInfo(MASTER, " Number of GPUs found", deviceCount) ;
	}

	startDevice = 0;
	endDevice = deviceCount - 1;

	for (int currentDevice = startDevice; currentDevice <= endDevice;
			currentDevice++) {
		cudaDeviceProp deviceProp;
		cudaError_t error_id = cudaGetDeviceProperties(&deviceProp, currentDevice);

		if (error_id == cudaSuccess) {
			string deviceStr = " Device #" + to_string(currentDevice) + "\t";
			printInfo(MASTER, deviceStr, deviceProp.name) ;

			if (deviceProp.computeMode == cudaComputeModeProhibited) {
				printError(" Error: device is running in <Compute Mode Prohibited>") ;
			}
		} else {
			printf("cudaGetDeviceProperties returned %d\n-> %s\n", (int)error_id,
					cudaGetErrorString(error_id));
		}
	}

	if (Config::Instance()->gpuMpiAware)
	{
		printInfo(MASTER, " MPI GPU-Aware Library", "ENABLED") ;
	}
	else
	{
		printInfo(MASTER, " MPI GPU-Aware Library", "DISABLED") ;
	}

	print_line5() ;

	printDebug(MID_DEBUG, "OUT Hardware_Cuda::info");
}

} // namespace hpcscan
