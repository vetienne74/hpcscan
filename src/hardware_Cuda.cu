
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
	printDebug(LIGHT_DEBUG, "IN Hardware_Cuda::Hardware_Cuda");

	supportGetPowerUsage = checkSupportGetPowerUsage() ;

	printDebug(LIGHT_DEBUG, "OUT Hardware_Cuda::Hardware_Cuda");
}

//-------------------------------------------------------------------------------------------------------

Hardware_Cuda::~Hardware_Cuda(void)
{
	printDebug(LIGHT_DEBUG, "IN Hardware_Cuda::~Hardware_Cuda");

	// TODO ~Hardware_Cuda

	printDebug(LIGHT_DEBUG, "OUT Hardware_Cuda::~Hardware_Cuda");
}

//-------------------------------------------------------------------------------------------------------

void Hardware_Cuda::info(void)
{
	printDebug(LIGHT_DEBUG, "IN Hardware_Cuda::info");

	print_blank();
	printInfo(MASTER, " Hardware information") ;

	// display host info
	Hardware::hostInfo() ;

	// display all available devices (GPUs)
	print_blank() ;
	printInfo(MASTER, " Target hardware:") ;
	{
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
	}

	// CUDA aware library
	if (Config::Instance()->gpuMpiAware)
	{
		printInfo(MASTER, " MPI GPU-Aware Library", "ENABLED") ;
	}
	else
	{
		printInfo(MASTER, " MPI GPU-Aware Library", "DISABLED") ;
	}

	// support for power usage
	if (supportGetPowerUsage)
	{
		printInfo(MASTER, " Read power usage", "SUPPORTED") ;
	}
	else
	{
		printInfo(MASTER, " Read power usage", "NOT SUPPORTED") ;
	}
	print_line2() ;

	printDebug(LIGHT_DEBUG, "OUT Hardware_Cuda::info");
}

//-------------------------------------------------------------------------------------------------------

bool Hardware_Cuda::checkSupportGetPowerUsage(void)
{
	printDebug(LIGHT_DEBUG, "IN Hardware_Cuda::checkSupportGetPowerUsage");

	bool retVal = false ;

	// get power consumption
	{
		//nvmlDeviceGetPowerUsage (nvmlDevice_t device, unsigned int* power)
		unsigned int power , i ;
		nvmlReturn_t result;

		// First initialize NVML library
		result = nvmlInit();
		if (NVML_SUCCESS != result)
		{
			printf("Failed to initialize NVML: %s\n", nvmlErrorString(result));
		}

		// TODO adapt to MPI multi-process
		i = 0 ; // get device for rank 0
		result = nvmlDeviceGetHandleByIndex(i, &myDevice);
		if (NVML_SUCCESS != result)
		{
			printf("Failed to get handle for device %u: %s\n", i, nvmlErrorString(result));
		}

		result = nvmlDeviceGetPowerUsage(myDevice, &power) ;
		if (result != NVML_SUCCESS)
		{
			printDebug(LIGHT_DEBUG, " nvmlDeviceGetPowerUsage", (const char*) nvmlErrorString(result)) ;
		}
		else
		{
			printInfo(MASTER, " Current power (mWatt)", (Myfloat) power / 1000.) ;
			retVal = true ;
		}
	}

	printDebug(LIGHT_DEBUG, "OUT Hardware_Cuda::checkSupportGetPowerUsage");
	return(retVal) ;
}

//-------------------------------------------------------------------------------------------------------

Myfloat Hardware_Cuda::measureCurrentPower(void)
{
	printDebug(LIGHT_DEBUG, "IN Hardware_Cuda::measureCurrentPower");

	Myfloat retVal = UNSPECIFIED ;
	if (supportGetPowerUsage)
	{
		unsigned int power ;
		nvmlDeviceGetPowerUsage(myDevice, &power) ;
		retVal = (Myfloat) power / 1000. ;
	}

	printDebug(LIGHT_DEBUG, "OUT Hardware_Cuda::measureCurrentPower");
	return(retVal) ;
}


} // namespace hpcscan
