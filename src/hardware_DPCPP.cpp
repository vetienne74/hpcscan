
//-------------------------------------------------------------------------------------------------------
// This class handles all queries related to the hardware
// Specialized class that targets INTEL CPUs, GPUs and FPGAs
// Derived class from Hardware
// Associated to the test mode DPCPP
//-------------------------------------------------------------------------------------------------------

#include "hardware_DPCPP.h"

#include <CL/sycl.hpp>

#include "config.h"
#include "constant.h"
#include "global.h"
#include "output_report.h"

using namespace std;

namespace hpcscan {

//-------------------------------------------------------------------------------------------------------

Hardware_DPCPP::Hardware_DPCPP(string gridMode) : Hardware(gridMode)
		{
	printDebug(LIGHT_DEBUG, "IN Hardware_DPCPP::Hardware_DPCPP");

	supportGetPowerUsage = checkSupportGetPowerUsage() ;

	printDebug(LIGHT_DEBUG, "OUT Hardware_DPCPP::Hardware_DPCPP");
		}

//-------------------------------------------------------------------------------------------------------

Hardware_DPCPP::~Hardware_DPCPP(void)
{
	printDebug(LIGHT_DEBUG, "IN Hardware_DPCPP::~Hardware_DPCPP");

	// TODO ~Hardware_DPCPP

	printDebug(LIGHT_DEBUG, "OUT Hardware_DPCPP::~Hardware_DPCPP");
}

//-------------------------------------------------------------------------------------------------------

void Hardware_DPCPP::info(void)
{
	printDebug(LIGHT_DEBUG, "IN Hardware_DPCPP::info");

	print_blank() ;
	printInfo(MASTER, " Hardware information") ;
	printInfo(MASTER, " Access devices via DPC++") ;

	// display all available devices
	// TODO display all available devices
	{
		// initialize device queue
		sycl::queue *myQ ;
		try {
			if (Config::Instance()->dpcppSelect.compare("Host") == 0)
			{
				myQ = new sycl::queue( sycl:: host_selector{} ) ;
			}
			else if (Config::Instance()->dpcppSelect.compare("CPU") == 0)
			{
				myQ = new sycl::queue( sycl:: cpu_selector{} ) ;
			}
			else if (Config::Instance()->dpcppSelect.compare("GPU") == 0)
			{
				myQ = new sycl::queue( sycl:: gpu_selector{} ) ;
			}
			else
			{
				printError("In Hardware_DPCPP::info, invalid dpcpp selector") ;
				return ;
			}
		} catch (exception const& ex) {
			printError("In Hardware_DPCPP::info, requested device is not available") ;
			return ;
		}
		printInfo(MASTER, " Selected device", myQ->get_device().get_info<sycl::info::device::name>() ) ;
		printInfo(MASTER, " Device vendor\t", myQ->get_device().get_info<sycl::info::device::vendor>() ) ;
		delete(myQ) ;
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

	printDebug(LIGHT_DEBUG, "OUT Hardware_DPCPP::info");
}

//-------------------------------------------------------------------------------------------------------

bool Hardware_DPCPP::checkSupportGetPowerUsage(void)
{
	printDebug(LIGHT_DEBUG, "IN Hardware_DPCPP::checkSupportGetPowerUsage");

	bool retVal = false ;

	// get power consumption
	{
		// TODO get power consumption
	}

	printDebug(LIGHT_DEBUG, "OUT Hardware_DPCPP::checkSupportGetPowerUsage");
	return(retVal) ;
}

//-------------------------------------------------------------------------------------------------------

Myfloat Hardware_DPCPP::measureCurrentPower(void)
{
	printDebug(LIGHT_DEBUG, "IN Hardware_DPCPP::measureCurrentPower");

	Myfloat retVal = UNSPECIFIED ;
	if (supportGetPowerUsage)
	{
		// TODO measureCurrentPower
	}

	printDebug(LIGHT_DEBUG, "OUT Hardware_DPCPP::measureCurrentPower");
	return(retVal) ;
}


} // namespace hpcscan
