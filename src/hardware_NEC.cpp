
//-------------------------------------------------------------------------------------------------------
// This class handles all queries related to the hardware
// Specialized class that targets NEC SX-Aurora
// Derived class from Hardware
// Associated to the test modes NEC & NEC_SCA
//-------------------------------------------------------------------------------------------------------

#include "hardware_NEC.h"

#include "config.h"
#include "constant.h"
#include "global.h"
#include "output_report.h"

using namespace std;

namespace hpcscan {

//-------------------------------------------------------------------------------------------------------

Hardware_NEC::Hardware_NEC(string gridMode) : Hardware(gridMode)
{
	printDebug(LIGHT_DEBUG, "IN Hardware_NEC::Hardware_NEC");

	supportGetPowerUsage = checkSupportGetPowerUsage() ;

	printDebug(LIGHT_DEBUG, "OUT Hardware_NEC::Hardware_NEC");
}

//-------------------------------------------------------------------------------------------------------

Hardware_NEC::~Hardware_NEC(void)
{
	printDebug(LIGHT_DEBUG, "IN Hardware_NEC::~Hardware_NEC");

	// TODO ~Hardware_NEC

	printDebug(LIGHT_DEBUG, "OUT Hardware_NEC::~Hardware_NEC");
}

//-------------------------------------------------------------------------------------------------------

void Hardware_NEC::info(void)
{
	printDebug(LIGHT_DEBUG, "IN Hardware_NEC::info");

	print_blank() ;
	printInfo(MASTER, " Hardware information") ;

	printInfo(MASTER, " Target hardware:") ;
	// TODO display target info

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

	printDebug(LIGHT_DEBUG, "OUT Hardware_NEC::info");
}

//-------------------------------------------------------------------------------------------------------

bool Hardware_NEC::checkSupportGetPowerUsage(void)
{
	printDebug(LIGHT_DEBUG, "IN Hardware_NEC::checkSupportGetPowerUsage");

	bool retVal = false ;

	// get power consumption
	{
		// TODO get power consumption
	}

	printDebug(LIGHT_DEBUG, "OUT Hardware_NEC::checkSupportGetPowerUsage");
	return(retVal) ;
}

//-------------------------------------------------------------------------------------------------------

Myfloat Hardware_NEC::measureCurrentPower(void)
{
	printDebug(LIGHT_DEBUG, "IN Hardware_NEC::measureCurrentPower");

	Myfloat retVal = UNSPECIFIED ;
	if (supportGetPowerUsage)
	{
		// TODO measureCurrentPower
	}

	printDebug(LIGHT_DEBUG, "OUT Hardware_NEC::measureCurrentPower");
	return(retVal) ;
}


} // namespace hpcscan
