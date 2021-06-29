
//-------------------------------------------------------------------------------------------------------
// This class handles all queries related to the hardware
// Generic class that targets CPUs
// Associated to the test modes Baseline & CacheBlk
//-------------------------------------------------------------------------------------------------------

#include "hardware.h"

#include "config.h"
#include "constant.h"
#include "global.h"
#include "output_report.h"

using namespace std;

namespace hpcscan {

//-------------------------------------------------------------------------------------------------------

Hardware::Hardware(string gridMode)
{
	printDebug(MID_DEBUG, "IN Hardware::Hardware");

	// TODO

	printDebug(MID_DEBUG, "OUT Hardware::Hardware");
}

//-------------------------------------------------------------------------------------------------------

Hardware::~Hardware(void)
{
	printDebug(MID_DEBUG, "IN Hardware::~Hardware");

	// TODO

	printDebug(MID_DEBUG, "OUT Hardware::~Hardware");
}

//-------------------------------------------------------------------------------------------------------

void Hardware::info(void)
{
	printDebug(MID_DEBUG, "IN Hardware::info");

	printInfo(MASTER, " Hardware information") ;
	printInfo(MASTER, " Generic CPU") ;
	printInfo(MASTER, " No information available") ;

	if (supportGetPowerUsage())
	{
		printInfo(MASTER, " Read power usage", "SUPPORTED") ;
	}
	else
	{
		printInfo(MASTER, " Read power usage", "NOT SUPPORTED") ;
	}

	print_line5() ;

	printDebug(MID_DEBUG, "OUT Hardware::info");
}

//-------------------------------------------------------------------------------------------------------

bool Hardware::supportGetPowerUsage(void)
{
	printDebug(MID_DEBUG, "IN Hardware::supportGetPowerUsage");

	bool retVal = false ;

	printDebug(MID_DEBUG, "OUT Hardware::supportGetPowerUsage");
	return(retVal) ;
}


} // namespace hpcscan
