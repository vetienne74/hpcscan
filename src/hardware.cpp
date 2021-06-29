
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


} // namespace hpcscan
