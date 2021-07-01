
//-------------------------------------------------------------------------------------------------------
// This class handles all queries related to the hardware
// Generic class that targets CPUs
// Associated to the test modes Baseline & CacheBlk
//-------------------------------------------------------------------------------------------------------

#include "hardware.h"

#include <cfloat>  // for FLT_MAX ;

#include "mpi.h"

#include "config.h"
#include "constant.h"
#include "global.h"
#include "output_report.h"

using namespace std;

namespace hpcscan {

//-------------------------------------------------------------------------------------------------------

Hardware::Hardware(string gridMode)
								{
	printDebug(LIGHT_DEBUG, "IN Hardware::Hardware");

	supportGetPowerUsage = checkSupportGetPowerUsage() ;

	printDebug(LIGHT_DEBUG, "OUT Hardware::Hardware");
								}

//-------------------------------------------------------------------------------------------------------

Hardware::~Hardware(void)
{
	printDebug(LIGHT_DEBUG, "IN Hardware::~Hardware");

	// TODO

	printDebug(LIGHT_DEBUG, "OUT Hardware::~Hardware");
}

//-------------------------------------------------------------------------------------------------------

void Hardware::info(void)
{
	printDebug(LIGHT_DEBUG, "IN Hardware::info");

	printInfo(MASTER, " Hardware information") ;
	printInfo(MASTER, " Generic CPU") ;
	printInfo(MASTER, " No information available") ;

	if (supportGetPowerUsage)
	{
		printInfo(MASTER, " Read power usage", "SUPPORTED") ;
	}
	else
	{
		printInfo(MASTER, " Read power usage", "NOT SUPPORTED") ;
	}

	print_line5() ;

	printDebug(LIGHT_DEBUG, "OUT Hardware::info");
}

//-------------------------------------------------------------------------------------------------------

bool Hardware::checkSupportGetPowerUsage(void)
{
	printDebug(LIGHT_DEBUG, "IN Hardware::checkSupportGetPowerUsage");

	bool retVal = false ;

	printDebug(LIGHT_DEBUG, "OUT Hardware::checkSupportGetPowerUsage");
	return(retVal) ;
}

//-------------------------------------------------------------------------------------------------------

void Hardware::updateHwCounter(void)
{
	printDebug(LIGHT_DEBUG, "IN Hardware::updateHwCounter");

	hwCounter_struct_type hwCounterTmp ;

	// get time
	hwCounterTmp.timeOfMeasure = MPI_Wtime() ;

	// get power usage
	hwCounterTmp.powerWatt = measureCurrentPower() ;

	// add entry
	hwCounterVec.push_back(hwCounterTmp);

	printDebug(LIGHT_DEBUG, "OUT Hardware::updateHwCounter");
}

//-------------------------------------------------------------------------------------------------------

void Hardware::watchTimeAndUpdateHwCounter(void)
{
	printDebug(LIGHT_DEBUG, "IN Hardware::watchTimeAndUpdateHwCounter");

	if (Config::Instance()->hwCounterDt > 0.)
	{
		// no entry found in hwCounterVec
		if (hwCounterVec.size() == 0)
		{
			printWarning("IN Hardware::watchTimeAndUpdateHwCounter, hwCounterVec.size() == 0") ;
		}
		else
		{
			// check it is time to update hardware counters
			auto lastTimeUpdate =hwCounterVec.back().timeOfMeasure ;
			if ((MPI_Wtime() - lastTimeUpdate) >= Config::Instance()->hwCounterDt)
			{
				updateHwCounter() ;
			}
		}
	}

	printDebug(LIGHT_DEBUG, "OUT Hardware::watchTimeAndUpdateHwCounter");
}

//-------------------------------------------------------------------------------------------------------

Myfloat Hardware::measureCurrentPower(void)
{
	printDebug(LIGHT_DEBUG, "IN Hardware::measureCurrentPower");

	Myfloat retVal = UNSPECIFIED ;

	printDebug(LIGHT_DEBUG, "OUT Hardware::measureCurrentPower");
	return(retVal) ;
}

//-------------------------------------------------------------------------------------------------------

void Hardware::displayCounterStat(void)
{
	printDebug(LIGHT_DEBUG, "IN Hardware::displayCounterStat");

	print_blank() ;
	print_line5() ;
	printInfo(MASTER, " Hardware counters statistics") ;

	if (Config::Instance()->hw->supportGetPowerUsage)
	{

		Myint nMeasure = hwCounterVec.size() ;
		printInfo(MASTER, " Number of measures", nMeasure) ;

		if (nMeasure >= 1)
		{

			double tFirst = hwCounterVec[0].timeOfMeasure ;
			double tLast = hwCounterVec[nMeasure-1].timeOfMeasure ;
			double totalTime = tLast - tFirst ;
			printInfo(MASTER, " Stat. over period (s)", totalTime) ;

			// min. power
			{
				Myfloat val = +FLT_MAX ;
				for (Myint ii = 0; ii < hwCounterVec.size(); ii++)
				{
					if (hwCounterVec[ii].powerWatt < val) val = hwCounterVec[ii].powerWatt ;
				}
				printInfo(MASTER, " Min. power (Watt)", val) ;
			}

			// max. power
			{
				Myfloat val = -FLT_MAX ;
				for (Myint ii = 0; ii < hwCounterVec.size(); ii++)
				{
					if (hwCounterVec[ii].powerWatt > val) val = hwCounterVec[ii].powerWatt ;
				}
				printInfo(MASTER, " Max. power (Watt)", val) ;
			}

			// aver. power and consumption
			{
				Myfloat val = 0 ;
				for (Myint ii = 0; ii < hwCounterVec.size(); ii++)
				{
					val += hwCounterVec[ii].powerWatt ;
				}
				val /= hwCounterVec.size() ;
				printInfo(MASTER, " Aver. power (Watt)", val) ;

				printInfo(MASTER, " Aver. consump. (W.h)", val *  totalTime / 3600.) ;
			}
			
		}
	}
	else
	{
		printInfo(MASTER, " No statistic available") ;
	}

	print_line5() ;

	printDebug(LIGHT_DEBUG, "OUT Hardware::displayCounterStat");
}

} // namespace hpcscan
