
//-------------------------------------------------------------------------------------------------------
// This class handles all queries related to the hardware
// Generic class that targets CPUs
// Associated to the test modes Baseline & CacheBlk
//-------------------------------------------------------------------------------------------------------

#include "hardware.h"

#include <cfloat>  // for FLT_MAX ;
#ifndef __NEC__
#include <cpuid.h>
#endif
#include <cstring> // needed for DPC++ (memset, memcpy)

#include "mpi.h"
#include <omp.h>

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

	// TODO ~Hardware

	printDebug(LIGHT_DEBUG, "OUT Hardware::~Hardware");
}

//-------------------------------------------------------------------------------------------------------

void Hardware::info(void)
{
	printDebug(LIGHT_DEBUG, "IN Hardware::info");

	print_blank();
	printInfo(MASTER, " Hardware information") ;
	printInfo(MASTER, " Target hardware:") ;
	hostInfo() ;

	if (supportGetPowerUsage)
	{
		printInfo(MASTER, " Read power usage", "SUPPORTED") ;
	}
	else
	{
		printInfo(MASTER, " Read power usage", "NOT SUPPORTED") ;
	}
	print_line2() ;

	printDebug(LIGHT_DEBUG, "OUT Hardware::info");
}

//-------------------------------------------------------------------------------------------------------

void Hardware::hostInfo(void)
{
	printDebug(LIGHT_DEBUG, "IN Hardware::hostInfo");
#ifndef __NEC__
	{
		// from https://stackoverflow.com/questions/850774/how-to-determine-the-hardware-cpu-and-ram-on-a-machine
		char CPUBrandString[0x40];
		unsigned int CPUInfo[4] = {0,0,0,0};

		__cpuid(0x80000000, CPUInfo[0], CPUInfo[1], CPUInfo[2], CPUInfo[3]);
		unsigned int nExIds = CPUInfo[0];

		memset(CPUBrandString, 0, sizeof(CPUBrandString));

		for (unsigned int i = 0x80000000; i <= nExIds; ++i)
		{
			__cpuid(i, CPUInfo[0], CPUInfo[1], CPUInfo[2], CPUInfo[3]);

			if (i == 0x80000002)
				memcpy(CPUBrandString, CPUInfo, sizeof(CPUInfo));
			else if (i == 0x80000003)
				memcpy(CPUBrandString + 16, CPUInfo, sizeof(CPUInfo));
			else if (i == 0x80000004)
				memcpy(CPUBrandString + 32, CPUInfo, sizeof(CPUInfo));
		}

		printInfo(MASTER, " CPU (Host) type", CPUBrandString) ;
	}
#endif
	printDebug(LIGHT_DEBUG, "OUT Hardware::hostInfo");
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

	if (supportGetPowerUsage)
	{

		Myint nMeasure = hwCounterVec.size() ;
		printInfo(MASTER, " Number of measures", nMeasure) ;

		if (nMeasure >= 1)
		{

			double tFirst = hwCounterVec[0].timeOfMeasure ;
			double tLast = hwCounterVec[nMeasure-1].timeOfMeasure ;
			double statTime = tLast - tFirst ;
			printInfo(MASTER, " Stat. over period (s)", statTime) ;

			// min. power
			Myfloat minPower = +FLT_MAX ;
			{

				for (Myint ii = 0; ii < hwCounterVec.size(); ii++)
				{
					if (hwCounterVec[ii].powerWatt < minPower) minPower = hwCounterVec[ii].powerWatt ;
				}
				printInfo(MASTER, " Min. power (Watt)", minPower) ;
			}

			// max. power
			Myfloat maxPower = -FLT_MAX ;
			{
				for (Myint ii = 0; ii < hwCounterVec.size(); ii++)
				{
					if (hwCounterVec[ii].powerWatt > maxPower) maxPower = hwCounterVec[ii].powerWatt ;
				}
				printInfo(MASTER, " Max. power (Watt)", maxPower) ;
			}

			// aver. power and consumption
			Myfloat averPower = 0. ;
			Myfloat wattHour  = 0. ;
			{
				for (Myint ii = 0; ii < hwCounterVec.size(); ii++)
				{
					averPower += hwCounterVec[ii].powerWatt ;
				}
				averPower /= hwCounterVec.size() ;
				printInfo(MASTER, " Aver. power (Watt)", averPower) ;
				wattHour = averPower *  statTime / 3600. ;
				printInfo(MASTER, " Aver. consump. (W.h)", wattHour) ;
			}

			// write measurements in log file
			if (myMpiRank == 0)
			{
				string file_name = "hpcscan.hwCounter." + Config::Instance()->testCaseName + ".log";
				ofstream hwCounterLogFile ;
				hwCounterLogFile.open(file_name, ios::app) ;

				// all strings first
				hwCounterLogFile << Config::Instance()->hostName << " " ;
				hwCounterLogFile << Config::Instance()->testCaseName << " " ;
				hwCounterLogFile << Config::Instance()->testMode << " " ;
				hwCounterLogFile << Config::Instance()->propagator << " " ;

				// numeric values follow
				hwCounterLogFile << nMpiProc << " " ; // 1
				hwCounterLogFile << Config::Instance()->nsub1 << " " ; // 2
				hwCounterLogFile << Config::Instance()->nsub2 << " " ; // 3
				hwCounterLogFile << Config::Instance()->nsub3 << " " ; // 4
				hwCounterLogFile << omp_get_max_threads() << " " ; // 5
				hwCounterLogFile << Config::Instance()->n1 << " " ; // 6
				hwCounterLogFile << Config::Instance()->n2 << " " ; // 7
				hwCounterLogFile << Config::Instance()->n3 << " " ; // 8
				hwCounterLogFile << Config::Instance()->fdOrder << " " ; // 9

				hwCounterLogFile << statTime << " " ; // 10
				hwCounterLogFile << hwCounterVec.size() << " " ; // 11
				hwCounterLogFile << minPower << " " ; // 12
				hwCounterLogFile << maxPower << " " ; // 13
				hwCounterLogFile << averPower << " " ; // 14
				hwCounterLogFile << wattHour << " " ; // 15
				hwCounterLogFile << "\n" ;

				hwCounterLogFile.close() ;
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
