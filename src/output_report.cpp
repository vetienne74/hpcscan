//-------------------------------------------------------------------------------------------------------
// Handle textual output to screen and log files with functions:
// printDebug    = for debugging purpose (output in log files, one per MPI process)
// printError    = for internal errors (output cout)
// printInfo     = for displaying general text output (output cout)
// printWarning  = for warning messages (output cout)
//-------------------------------------------------------------------------------------------------------

#include "output_report.h"

#include <cassert>
#include <fstream>
#include <iostream>
#include <string>

#include "config.h"
#include "global.h"
#include "type_def.h"
#include "version_hpcscan.h"

using namespace std;

namespace hpcscan {

static const string CURRENT_VERSION = "1.1" ;

static const string LINE_REPORT_T1 = "================================================================================" ;
static const string LINE_REPORT_T2 = "--------------------------------------------------------------------------------" ;
static const string LINE_REPORT_T3 = "................................................................................" ;
static const string LINE_REPORT_T4 = "XxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXx" ;
static const string LINE_REPORT_T5 = "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*" ;

//-------------------------------------------------------------------------------------------------------

Rtn_code print_header_of_output_report(void)
{
	if (myid_world == 0)
	{
		print_line1() ;
		cout << "\t\t\t H P C S C A N - ver " << CURRENT_VERSION << "\n\n" ;
		printInfo(MASTER, " Git version", HPCSCAN_GIT_COMMIT) ;
		printInfo(MASTER, "\t", HPCSCAN_GIT_AUTHOR) ;
		printInfo(MASTER, "\t", HPCSCAN_GIT_DATE) ;
		printInfo(MASTER, "") ;

		printInfo(MASTER, " Host name", Config::Instance()->hostName) ;
		printInfo(MASTER, " User name", Config::Instance()->userName) ;

		// single or double precision
#ifdef _DOUBLE_PRECISION_
		printInfo(MASTER, " Computation", "DOUBLE PRECISION") ;
#else
		printInfo(MASTER, " Computation", "SINGLE PRECISION") ;
#endif

		printInfo(MASTER, " Mode Baseline", "ENABLED") ;
		printInfo(MASTER, " Mode CacheBlk", "ENABLED") ;

#ifdef __CUDA__
		printInfo(MASTER, " Mode CUDA", "ENABLED") ;
#else
		printInfo(MASTER, " Mode CUDA", "DISABLED") ;
#endif

#ifdef __HIP__
		printInfo(MASTER, " Mode HIP", "ENABLED") ;
#else
		printInfo(MASTER, " Mode HIP", "DISABLED") ;
#endif

#ifdef __OPENACC__
		printInfo(MASTER, " Mode OpenAcc", "ENABLED") ;
#else
		printInfo(MASTER, " Mode OpenAcc", "DISABLED") ;
#endif

#ifdef __NEC__
		printInfo(MASTER, " Mode NEC", "ENABLED") ;
		printInfo(MASTER, " Mode NEC_SCA", "ENABLED") ;
#else
		printInfo(MASTER, " Mode NEC", "DISABLED") ;
		printInfo(MASTER, " Mode NEC_SCA", "DISABLED") ;
#endif
		print_line1() ;
	}

	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code print_end_of_output_report(Rtn_code rtnCode)
{
	if (myid_world == 0)
	{
		cout << "\n" ;
		if (rtnCode == RTN_CODE_OK)
		{
			print_line1() ;
			cout << "\t\t\t HPCSCAN TERMINATED SUCCESSFULLY\n" ;
			print_line1() ;
		}
		else if (rtnCode == RTN_CODE_KO)
		{
			print_line1() ;
			cout << "\t\t\t HPCSCAN TERMINATED WITH ERROR(S)\n" ;
			print_line1() ;
		}
	}

	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

void printDebug(Debug_level debug_l, const char* text)
{
	if (debug_l <= debug)
	{
		Config::Instance()->debugLogFile << text << "\n" << flush ;
	}
}

//-------------------------------------------------------------------------------------------------------

void printDebug(Debug_level debug_l, const char* text, string text2)
{
	if (debug_l <= debug)
	{
		Config::Instance()->debugLogFile << text << " " << text2 << "\n" << flush ;
	}
}

//-------------------------------------------------------------------------------------------------------

void printDebug(Debug_level debug_l, const char* text, Myint64 val)
{
	if (debug_l <= debug)
	{
		Config::Instance()->debugLogFile << text << " " << val << "\n" << flush ;
	}
}

//-------------------------------------------------------------------------------------------------------

void printDebug(Debug_level debug_l, const char* text, Myint nb)
{
	if (debug_l <= debug)
	{
		Config::Instance()->debugLogFile << text << " " << nb << "\n" << flush ;
	}
}

//-------------------------------------------------------------------------------------------------------

void printDebug(Debug_level debug_l, const char* text, Myint nb, Myint nb2)
{
	if (debug_l <= debug)
	{
		Config::Instance()->debugLogFile << text << " " << nb << " " << nb2 << "\n" << flush ;
	}
}

//-------------------------------------------------------------------------------------------------------

void printDebug(Debug_level debug_l, const char* text, Myfloat32 nb)
{
	if (debug_l <= debug)
	{
		Config::Instance()->debugLogFile << text << " " << nb << "\n" << flush ;
	}
}

//-------------------------------------------------------------------------------------------------------

void printDebug(Debug_level debug_l, const char* text, Myfloat64 nb)
{
	if (debug_l <= debug)
	{
		Config::Instance()->debugLogFile << text << " " << nb << "\n" << flush ;
	}
}

//-------------------------------------------------------------------------------------------------------

void printDebug(Debug_level debug_l, const char* text, const char* text2)
{
	if (debug_l <= debug)
	{
		Config::Instance()->debugLogFile << text << "\t" << text2 << "\n" << flush ;
	}
}

//-------------------------------------------------------------------------------------------------------

void print_line1(void)
{
	if (myid_world == 0) cout << LINE_REPORT_T1 << "\n" << flush ;
}
void print_line2(void)
{
	if (myid_world == 0) cout << LINE_REPORT_T2 << "\n" << flush ;
}
void print_line3(void)
{
	if (myid_world == 0) cout << LINE_REPORT_T3 << "\n" << flush ;
}
void print_line4(void)
{
	if (myid_world == 0) cout << LINE_REPORT_T4 << "\n" << flush ;
}
void print_line5(void)
{
	if (myid_world == 0) cout << LINE_REPORT_T5 << "\n" << flush ;
}
void print_blank(void)
{
	if (myid_world == 0) cout << "\n" << flush ;
}

//-------------------------------------------------------------------------------------------------------

void printInfo(Display_type display_t, const char* text)
{
	if ((display_t == ALL) || ((display_t == MASTER) && (myid_world == 0)))
	{
		cout << text << "\n" << flush ;
	}
}

//-------------------------------------------------------------------------------------------------------

void printInfo(Display_type display_t, const char* text, Myint nb)
{
	if ((display_t == ALL) || ((display_t == MASTER) && (myid_world == 0)))
	{
		cout << text << "\t" << nb << "\n" << flush ;
	}
}

//-------------------------------------------------------------------------------------------------------

void printInfo(Display_type display_t, const char* text, Myint nb, Myint nb2)
{
	if ((display_t == ALL) || ((display_t == MASTER) && (myid_world == 0)))
	{
		cout << text << "\t" << nb << " " << nb2 << "\n" << flush ;
	}
}

//-------------------------------------------------------------------------------------------------------

void printInfo(Display_type display_t, const char* text, Myint64 nb)
{
	if ((display_t == ALL) || ((display_t == MASTER) && (myid_world == 0)))
	{
		cout << text << "\t" << nb << "\n" << flush ;
	}
}

//-------------------------------------------------------------------------------------------------------

void printInfo(Display_type display_t, const char* text, Myfloat32 nb)
{
	if ((display_t == ALL) || ((display_t == MASTER) && (myid_world == 0)))
	{
		cout << text << "\t" << nb << "\n" << flush ;
	}
}

//-------------------------------------------------------------------------------------------------------

void printInfo(Display_type display_t, const char* text, Myfloat64 nb)
{
	if ((display_t == ALL) || ((display_t == MASTER) && (myid_world == 0)))
	{
		cout << text << "\t" << nb << "\n" << flush ;
	}
}

//-------------------------------------------------------------------------------------------------------

void printInfo(Display_type display_t, const char* text, string text2)
{
	if ((display_t == ALL) || ((display_t == MASTER) && (myid_world == 0)))
	{
		cout << text << "\t" << text2 << "\n" << flush ;
	}
}

//-------------------------------------------------------------------------------------------------------

void printInfo(Display_type display_t, const char* text, const char* text2)
{
	if ((display_t == ALL) || ((display_t == MASTER) && (myid_world == 0)))
	{
		cout << text << "\t" << text2 << "\n" << flush ;
	}
}

//-------------------------------------------------------------------------------------------------------

void printInfo(Display_type display_t, string text, const char* text2)
{
	if ((display_t == ALL) || ((display_t == MASTER) && (myid_world == 0)))
	{
		cout << text << "\t" << text2 << "\n" << flush ;
	}
}

//-------------------------------------------------------------------------------------------------------

void printWarning(const char* text, Myint nb)
{
	cout << "\n" ;
	cout << LINE_REPORT_T3 ;
	cout << "\n\t\t\t W A R N I N G \n" ;
	cout << text << "\t" << nb << "\n\n" ;
	cout << LINE_REPORT_T3 ;
	cout << "\n" ;
}

//-------------------------------------------------------------------------------------------------------

void printWarning(const char* text)
{
	cout << "\n" ;
	cout << LINE_REPORT_T3 ;
	cout << "\n\t\t\t W A R N I N G \n" ;
	cout << text << "\n\n" ;
	cout << LINE_REPORT_T3 ;
	cout << "\n" ;
}

//-------------------------------------------------------------------------------------------------------

void printError(string* text)
{
	cout << "\n" ;
	cout << LINE_REPORT_T4 ;
	cout << "\n\t\t\t F A T A L    E R R O R\n\n" ;
	cout << *text << "\n\n" ;
	cout << LINE_REPORT_T4 ;
	cout << "\n" ;
}

//-------------------------------------------------------------------------------------------------------

void printError(const char* text, Myint nb)
{
	cout << "\n" ;
	cout << LINE_REPORT_T4 ;
	cout << "\n\t\t\t F A T A L    E R R O R\n\n" ;
	cout << text << "\t" << nb << "\n\n" ;
	cout << LINE_REPORT_T4 ;
	cout << "\n" ;
}

//-------------------------------------------------------------------------------------------------------

void printError(const char* text, string text2)
{
	cout << "\n" ;
	cout << LINE_REPORT_T4 ;
	cout << "\n\t\t\t F A T A L    E R R O R\n\n" ;
	cout << text << "\t" << text2 << "\n\n" ;
	cout << LINE_REPORT_T4 ;
	cout << "\n" ;
}

//-------------------------------------------------------------------------------------------------------

void printError(const char* text, Myint nb1, Myint nb2)
{
	cout << "\n" ;
	cout << LINE_REPORT_T4 ;
	cout << "\n\t\t\t F A T A L    E R R O R\n\n" ;
	cout << text << "\t" << nb1 << "\t" << nb2 << "\n\n" ;
	cout << LINE_REPORT_T4 ;
	cout << "\n" ;
}

//-------------------------------------------------------------------------------------------------------

void printError(const char* text)
{
	cout << "\n" ;
	cout << LINE_REPORT_T4 ;
	cout << "\n\t\t\t F A T A L    E R R O R\n\n" ;
	cout << text << "\n\n" ;
	cout << LINE_REPORT_T4 ;
	cout << "\n" ;
}

//-------------------------------------------------------------------------------------------------------

void printError(const char* text1, const char* text2)
{
	cout << "\n" ;
	cout << LINE_REPORT_T4 ;
	cout << "\n\t\t\t F A T A L    E R R O R\n\n" ;
	cout << text1 << "\t" << text2 << "\n\n" ;
	cout << LINE_REPORT_T4 ;
	cout << "\n" ;
}

} // namespace hpcscan
