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
		cout << "\t\t\t H P C S C A N - ver " << CURRENT_VERSION << " (2020)\n\n" ;
		cout << " GIT " << HPCSCAN_GIT_COMMIT << "\n" ;
		cout << " " << HPCSCAN_GIT_AUTHOR << "\n" ;
		cout << " " << HPCSCAN_GIT_DATE << "\n" ;

		// single or double precision
#ifdef _DOUBLE_PRECISION_
		cout << " Computations in DOUBLE PRECISION\n" ;
#else
		cout << " Computations in SINGLE PRECISION\n" ;
#endif

#ifdef __CUDA__
		cout << " CUDA enabled\n" ;
#endif
#ifdef __NEC__
		cout << " NEC enabled\n" ;
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
		print_line1() ;
		if (rtnCode == RTN_CODE_OK) cout << "\t\t\t HPCSCAN TERMINATED SUCCESSFULLY\n" ;
		else cout << "\t\t\t HPCSCAN TERMINATED WITH ERROR(S)\n" ;
		print_line1() ;
	}

	return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

void printDebug(Debug_level debug_l, char* text)
{
	if (debug_l <= debug)
	{
		Config::Instance()->debugLogFile << text << "\n" << flush ;
	}
}

//-------------------------------------------------------------------------------------------------------

void printDebug(Debug_level debug_l, char* text, string text2)
{
	if (debug_l <= debug)
	{
		Config::Instance()->debugLogFile << text << " " << text2 << "\n" << flush ;
	}
}

//-------------------------------------------------------------------------------------------------------

void printDebug(Debug_level debug_l, char* text, Myint64 val)
{
	if (debug_l <= debug)
	{
		Config::Instance()->debugLogFile << text << " " << val << "\n" << flush ;
	}
}

//-------------------------------------------------------------------------------------------------------

void printDebug(Debug_level debug_l, char* text, Myint nb)
{
	if (debug_l <= debug)
	{
		Config::Instance()->debugLogFile << text << " " << nb << "\n" << flush ;
	}
}

//-------------------------------------------------------------------------------------------------------

void printDebug(Debug_level debug_l, char* text, Myint nb, Myint nb2)
{
	if (debug_l <= debug)
	{
		Config::Instance()->debugLogFile << text << " " << nb << " " << nb2 << "\n" << flush ;
	}
}

//-------------------------------------------------------------------------------------------------------

void printDebug(Debug_level debug_l, char* text, Myfloat32 nb)
{
	if (debug_l <= debug)
	{
		Config::Instance()->debugLogFile << text << " " << nb << "\n" << flush ;
	}
}

//-------------------------------------------------------------------------------------------------------

void printDebug(Debug_level debug_l, char* text, Myfloat64 nb)
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

void printInfo(Display_type display_t, char* text)
{
	if ((display_t == ALL) || ((display_t == MASTER) && (myid_world == 0)))
	{
		cout << text << "\n" << flush ;
	}
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

void printInfo(Display_type display_t, char* text, Myint nb)
{
	if ((display_t == ALL) || ((display_t == MASTER) && (myid_world == 0)))
	{
		cout << text << "\t" << nb << "\n" << flush ;
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

void printInfo(Display_type display_t, char* text, Myfloat nb)
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

void printInfo(Display_type display_t, char* text, char* text2)
{
	if ((display_t == ALL) || ((display_t == MASTER) && (myid_world == 0)))
	{
		cout << text << "\t" << text2 << "\n" << flush ;
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

void printInfo(Display_type display_t, const char* text, char* text2)
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

void printWarning(char* text, Myint nb)
{
	cout << "\n" ;
	cout << LINE_REPORT_T3 ;
	cout << "\n\t\t\t W A R N I N G \n" ;
	cout << text << "\t" << nb << "\n\n" ;
	cout << LINE_REPORT_T3 ;
	cout << "\n" ;
}

//-------------------------------------------------------------------------------------------------------

void printWarning(char* text)
{
	cout << "\n" ;
	cout << LINE_REPORT_T3 ;
	cout << "\n\t\t\t W A R N I N G \n" ;
	cout << text << "\n\n" ;
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

void printError(char* text, Myint nb)
{
	cout << "\n" ;
	cout << LINE_REPORT_T4 ;
	cout << "\n\t\t\t F A T A L    E R R O R\n\n" ;
	cout << text << "\t" << nb << "\n\n" ;
	cout << LINE_REPORT_T4 ;
	cout << "\n" ;
}

//-------------------------------------------------------------------------------------------------------

void printError(char* text, string text2)
{
	cout << "\n" ;
	cout << LINE_REPORT_T4 ;
	cout << "\n\t\t\t F A T A L    E R R O R\n\n" ;
	cout << text << "\t" << text2 << "\n\n" ;
	cout << LINE_REPORT_T4 ;
	cout << "\n" ;
}

//-------------------------------------------------------------------------------------------------------

void printError(char* text, Myint nb1, Myint nb2)
{
	cout << "\n" ;
	cout << LINE_REPORT_T4 ;
	cout << "\n\t\t\t F A T A L    E R R O R\n\n" ;
	cout << text << "\t" << nb1 << "\t" << nb2 << "\n\n" ;
	cout << LINE_REPORT_T4 ;
	cout << "\n" ;
}

//-------------------------------------------------------------------------------------------------------

void printError(char* text)
{
	cout << "\n" ;
	cout << LINE_REPORT_T4 ;
	cout << "\n\t\t\t F A T A L    E R R O R\n\n" ;
	cout << text << "\n\n" ;
	cout << LINE_REPORT_T4 ;
	cout << "\n" ;
}

//-------------------------------------------------------------------------------------------------------

void printError(char* text1, char* text2)
{
	cout << "\n" ;
	cout << LINE_REPORT_T4 ;
	cout << "\n\t\t\t F A T A L    E R R O R\n\n" ;
	cout << text1 << "\t" << text2 << "\n\n" ;
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

} // namespace hpcscan
