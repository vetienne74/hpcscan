
//-------------------------------------------------------------------------------------------------------
// Singleton class handling all configuration parameters
// Parameters are set to default values and overwritten by command line arguments
//-------------------------------------------------------------------------------------------------------

#include "config.h"

#include <cstddef> // for NULL
#include <stdio.h>
#include <unistd.h> // for hostname and username

#include <omp.h>

#include "constant.h"
#include "global.h"
#include "grid.h"
#include "output_report.h"

using namespace std;

namespace hpcscan {

//-------------------------------------------------------------------------------------------------------
// Constant definitions
//-------------------------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------------------------
// Default values
//-------------------------------------------------------------------------------------------------------

static const bool      DEFAULT_AUTO_PAD       = false ;
static const string    DEFAULT_BOUNDARY       = "FreeSurf" ;
static const Myint     DEFAULT_CB1            = 9999 ;
static const Myint     DEFAULT_CB2            = 4 ;
static const Myint     DEFAULT_CB3            = 16 ;
static const Dim_type  DEFAULT_DIM            = DIM3 ;
static const Myfloat64 DEFAULT_DT             = 0.0 ; // stable dt will be used
static const Myint     DEFAULT_FD_ORDER       = 8 ;
static const Myint     DEFAULT_GPU_BLKSIZE    = 256 ;
static const Myint     DEFAULT_GPU_GRIDSIZE   = 512 ;
static const Myfloat64 DEFAULT_H              = PI / 30 ;
static const Myint     DEFAULT_INNER_N1       = 61 ;
static const Myint     DEFAULT_INNER_N2       = 61 ;
static const Myint     DEFAULT_INNER_N3       = 61 ;
static const Myint     DEFAULT_N1_ADD_PAD     = UNSPECIFIED ;
static const Myint     DEFAULT_N2_ADD_PAD     = UNSPECIFIED ;
static const Myint     DEFAULT_N3_ADD_PAD     = UNSPECIFIED ;
static const Myint     DEFAULT_N1_MUL_PAD     = UNSPECIFIED ;
static const Myint     DEFAULT_N2_MUL_PAD     = UNSPECIFIED ;
static const Myint     DEFAULT_N3_MUL_PAD     = UNSPECIFIED ;
static const Myint     DEFAULT_NLAYER         = 0 ;
static const Myint     DEFAULT_NSUB1          = 1 ;
static const Myint     DEFAULT_NSUB2          = 1 ;
static const Myint     DEFAULT_NSUB3          = 1 ;
static const Myint     DEFAULT_NT             = 20 ;
static const Myint     DEFAULT_NTRY           = 10 ;
static const Myfloat64 DEFAULT_PARAM1         = 0.8 ;
static const Myfloat64 DEFAULT_PARAM2         = 0.9 ;
static const Myfloat64 DEFAULT_PARAM3         = 1.0 ;
static const Myfloat64 DEFAULT_PARAM4         = 1.1 ;
static const string    DEFAULT_PROPAGATOR     = PROPA_TYPE_AC2STANDARD ;
static const Myfloat64 DEFAULT_RATIO_CFL      = 1.0 ;
static const Myfloat64 DEFAULT_SNAP_DT        = 0.0 ;
static const Myint     DEFAULT_SNAP_INC       = 1 ;
static const string    DEFAULT_TEST_CASE      = "UNSPECIFIED" ;
static const string    DEFAULT_TEST_MODE      = GRID_MODE_BASELINE ;
static const Myfloat64 DEFAULT_TMAX           = UNSPECIFIED ;
static const bool      DEFAULT_WRITE_GRID     = false ;

//-------------------------------------------------------------------------------------------------------
// Min and max
//-------------------------------------------------------------------------------------------------------
static const Myint     MIN_NPOINT_NSUB        = 20 ;

//-------------------------------------------------------------------------------------------------------

// Global static pointer used to ensure a single instance of the class.
Config* Config::m_pInstance = NULL;

//-------------------------------------------------------------------------------------------------------

// This function is called to create an instance of the class.
// Calling the constructor publicly is not allowed. The constructor
// is private and is only called by this Instance function.

Config* Config::Instance()
{
	if (!m_pInstance)   // Only allow one instance of class to be generated.
		m_pInstance = new Config;

	return m_pInstance;
}

//-------------------------------------------------------------------------------------------------------

// constructor

Config::Config(void)
{
	autoPad      = DEFAULT_AUTO_PAD ;
	boundary     = DEFAULT_BOUNDARY ;
	cb1          = DEFAULT_CB1 ;
	cb2          = DEFAULT_CB2 ;
	cb3          = DEFAULT_CB3 ;
	dim          = DEFAULT_DIM ;
	dt           = DEFAULT_DT ;
	fdOrder      = DEFAULT_FD_ORDER ;
	gpuBlkSize   = DEFAULT_GPU_BLKSIZE ;
	gpuGridSize  = DEFAULT_GPU_GRIDSIZE ;
	h            = DEFAULT_H ;
	n1           = DEFAULT_INNER_N1 ;
	n2           = DEFAULT_INNER_N2 ;
	n3           = DEFAULT_INNER_N3 ;
	n1AddPad     = DEFAULT_N1_ADD_PAD ;
	n2AddPad     = DEFAULT_N2_ADD_PAD ;
	n3AddPad     = DEFAULT_N3_ADD_PAD ;
	n1MulPad     = DEFAULT_N1_MUL_PAD ;
	n2MulPad     = DEFAULT_N2_MUL_PAD ;
	n3MulPad     = DEFAULT_N3_MUL_PAD ;
	nlayer       = DEFAULT_NLAYER ;
	nt           = DEFAULT_NT ;
	ntry         = DEFAULT_NTRY ;
	nsub1        = DEFAULT_NSUB1 ;
	nsub2        = DEFAULT_NSUB2 ;
	nsub3        = DEFAULT_NSUB3 ;
	param1       = DEFAULT_PARAM1 ;
	param2       = DEFAULT_PARAM2 ;
	param3       = DEFAULT_PARAM3 ;
	param4       = DEFAULT_PARAM4 ;
	propagator   = DEFAULT_PROPAGATOR ;
	ratioCFL     = DEFAULT_RATIO_CFL ;
	snapDt       = DEFAULT_SNAP_DT ;
	snapInc      = DEFAULT_SNAP_INC ;
	testCaseName = DEFAULT_TEST_CASE ;
	testMode     = DEFAULT_TEST_MODE ;
	tmax         = DEFAULT_TMAX ;
	writeGrid    = DEFAULT_WRITE_GRID ;

	gethostname(hostName, HOST_NAME_MAX);
	getlogin_r(userName, LOGIN_NAME_MAX);
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Config::parse_argument(int argc, char* argv[])
{
	printDebug(MID_DEBUG, "IN Config::parse_argument");

	// check for -version or -help or -h
	Myint ii = 0 ;
	while (ii < argc)
	{
		if ((string(argv[ii]) == "-version") ||  (string(argv[ii]) == "-v"))
		{
			print_header_of_output_report() ;
			// exit
			return(RTN_CODE_EXIT);
		}
		else if ((argc == 1) || (string(argv[ii]) == "-help") || (string(argv[ii]) == "-h"))
		{
			printInfo(MASTER, " List of command line parameters:") ;
			printInfo(MASTER, " -autoPad             = automatic grid padding on all axis") ;
			printInfo(MASTER, "     (if specified autoPad overrides all other padding options)") ;
			printInfo(MASTER, " -boundary <string>   = boundary condition type") ;
			printInfo(MASTER, "     FreeSurf         * Free surface (DEFAULT)") ;
			printInfo(MASTER, "     None             * No boundary condition") ;
			printInfo(MASTER, " -cb1 <int>           = cache block size axis 1 [grid pts]") ;
			printInfo(MASTER, " -cb2 <int>           = cache block size axis 2 [grid pts]") ;
			printInfo(MASTER, " -cb3 <int>           = cache block size axis 3 [grid pts}") ;
			printInfo(MASTER, " -debug <OPT>         = debug trace [OPT=none/light/mid/full]") ;
			printInfo(MASTER, " -dim <int>           = space dimension [1,2 or 3]") ;
			printInfo(MASTER, " -dt <float>          = time step (s) for propagator") ;
			printInfo(MASTER, " -fdOrder <int>       = spatial FD order [2, 4, 8, 12, 16]") ;
			printInfo(MASTER, " -gpuBlkSize <int>    = GPU, number of threads per block") ;
			printInfo(MASTER, " -gpuGridSize <int>   = GPU, number of blocks per grid") ;
			printInfo(MASTER, " -help or -h          = list of command line parameters") ;
			printInfo(MASTER, " -n1 <int>            = inner domain size axis 1 [grid pts]") ;
			printInfo(MASTER, " -n2 <int>            = inner domain size axis 2 [grid pts]") ;
			printInfo(MASTER, " -n3 <int>            = inner domain size axis 3 [grid pts]") ;
			printInfo(MASTER, " -n1AddPad <int>      = add N points to grid along axis 1") ;
			printInfo(MASTER, " -n2AddPad <int>      = add N points to grid along axis 2") ;
			printInfo(MASTER, " -n3AddPad <int>      = add N points to grid along axis 3") ;
			printInfo(MASTER, " -n1MulPad <int>      = grid size multiple of N axis 1") ;
			printInfo(MASTER, " -n2MulPad <int>      = grid size multiple of N axis 2") ;
			printInfo(MASTER, " -n3MulPad <int>      = grid size multiple of N axis 3") ;
			printInfo(MASTER, " -nsub1 <int>         = no. of subdomains axis 1") ;
			printInfo(MASTER, " -nsub2 <int>         = no. of subdomains axis 2") ;
			printInfo(MASTER, " -nsub3 <int>         = no. of subdomains axis 3") ;
			printInfo(MASTER, " -nt <int>            = no. of time steps for propagator") ;
			printInfo(MASTER, " -ntry <int>          = no. of tries for each testCase") ;
			printInfo(MASTER, " -param1 <float>      = parameter 1 used in testCases") ;
			printInfo(MASTER, " -param2 <float>      = parameter 2 used in testCases") ;
			printInfo(MASTER, " -param3 <float>      = parameter 3 used in testCases") ;
			printInfo(MASTER, " -param4 <float>      = parameter 4 used in testCases") ;
			printInfo(MASTER, " -propagator <string> = propagator type") ;
			printInfo(MASTER, "     Ac2Standard      * Acoustic 2nd Standard (DEFAULT)") ;
			printInfo(MASTER, "     Ac2SplitComp     * Acoustic 2nd Separate Laplacian computation") ;
			printInfo(MASTER, " -ratioCFL <float>    = ratio of stability dt for propagator") ;
			printInfo(MASTER, " -snapDt <float>      = snaphots increment (time in sec.)") ;
			printInfo(MASTER, "     (if specified snapDt overrides snapInc)") ;
			printInfo(MASTER, " -snapInc <int>       = snaphots increment (no. of time steps)") ;
			printInfo(MASTER, " -testCase <string>   = run specific testCase by name") ;
			printInfo(MASTER, "     All              * All test cases (DEFAULT)") ;
			printInfo(MASTER, "     Comm             * MPI communication") ;
			printInfo(MASTER, "     FD_D2            * Finite-difference computation") ;
			printInfo(MASTER, "     Grid             * Grid operation") ;
			printInfo(MASTER, "     Memory           * Memory") ;
			printInfo(MASTER, "     Propa            * Propagator") ;
			printInfo(MASTER, "     Util             * Utility for developers") ;
			printInfo(MASTER, " -testMode <string>   = test mode") ;
			printInfo(MASTER, "     Baseline         * CPU without optimization (DEFAULT)") ;
			printInfo(MASTER, "     CacheBlk         * CPU with cache blocking techniques") ;
			printInfo(MASTER, "     CUDA             * NVIDIA GPU with CUDA without optimization") ;
			printInfo(MASTER, "     DPC++            * CPU/GPU/FPGA with DPC++ without optimization") ;
			printInfo(MASTER, "     HIP              * AMD GPU with HIP without optimization") ;
			printInfo(MASTER, "     NEC              * NEC SX-Aurora with compiler directives") ;
			printInfo(MASTER, "     NEC_SCA          * NEC SX-Aurora with Stencil Code Accelerator") ;
//			printInfo(MASTER, "     OpenAcc          * GPU with OpenAcc without optimization") ;
			printInfo(MASTER, " -tmax <float>        = max. time (s) for propagator") ;
			printInfo(MASTER, "     (if specified tmax overrides nt)") ;
			printInfo(MASTER, " -version or -v       = print version information") ;
			printInfo(MASTER, " -writeGrid           = write grids on disk") ;

			// exit
			return(RTN_CODE_EXIT);
		}
		ii++ ;
	}

	print_line2() ;
	printInfo(MASTER, " COMMAND LINE ARGUMENTS\n") ;

	if (argc <= 1)
	{
		printInfo(MASTER, " No input arguments") ;
		return(RTN_CODE_OK);
	}
	else
	{
		printInfo(MASTER, " Executable\t ", argv[0]) ;
		ii = 1 ;

		// loop on all arguments
		while (ii < argc)
		{

			if (string(argv[ii]) == "-autoPad")
			{
				printInfo(MASTER, " Automatic padding","YES") ;
				autoPad = true ;
			}

			else if (string(argv[ii]) == "-boundary")
			{
				ii++ ;
				if (ii >= argc)
				{
					printError(" parameter is needed after -boundary") ;
					return(RTN_CODE_KO) ;
				}
				boundary = argv[ii];
				printInfo(MASTER, " Boundary\t", boundary) ;
			}

			else if (string(argv[ii]) == "-cb1")
			{
				ii++ ;
				if (ii >= argc)
				{
					printError(" parameter is needed after -cb1") ;
					return(RTN_CODE_KO) ;
				}
				cb1 = atoi(argv[ii]);
				printInfo(MASTER, " cb1\t\t", cb1) ;
				if (cb1 <= 0)
				{
					printError(" cb1 should be > 0") ;
					return(RTN_CODE_KO) ;
				}
			}

			else if (string(argv[ii]) == "-cb2")
			{
				ii++ ;
				if (ii >= argc)
				{
					printError(" parameter is needed after -cb2") ;
					return(RTN_CODE_KO) ;
				}
				cb2 = atoi(argv[ii]);
				printInfo(MASTER, " cb2\t\t", cb2) ;
				if (cb2 <= 0)
				{
					printError(" cb2 should be > 0") ;
					return(RTN_CODE_KO) ;
				}
			}

			else if (string(argv[ii]) == "-cb3")
			{
				ii++ ;
				if (ii >= argc)
				{
					printError(" parameter is needed after -cb3") ;
					return(RTN_CODE_KO) ;
				}
				cb3 = atoi(argv[ii]);
				printInfo(MASTER, " cb3\t\t", cb3) ;
				if (cb3 <= 0)
				{
					printError(" cb3 should be > 0") ;
					return(RTN_CODE_KO) ;
				}
			}

			else if (string(argv[ii]) == "-debug")
			{
				ii++ ;
				if (ii >= argc)
				{
					printError(" parameter is needed after -debug") ;
					return(RTN_CODE_KO) ;
				}
				if (string(argv[ii]) == "none")
				{
					debug = NO_DEBUG ;
					printInfo(MASTER, " Debug info level","NONE") ;
				}
				else if (string(argv[ii]) == "light")
				{
					debug = LIGHT_DEBUG ;
					printInfo(MASTER, " Debug info level","LIGHT") ;
				}
				else if (string(argv[ii]) == "mid")
				{
					debug = MID_DEBUG ;
					printInfo(MASTER, " Debug info level","MEDIUM") ;
				}
				else if (string(argv[ii]) == "full")
				{
					debug = FULL_DEBUG ;
					printInfo(MASTER, " Debug info level","FULL") ;
				}
				else
				{
					printError(" Invalid debug level", argv[ii]) ;
					return(RTN_CODE_KO) ;
				}
			}

			else if (string(argv[ii]) == "-dim")
			{
				ii++ ;
				if (ii >= argc)
				{
					printError(" parameter is needed after -dim") ;
					return(RTN_CODE_KO) ;
				}
				dim = (Dim_type) atoi(argv[ii]);
				printInfo(MASTER, " Space dimension", dim) ;

				if ((dim != DIM1)
						&& (dim != DIM2)
						&& (dim != DIM3))
				{
					printError(" Space dimension should be 1, 2 or 3") ;
					return(RTN_CODE_KO) ;
				}
			}

			else if (string(argv[ii]) == "-dt")
			{
				ii++ ;
				if (ii >= argc)
				{
					printError(" parameter is needed after -dt") ;
					return(RTN_CODE_KO) ;
				}
				dt = atof(argv[ii]);
				printInfo(MASTER, " Time step (s)\t", dt) ;
			}

			else if (string(argv[ii]) == "-fdOrder")
			{
				ii++ ;
				if (ii >= argc)
				{
					printError(" parameter is needed after -fdOrder") ;
					return(RTN_CODE_KO) ;
				}
				fdOrder = atoi(argv[ii]);
				printInfo(MASTER, " FD space Order\t", fdOrder) ;

				if ((fdOrder != 2)
						&& (fdOrder != 4)
						&& (fdOrder != 8)
						&& (fdOrder != 12)
						&& (fdOrder != 16))
				{
					printError(" FD space order should be 2, 4, 8, 12 or 16") ;
					return(RTN_CODE_KO) ;
				}
			}

			else if (string(argv[ii]) == "-gpuBlkSize")
			{
				ii++ ;
				if (ii >= argc)
				{
					printError(" parameter is needed after -gpuBlkSize") ;
					return(RTN_CODE_KO) ;
				}
				gpuBlkSize = atoi(argv[ii]);
				printInfo(MASTER, " gpuBlkSize\t", n1) ;
				if (gpuBlkSize <= 0)
				{
					printError(" gpuBlkSize should be > 0") ;
					return(RTN_CODE_KO) ;
				}
			}

			else if (string(argv[ii]) == "-gpuGridSize")
			{
				ii++ ;
				if (ii >= argc)
				{
					printError(" parameter is needed after -gpuGridSize") ;
					return(RTN_CODE_KO) ;
				}
				gpuGridSize = atoi(argv[ii]);
				printInfo(MASTER, " gpuGridSize\t", n1) ;
				if (gpuGridSize <= 0)
				{
					printError(" gpuGridSize should be > 0") ;
					return(RTN_CODE_KO) ;
				}
			}

			else if (string(argv[ii]) == "-n1")
			{
				ii++ ;
				if (ii >= argc)
				{
					printError(" parameter is needed after -n1") ;
					return(RTN_CODE_KO) ;
				}
				n1 = atoi(argv[ii]);
				printInfo(MASTER, " n1\t\t", n1) ;
				if (n1 <= 0)
				{
					printError(" n1 should be > 0") ;
					return(RTN_CODE_KO) ;
				}
			}

			else if (string(argv[ii]) == "-n2")
			{
				ii++ ;
				if (ii >= argc)
				{
					printError(" parameter is needed after -n2") ;
					return(RTN_CODE_KO) ;
				}
				n2 = atoi(argv[ii]);
				printInfo(MASTER, " n2\t\t", n2) ;
				if (n2 <= 0)
				{
					printError(" n2 should be > 0") ;
					return(RTN_CODE_KO) ;
				}
			}

			else if (string(argv[ii]) == "-n3")
			{
				ii++ ;
				if (ii >= argc)
				{
					printError(" parameter is needed after -n3") ;
					return(RTN_CODE_KO) ;
				}
				n3 = atoi(argv[ii]);
				printInfo(MASTER, " n3\t\t", n3) ;
				if (n3 <= 0)
				{
					printError(" n3 should be > 0") ;
					return(RTN_CODE_KO) ;
				}
			}

			else if (string(argv[ii]) == "-n1AddPad")
			{
				ii++ ;
				if (ii >= argc)
				{
					printError(" parameter is needed after -n1AddPad") ;
					return(RTN_CODE_KO) ;
				}
				n1AddPad = atoi(argv[ii]);
				printInfo(MASTER, " n1AddPad\t", n1AddPad) ;
				if (n1AddPad < 0)
				{
					printError(" n1AddPad should be >= 0") ;
					return(RTN_CODE_KO) ;
				}
			}

			else if (string(argv[ii]) == "-n2AddPad")
			{
				ii++ ;
				if (ii >= argc)
				{
					printError(" parameter is needed after -n2AddPad") ;
					return(RTN_CODE_KO) ;
				}
				n2AddPad = atoi(argv[ii]);
				printInfo(MASTER, " n2AddPad\t", n2AddPad) ;
				if (n2AddPad < 0)
				{
					printError(" n2AddPad should be >= 0") ;
					return(RTN_CODE_KO) ;
				}
			}

			else if (string(argv[ii]) == "-n3AddPad")
			{
				ii++ ;
				if (ii >= argc)
				{
					printError(" parameter is needed after -n3AddPad") ;
					return(RTN_CODE_KO) ;
				}
				n3AddPad = atoi(argv[ii]);
				printInfo(MASTER, " n3AddPad\t", n3AddPad) ;
				if (n3AddPad < 0)
				{
					printError(" n3AddPad should be >= 0") ;
					return(RTN_CODE_KO) ;
				}
			}

			else if (string(argv[ii]) == "-n1MulPad")
			{
				ii++ ;
				if (ii >= argc)
				{
					printError(" parameter is needed after -n1MulPad") ;
					return(RTN_CODE_KO) ;
				}
				n1MulPad = atoi(argv[ii]);
				printInfo(MASTER, " n1MulPad\t", n1MulPad) ;
				if (n1MulPad <= 0)
				{
					printError(" n1MulPad should be > 0") ;
					return(RTN_CODE_KO) ;
				}
			}

			else if (string(argv[ii]) == "-n2MulPad")
			{
				ii++ ;
				if (ii >= argc)
				{
					printError(" parameter is needed after -n2MulPad") ;
					return(RTN_CODE_KO) ;
				}
				n2MulPad = atoi(argv[ii]);
				printInfo(MASTER, " n2MulPad\t", n2MulPad) ;
				if (n2MulPad <= 0)
				{
					printError(" n2MulPad should be > 0") ;
					return(RTN_CODE_KO) ;
				}
			}

			else if (string(argv[ii]) == "-n3MulPad")
			{
				ii++ ;
				if (ii >= argc)
				{
					printError(" parameter is needed after -n3MulPad") ;
					return(RTN_CODE_KO) ;
				}
				n3MulPad = atoi(argv[ii]);
				printInfo(MASTER, " n3MulPad\t", n3MulPad) ;
				if (n3MulPad <= 0)
				{
					printError(" n3MulPad should be > 0") ;
					return(RTN_CODE_KO) ;
				}
			}

			else if (string(argv[ii]) == "-nt")
			{
				ii++ ;
				if (ii >= argc)
				{
					printError(" parameter is needed after -nt") ;
					return(RTN_CODE_KO) ;
				}
				nt = atoi(argv[ii]);
				if (nt <= 0)
				{
					printError(" nt should be > 0") ;
					return(RTN_CODE_KO) ;
				}
				printInfo(MASTER, " No. of time steps", nt) ;
			}

			else if (string(argv[ii]) == "-ntry")
			{
				ii++ ;
				if (ii >= argc)
				{
					printError(" parameter is needed after -ntry") ;
					return(RTN_CODE_KO) ;
				}
				ntry = atoi(argv[ii]);
				printInfo(MASTER, " No. of tries\t", ntry) ;
			}

			else if (string(argv[ii]) == "-nsub1")
			{
				ii++ ;
				if (ii >= argc)
				{
					printError(" parameter is needed after -nsub1") ;
					return(RTN_CODE_KO) ;
				}
				nsub1 = atoi(argv[ii]);
				printInfo(MASTER, " No. of nsub1\t", nsub1) ;
			}

			else if (string(argv[ii]) == "-nsub2")
			{
				ii++ ;
				if (ii >= argc)
				{
					printError(" parameter is needed after -nsub2") ;
					return(RTN_CODE_KO) ;
				}
				nsub2 = atoi(argv[ii]);
				printInfo(MASTER, " No. of nsub2\t", nsub2) ;
			}

			else if (string(argv[ii]) == "-nsub3")
			{
				ii++ ;
				if (ii >= argc)
				{
					printError(" parameter is needed after -nsub3") ;
					return(RTN_CODE_KO) ;
				}
				nsub3 = atoi(argv[ii]);
				printInfo(MASTER, " No. of nsub3\t", nsub3) ;
			}

			else if (string(argv[ii]) == "-param1")
			{
				ii++ ;
				if (ii >= argc)
				{
					printError(" parameter is needed after -param1") ;
					return(RTN_CODE_KO) ;
				}
				param1 = atof(argv[ii]);
				printInfo(MASTER, " Parameter 1\t", param1) ;
			}

			else if (string(argv[ii]) == "-param2")
			{
				ii++ ;
				if (ii >= argc)
				{
					printError(" parameter is needed after -param2") ;
					return(RTN_CODE_KO) ;
				}
				param2 = atof(argv[ii]);
				printInfo(MASTER, " Parameter 2\t", param2) ;
			}

			else if (string(argv[ii]) == "-param3")
			{
				ii++ ;
				if (ii >= argc)
				{
					printError(" parameter is needed after -param3") ;
					return(RTN_CODE_KO) ;
				}
				param3 = atof(argv[ii]);
				printInfo(MASTER, " Parameter 3\t", param3) ;
			}

			else if (string(argv[ii]) == "-param4")
			{
				ii++ ;
				if (ii >= argc)
				{
					printError(" parameter is needed after -param4") ;
					return(RTN_CODE_KO) ;
				}
				param4 = atof(argv[ii]);
				printInfo(MASTER, " Parameter 4\t", param4) ;
			}

			else if (string(argv[ii]) == "-propagator")
			{
				ii++ ;
				if (ii >= argc)
				{
					printError(" parameter is needed after -propagator") ;
					return(RTN_CODE_KO) ;
				}
				propagator = argv[ii];
				printInfo(MASTER, " Propagator\t", propagator) ;
			}

			else if (string(argv[ii]) == "-ratioCFL")
			{
				ii++ ;
				if (ii >= argc)
				{
					printError(" parameter is needed after -ratioCFL") ;
					return(RTN_CODE_KO) ;
				}
				ratioCFL = atof(argv[ii]);
				printInfo(MASTER, " Ratio CFL\t", ratioCFL) ;
			}

			else if (string(argv[ii]) == "-snapInc")
			{
				ii++ ;
				if (ii >= argc)
				{
					printError(" parameter is needed after -snapInc") ;
					return(RTN_CODE_KO) ;
				}
				snapInc = atoi(argv[ii]);
				printInfo(MASTER, " Snapshots inc. (steps)", snapInc) ;
				if (snapInc <= 0)
				{
					printError(" snapInc should be > 0") ;
					return(RTN_CODE_KO) ;
				}
			}

			else if (string(argv[ii]) == "-snapDt")
			{
				ii++ ;
				if (ii >= argc)
				{
					printError(" parameter is needed after -snapDt") ;
					return(RTN_CODE_KO) ;
				}
				snapDt = atof(argv[ii]);
				printInfo(MASTER, " Snapshots inc. (sec.)", snapDt) ;
				if (snapDt <= 0.0)
				{
					printError(" snapDt should be > 0") ;
					return(RTN_CODE_KO) ;
				}
			}

			else if (string(argv[ii]) == "-testCase")
			{
				ii++ ;
				if (ii >= argc)
				{
					printError(" parameter is needed after -testCase") ;
					return(RTN_CODE_KO) ;
				}
				testCaseName = argv[ii];
				printInfo(MASTER, " TestCase name\t", testCaseName) ;
			}

			else if (string(argv[ii]) == "-testMode")
			{
				ii++ ;
				if (ii >= argc)
				{
					printError(" parameter is needed after -testMode") ;
					return(RTN_CODE_KO) ;
				}
				testMode = argv[ii];
				printInfo(MASTER, " Test mode\t", testMode) ;
			}

			else if (string(argv[ii]) == "-tmax")
			{
				ii++ ;
				if (ii >= argc)
				{
					printError(" parameter is needed after -tmax") ;
					return(RTN_CODE_KO) ;
				}
				tmax = atof(argv[ii]);
				printInfo(MASTER, " Tmax\t\t", tmax) ;
				if (tmax <= 0.0)
				{
					printError("tmax should be > 0") ;
					return(RTN_CODE_KO);
				}
			}


			else if (string(argv[ii]) == "-writeGrid")
			{
				printInfo(MASTER, " Write grids on disk","YES") ;
				writeGrid = true ;
			}

			else
			{
				printError(" Unknown argument ", argv[ii]) ;
				return(RTN_CODE_KO);
			}
			ii++ ;
		}
	}
	printDebug(MID_DEBUG, "OUT Config::parse_argument");
	return(RTN_CODE_OK);
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Config::info(void)
{
	printDebug(MID_DEBUG, "IN Config::info");

	// print main config. param.
	print_blank() ;
	print_line5() ;
	Myint nproc = hpcscan::nproc_world ;
	printInfo(MASTER, " Configuration parameters") ;
	printInfo(MASTER, " TestCase name\t", testCaseName) ;

	if (testCaseName.compare(DEFAULT_TEST_CASE) == 0)
	{
		printWarning(" Please specify test case name with -testCase") ;
		return(RTN_CODE_KO);
	}

	printInfo(MASTER, " No. of MPI process", nproc) ;
	printInfo(MASTER, " No. subdomains axis1", nsub1) ;
	printInfo(MASTER, " No. subdomains axis2", nsub2) ;
	printInfo(MASTER, " No. subdomains axis3", nsub3) ;
	printInfo(MASTER, " OpenMP threads / MPI", omp_get_max_threads()) ;
	printInfo(MASTER, " Number of tries", ntry) ;
	printInfo(MASTER, " FD space order\t", fdOrder) ;
	print_line5() ;

	// check sub-domain decomposition
	if (nproc != (nsub1*nsub2*nsub3))
	{
		printError("Number subdomains inconsistent with nproc") ;
		return(RTN_CODE_KO);
	}
	if (dim == DIM1)
	{
		if (nsub2 > 1)
		{
			printError("nsub2 > 1 not allowed for 1D") ;
			return(RTN_CODE_KO);
		}
		if (nsub3 > 1)
		{
			printError("nsub3 > 1 not allowed for 1D") ;
			return(RTN_CODE_KO);
		}
		if ((n1 / nsub1) < MIN_NPOINT_NSUB)
		{
			printError("(n1 / nsub1) < MIN_NPOINT_NSUB, enlarge n1") ;
			printInfo(MASTER, "MIN_NPOINT_NSUB", MIN_NPOINT_NSUB) ;
			return(RTN_CODE_KO);
		}
	}
	else if (dim == DIM2)
	{
		if (nsub3 > 1)
		{
			printError("nsub3 > 1 not allowed for 2D") ;
			return(RTN_CODE_KO);
		}
		if ((n1 / nsub1) < MIN_NPOINT_NSUB)
		{
			printError("(n1 / nsub1) < MIN_NPOINT_NSUB, enlarge n1") ;
			printInfo(MASTER, "MIN_NPOINT_NSUB", MIN_NPOINT_NSUB) ;
			return(RTN_CODE_KO);
		}
		if ((n2 / nsub2) < MIN_NPOINT_NSUB)
		{
			printError("(n2 / nsub2) < MIN_NPOINT_NSUB, enlarge n2") ;
			printInfo(MASTER, "MIN_NPOINT_NSUB", MIN_NPOINT_NSUB) ;
			return(RTN_CODE_KO);
		}
	}
	else if (dim == DIM3)
	{
		if ((n1 / nsub1) < MIN_NPOINT_NSUB)
		{
			printError("(n1 / nsub1) < MIN_NPOINT_NSUB, enlarge n1") ;
			printInfo(MASTER, "MIN_NPOINT_NSUB", MIN_NPOINT_NSUB) ;
			return(RTN_CODE_KO);
		}
		if ((n2 / nsub2) < MIN_NPOINT_NSUB)
		{
			printError("(n2 / nsub2) < MIN_NPOINT_NSUB, enlarge n2") ;
			printInfo(MASTER, "MIN_NPOINT_NSUB", MIN_NPOINT_NSUB) ;
			return(RTN_CODE_KO);
		}
		if ((n3 / nsub3) < MIN_NPOINT_NSUB)
		{
			printError("(n3 / nsub3) < MIN_NPOINT_NSUB, enlarge n3") ;
			printInfo(MASTER, "MIN_NPOINT_NSUB", MIN_NPOINT_NSUB) ;
			return(RTN_CODE_KO);
		}
	}

	printDebug(MID_DEBUG, "OUT Config::info");
	return(RTN_CODE_OK);
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Config::initialize(void)
{
	printDebug(MID_DEBUG, "IN Config::initialize");

	// open log files
	if (debug > NO_DEBUG)
	{
		string file_name = "hpcscan.debug.proc" + to_string(myid_world) + ".log";
		debugLogFile.open(file_name);
	}

	printDebug(MID_DEBUG, "OUT Config::initialize");
	return(RTN_CODE_OK);
}

//-------------------------------------------------------------------------------------------------------

Rtn_code Config::finalize(Rtn_code rtnCode)
{
	printDebug(MID_DEBUG, "IN Config::finalize");

	print_end_of_output_report(rtnCode) ;

	// close log files
	if (debug > NO_DEBUG)
	{
		debugLogFile.close() ;
	}

	//MPI_Barrier(MPI_COMM_WORLD) ;

	// finalize MPI environment
	MPI_Finalize() ;

	printDebug(MID_DEBUG, "OUT Config::finalize");
	return(RTN_CODE_OK);
}


} // namespace hpcscan
