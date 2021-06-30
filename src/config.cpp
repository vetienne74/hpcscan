
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
static const string    DEFAULT_DPCPP_SELECTOR = "Host" ;
static const Myfloat64 DEFAULT_DT             = 0.0 ; // stable dt will be used
static const Myint     DEFAULT_FD_ORDER       = 8 ;
static const Myint     DEFAULT_GPU_BLKSIZE    = 256 ;
static const Myint     DEFAULT_GPU_BLKSIZE1   = 32 ;
static const Myint     DEFAULT_GPU_BLKSIZE2   = 2 ;
static const Myint     DEFAULT_GPU_BLKSIZE3   = 16 ;
static const Myint     DEFAULT_GPU_GRIDSIZE   = 512 ;
static const bool      DEFAULT_GPU_MPI_AWARE  = false ;
static const Myfloat64 DEFAULT_H              = PI / 30 ;
static const Myfloat64 DEFAULT_HW_COUNTER_DT  = 0 ; // no update hardware counters
static const Myint     DEFAULT_INNER_N1       = 61 ;
static const Myint     DEFAULT_INNER_N2       = 61 ;
static const Myint     DEFAULT_INNER_N3       = 61 ;
static const Myint     DEFAULT_N1_ADD_PAD     = UNSPECIFIED ;
static const Myint     DEFAULT_N2_ADD_PAD     = UNSPECIFIED ;
static const Myint     DEFAULT_N3_ADD_PAD     = UNSPECIFIED ;
static const Myint     DEFAULT_N1_MUL_PAD     = UNSPECIFIED ;
static const Myint     DEFAULT_N2_MUL_PAD     = UNSPECIFIED ;
static const Myint     DEFAULT_N3_MUL_PAD     = UNSPECIFIED ;
static const Myint     DEFAULT_N1_OFFSET      = 0 ;
static const Myint     DEFAULT_N2_OFFSET      = 0 ;
static const Myint     DEFAULT_N3_OFFSET      = 0 ;
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
	dpcppSelect  = DEFAULT_DPCPP_SELECTOR ;
	fdOrder      = DEFAULT_FD_ORDER ;
	gpuBlkSize   = DEFAULT_GPU_BLKSIZE ;
	gpuBlkSize1  = DEFAULT_GPU_BLKSIZE1 ;
	gpuBlkSize2  = DEFAULT_GPU_BLKSIZE2 ;
	gpuBlkSize3  = DEFAULT_GPU_BLKSIZE3 ;
	gpuGridSize  = DEFAULT_GPU_GRIDSIZE ;
	gpuMpiAware  = DEFAULT_GPU_MPI_AWARE ;
	h            = DEFAULT_H ;
	hw           = nullptr ; // hardware is set in initialize()
	hwCounterDt  = DEFAULT_HW_COUNTER_DT ;
	n1           = DEFAULT_INNER_N1 ;
	n2           = DEFAULT_INNER_N2 ;
	n3           = DEFAULT_INNER_N3 ;
	n1AddPad     = DEFAULT_N1_ADD_PAD ;
	n2AddPad     = DEFAULT_N2_ADD_PAD ;
	n3AddPad     = DEFAULT_N3_ADD_PAD ;
	n1MulPad     = DEFAULT_N1_MUL_PAD ;
	n2MulPad     = DEFAULT_N2_MUL_PAD ;
	n3MulPad     = DEFAULT_N3_MUL_PAD ;
	n1Offset     = DEFAULT_N1_OFFSET ;
	n2Offset     = DEFAULT_N2_OFFSET ;
	n3Offset     = DEFAULT_N3_OFFSET ;
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
			printInfo(MASTER, " -dpcppSelect <char>  = device selector for DPC++") ;
			printInfo(MASTER, "     Host             * Launch kernels on Host (DEFAULT)") ;
			printInfo(MASTER, "     CPU              * Launch kernels on CPU") ;
			printInfo(MASTER, "     GPU              * Launch kernels on GPU") ;
			printInfo(MASTER, "     FPGA             * Launch kernels on FPGA") ;
			printInfo(MASTER, " -dt <float>          = time step (s) for propagator") ;
			printInfo(MASTER, " -fdOrder <int>       = spatial FD order [2,4,6,8,10,12,14 or 16]") ;
			printInfo(MASTER, " -gpuBlkSize <int>    = GPU, no. of threads / 1D block ") ;
			printInfo(MASTER, " -gpuBlkSize1 <int>   = GPU, no. of threads / 3D block axis 1") ;
			printInfo(MASTER, " -gpuBlkSize2 <int>   = GPU, no. of threads / 3D block axis 2") ;
			printInfo(MASTER, " -gpuBlkSize3 <int>   = GPU, no. of threads / 3D block axis 3") ;
			printInfo(MASTER, " -gpuGridSize <int>   = GPU, no. of 1D blocks per grid") ;
			printInfo(MASTER, " -gpuMpiAware         = use MPI GPU-aware library") ;
			printInfo(MASTER, " -help or -h          = list of command line parameters") ;
			printInfo(MASTER, " -hwCounterDt <float> = interval between hardware counter (DEFAULT 0=no update)") ;
			printInfo(MASTER, " -n1 <int>            = inner domain size axis 1 [grid pts]") ;
			printInfo(MASTER, " -n2 <int>            = inner domain size axis 2 [grid pts]") ;
			printInfo(MASTER, " -n3 <int>            = inner domain size axis 3 [grid pts]") ;
			printInfo(MASTER, " -n1AddPad <int>      = add N points to grid at the end of axis 1") ;
			printInfo(MASTER, " -n2AddPad <int>      = add N points to grid at the end of axis 2") ;
			printInfo(MASTER, " -n3AddPad <int>      = add N points to grid at the end of axis 3") ;
			printInfo(MASTER, " -n1MulPad <int>      = add some points to grid at end of axis 1 to be multiple of N") ;
			printInfo(MASTER, " -n2MulPad <int>      = add some points to grid at end of axis 1 to be multiple of N") ;
			printInfo(MASTER, " -n3MulPad <int>      = add some points to grid at end of axis 1 to be multiple of N") ;
			printInfo(MASTER, " -n1Offset <int>      = add N points to grid at the beginning of axis 1") ;
			printInfo(MASTER, " -n2Offset <int>      = add N points to grid at the beginning of axis 2") ;
			printInfo(MASTER, " -n3Offset <int>      = add N points to grid at the beginning of axis 3") ;
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
			printInfo(MASTER, "     Baseline         * Generic CPU reference implementation (DEFAULT test mode)") ;
			printInfo(MASTER, "     CacheBlk         * Generic CPU optimized with cache blocking techniques") ;
			printInfo(MASTER, "     CUDA             * NVIDIA GPU CUDA regular implementation") ;
			printInfo(MASTER, "     CUDA_Opt         * NVIDIA GPU CUDA optimized implementation") ;
			printInfo(MASTER, "     CUDA_ref         * NVIDIA GPU CUDA reference implementation (for developers)") ;
			printInfo(MASTER, "     DPC++            * Intel CPU/GPU/FPGA regular DPC++ implementation") ;
			printInfo(MASTER, "     HIP              * AMD GPU HIP regular implementation") ;
			printInfo(MASTER, "     HIP_Opt          * AMD GPU HIP optimized implementation") ;
			printInfo(MASTER, "     NEC              * NEC SX-Aurora with NEC compiler directives") ;
			printInfo(MASTER, "     NEC_SCA          * NEC SX-Aurora with NEC Stencil Code Accelerator Library") ;
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

			else if (string(argv[ii]) == "-dpcppSelect")
			{
				ii++ ;
				if (ii >= argc)
				{
					printError(" parameter is needed after -dpcppSelect") ;
					return(RTN_CODE_KO) ;
				}
				dpcppSelect = argv[ii];
				printInfo(MASTER, " DPC++ selector\t", dpcppSelect) ;
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
						&& (fdOrder != 6)
						&& (fdOrder != 8)
						&& (fdOrder != 10)
						&& (fdOrder != 12)
						&& (fdOrder != 14)
						&& (fdOrder != 16))
				{
					printError(" FD space order should be 2,4,8,10,12,14 or 16") ;
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

			else if (string(argv[ii]) == "-gpuBlkSize1")
			{
				ii++ ;
				if (ii >= argc)
				{
					printError(" parameter is needed after -gpuBlkSize1") ;
					return(RTN_CODE_KO) ;
				}
				gpuBlkSize1 = atoi(argv[ii]);
				printInfo(MASTER, " gpuBlkSize1\t", n1) ;
				if (gpuBlkSize1 <= 0)
				{
					printError(" gpuBlkSize1 should be > 0") ;
					return(RTN_CODE_KO) ;
				}
			}

			else if (string(argv[ii]) == "-gpuBlkSize2")
			{
				ii++ ;
				if (ii >= argc)
				{
					printError(" parameter is needed after -gpuBlkSize2") ;
					return(RTN_CODE_KO) ;
				}
				gpuBlkSize2 = atoi(argv[ii]);
				printInfo(MASTER, " gpuBlkSize2\t", n1) ;
				if (gpuBlkSize2 <= 0)
				{
					printError(" gpuBlkSize2 should be > 0") ;
					return(RTN_CODE_KO) ;
				}
			}

			else if (string(argv[ii]) == "-gpuBlkSize3")
			{
				ii++ ;
				if (ii >= argc)
				{
					printError(" parameter is needed after -gpuBlkSize3") ;
					return(RTN_CODE_KO) ;
				}
				gpuBlkSize3 = atoi(argv[ii]);
				printInfo(MASTER, " gpuBlkSize3\t", n1) ;
				if (gpuBlkSize3 <= 0)
				{
					printError(" gpuBlkSize3 should be > 0") ;
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

			else if (string(argv[ii]) == "-gpuMpiAware")
			{
				gpuMpiAware = true ;
				printInfo(MASTER, " gpuMpiAware\t", "ON") ;
			}

			else if (string(argv[ii]) == "-hwCounterDt")
			{
				ii++ ;
				if (ii >= argc)
				{
					printError(" parameter is needed after -hwCounterDt") ;
					return(RTN_CODE_KO) ;
				}
				hwCounterDt = atof(argv[ii]);
				printInfo(MASTER, " Dt hw counter (s)", hwCounterDt) ;
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

			else if (string(argv[ii]) == "-n1Offset")
			{
				ii++ ;
				if (ii >= argc)
				{
					printError(" parameter is needed after -n1Offset") ;
					return(RTN_CODE_KO) ;
				}
				n1Offset = atoi(argv[ii]);
				printInfo(MASTER, " n1Offset\t", n1Offset) ;
				if (n1Offset < 0)
				{
					printError(" n1Offset should be >= 0") ;
					return(RTN_CODE_KO) ;
				}
			}

			else if (string(argv[ii]) == "-n2Offset")
			{
				ii++ ;
				if (ii >= argc)
				{
					printError(" parameter is needed after -n2Offset") ;
					return(RTN_CODE_KO) ;
				}
				n2Offset = atoi(argv[ii]);
				printInfo(MASTER, " n2Offset\t", n2Offset) ;
				if (n2Offset < 0)
				{
					printError(" n2Offset should be >= 0") ;
					return(RTN_CODE_KO) ;
				}
			}

			else if (string(argv[ii]) == "-n3Offset")
			{
				ii++ ;
				if (ii >= argc)
				{
					printError(" parameter is needed after -n3Offset") ;
					return(RTN_CODE_KO) ;
				}
				n3Offset = atoi(argv[ii]);
				printInfo(MASTER, " n3Offset\t", n3Offset) ;
				if (n3Offset < 0)
				{
					printError(" n3Offset should be >= 0") ;
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
	Myint nproc = hpcscan::nMpiProc ;
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
		string file_name = "hpcscan.debug.proc" + to_string(myMpiRank) + ".log";
		debugLogFile.open(file_name);
	}

	// initialize Hardware object
	{
		hw = Hardware_Factory::create(Config::Instance()->testMode) ;
		if (hw == nullptr)
		{
			printError("IN Config::initialize, can not initialize hardware") ;
			return(RTN_CODE_KO);
		}
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
