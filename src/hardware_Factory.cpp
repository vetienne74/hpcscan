
//-------------------------------------------------------------------------------------------------------
// Hardware factory used to create specialized Hardware objects
// All hardware implementations are referenced
//-------------------------------------------------------------------------------------------------------

#include "hardware_Factory.h"

#include "constant.h"
#ifdef __CUDA__
#include "hardware_Cuda.h"
#endif
#ifdef __DPCPP__
#include "hardware_DPCPP.h"
#endif
#ifdef __HIP__
#include "hardware_Hip.h"
#endif
#ifdef __NEC__
#include "hardware_NEC.h"
#endif
//#ifdef __OPENACC__
//#include "hardware_OpenAcc.h"
//#endif
#include "output_report.h"

using namespace std;

namespace hpcscan {

//-------------------------------------------------------------------------------------------------------

shared_ptr<Hardware> Hardware_Factory::create(string gridMode)
{
	printDebug(MID_DEBUG, "IN Hardware_Factory::create");

	Hardware *Rhardware = nullptr ;

	if (gridMode.compare(GRID_MODE_BASELINE) == 0)
	{
		Rhardware = new Hardware(gridMode) ;
	}
	else if (gridMode.compare(GRID_MODE_CACHEBLK) == 0)
	{
		Rhardware = new Hardware(gridMode) ;
	}

#ifdef __CUDA__
	else if (gridMode.compare(GRID_MODE_CUDA) == 0)
	{
		Rhardware = new Hardware_Cuda(gridMode) ;
	}
	else if (gridMode.compare(GRID_MODE_CUDA_OPTIM) == 0)
	{
		Rhardware = new Hardware_Cuda(gridMode) ;
	}
	else if (gridMode.compare(GRID_MODE_CUDA_REF) == 0)
	{
		Rhardware = new Hardware_Cuda(gridMode) ;
	}
#endif

#ifdef __DPCPP__
	else if (gridMode.compare(GRID_MODE_DPCPP) == 0)
	{
		Rhardware = new Hardware_DPCPP(gridMode) ;
	}
#endif

#ifdef __HIP__
	else if (gridMode.compare(GRID_MODE_HIP) == 0)
	{
		Rhardware = new Hardware_Hip(gridMode) ;
	}
	else if (gridMode.compare(GRID_MODE_HIP_OPTIM) == 0)
	{
		Rhardware = new Hardware_Hip(gridMode) ;
	}
#endif

#ifdef __NEC__
	else if (gridMode.compare(GRID_MODE_NEC) == 0)
	{
		Rhardware = new Hardware_NEC(gridMode) ;
	}
	else if (gridMode.compare(GRID_MODE_NEC_SCA) == 0)
	{
		Rhardware = new Hardware_NEC(gridMode) ;
	}
#endif

//#ifdef __OPENACC__
//	else if (gridMode.compare(GRID_MODE_OPENACC) == 0)
//	{
//		Rhardware = new Hardware_OpenAcc(gridMode) ;
//	}
//#endif

	else
	{
		printError("IN Hardware_Factory::create, not supported or invalid gridMode", gridMode) ;
	}

	printDebug(MID_DEBUG, "OUT Hardware_Factory::create");
	if (Rhardware != nullptr)
	{
		return shared_ptr<Hardware>(Rhardware) ;
	}
	else
	{
		printError("IN Hardware_Factory::create, return nullptr", gridMode) ;
		return nullptr ;
	}
}

} // namespace hpcscan
