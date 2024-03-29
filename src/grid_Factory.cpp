
//-------------------------------------------------------------------------------------------------------
// Grid factory used to create specialized Grid objects
// Unique location where all grid implementations are referenced
//-------------------------------------------------------------------------------------------------------

#include "grid_Factory.h"

#include "constant.h"
#include "grid_CacheBlk.h"
#ifdef __CUDA__
#include "grid_Cuda.h"
#include "grid_Cuda_Optim.h"
#include "grid_Cuda_Ref.h"
#endif
#ifdef __DPCPP__
#include "grid_DPCPP.h"
#endif
#ifdef __HIP__
#include "grid_Hip.h"
#include "grid_Hip_Optim.h"
#endif
#ifdef __NEC__
#include "grid_NEC.h"
#include "grid_NEC_SCA.h"
#endif
#ifdef __OPENACC__
#include "grid_OpenAcc.h"
#endif
#include "output_report.h"

using namespace std;

namespace hpcscan {

//-------------------------------------------------------------------------------------------------------

shared_ptr<Grid> Grid_Factory::create(string gridMode, Grid_type gridType)
{
	printDebug(MID_DEBUG, "IN Grid_Factory::create");

	Grid *retGrid = nullptr ;

	if (gridMode.compare(GRID_MODE_BASELINE) == 0)
	{
		retGrid = new Grid(gridType) ;
	}
	else if (gridMode.compare(GRID_MODE_CACHEBLK) == 0)
	{
		retGrid = new Grid_CacheBlk(gridType) ;
	}

#ifdef __CUDA__
	else if (gridMode.compare(GRID_MODE_CUDA) == 0)
	{
		retGrid = new Grid_Cuda(gridType) ;
	}
	else if (gridMode.compare(GRID_MODE_CUDA_OPTIM) == 0)
	{
		retGrid = new Grid_Cuda_Optim(gridType) ;
	}
	else if (gridMode.compare(GRID_MODE_CUDA_REF) == 0)
	{
		retGrid = new Grid_Cuda_Ref(gridType) ;
	}
#endif

#ifdef __DPCPP__
	else if (gridMode.compare(GRID_MODE_DPCPP) == 0)
	{
		retGrid = new Grid_DPCPP(gridType) ;
	}
#endif

#ifdef __HIP__
	else if (gridMode.compare(GRID_MODE_HIP) == 0)
	{
		retGrid = new Grid_Hip(gridType) ;
	}
	else if (gridMode.compare(GRID_MODE_HIP_OPTIM) == 0)
	{
		retGrid = new Grid_Hip_Optim(gridType) ;
	}
#endif

#ifdef __NEC__
	else if (gridMode.compare(GRID_MODE_NEC) == 0)
	{
		retGrid = new Grid_NEC(gridType) ;
	}
	else if (gridMode.compare(GRID_MODE_NEC_SCA) == 0)
	{
		retGrid = new Grid_NEC_SCA(gridType) ;
	}
#endif

#ifdef __OPENACC__
	else if (gridMode.compare(GRID_MODE_OPENACC) == 0)
	{
		retGrid = new Grid_OpenAcc(gridType) ;
	}
#endif

	else
	{
		printError("IN Grid_Factory::create, not supported or invalid gridMode", gridMode) ;
	}

	printDebug(MID_DEBUG, "OUT Grid_Factory::create");
	if (retGrid != nullptr)
	{
		return shared_ptr<Grid>(retGrid) ;
	}
	else
	{
		printError("IN Grid_Factory::create, return nullptr", gridMode) ;
		return nullptr ;
	}
}

//-------------------------------------------------------------------------------------------------------

shared_ptr<Grid> Grid_Factory::create(string gridMode, Grid_type gridType, Dim_type dim,
		Myint64 n1, Myint64 n2, Myint64 n3)
{
	printDebug(MID_DEBUG, "IN Grid_Factory::create");

	Grid *retGrid = nullptr ;

	if (gridMode.compare(GRID_MODE_BASELINE) == 0)
	{
		retGrid = new Grid(gridType, dim, n1, n2, n3) ;
	}
	else if (gridMode.compare(GRID_MODE_CACHEBLK) == 0)
	{
		retGrid = new Grid_CacheBlk(gridType, dim, n1, n2, n3) ;
	}

#ifdef __CUDA__
	else if (gridMode.compare(GRID_MODE_CUDA) == 0)
	{
		retGrid = new Grid_Cuda(gridType, dim, n1, n2, n3) ;
	}
	else if (gridMode.compare(GRID_MODE_CUDA_OPTIM) == 0)
	{
		retGrid = new Grid_Cuda_Optim(gridType, dim, n1, n2, n3) ;
	}
	else if (gridMode.compare(GRID_MODE_CUDA_REF) == 0)
	{
		retGrid = new Grid_Cuda_Ref(gridType, dim, n1, n2, n3) ;
	}
#endif

#ifdef __DPCPP__
	else if (gridMode.compare(GRID_MODE_DPCPP) == 0)
	{
		retGrid = new Grid_DPCPP(gridType, dim, n1, n2, n3) ;
	}
#endif

#ifdef __HIP__
	else if (gridMode.compare(GRID_MODE_HIP) == 0)
	{
		retGrid = new Grid_Hip(gridType, dim, n1, n2, n3) ;
	}
	else if (gridMode.compare(GRID_MODE_HIP_OPTIM) == 0)
	{
		retGrid = new Grid_Hip_Optim(gridType, dim, n1, n2, n3) ;
	}
#endif

#ifdef __NEC__
	else if (gridMode.compare(GRID_MODE_NEC) == 0)
	{
		retGrid = new Grid_NEC(gridType, dim, n1, n2, n3) ;
	}
	else if (gridMode.compare(GRID_MODE_NEC_SCA) == 0)
	{
		retGrid = new Grid_NEC_SCA(gridType, dim, n1, n2, n3) ;
	}
#endif

#ifdef __OPENACC__
	else if (gridMode.compare(GRID_MODE_OPENACC) == 0)
	{
		retGrid = new Grid_OpenAcc(gridType, dim, n1, n2, n3) ;
	}
#endif

	else
	{
		printError("IN Grid_Factory::create, not supported or invalid gridMode") ;
	}

	printDebug(MID_DEBUG, "OUT Grid_Factory::create");
	if (retGrid != nullptr)
	{
		return shared_ptr<Grid>(retGrid) ;
	}
	else
	{
		return nullptr ;
	}
}


} // namespace hpcscan
