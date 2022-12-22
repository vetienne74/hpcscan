
//-------------------------------------------------------------------------------------------------------
// Data acquisition factory used to create specialized Data_acquisition objects
//-------------------------------------------------------------------------------------------------------

#include "data_Acquisition_Factory.h"

#include "constant.h"
#include "data_Acquisition.h"
#include "data_Acquisition_Cuda.h"
#include "data_Acquisition_NEC.h"
#include "data_Acquisition_NEC_SCA.h"
#include "output_report.h"

using namespace std;

namespace hpcscan {

//-------------------------------------------------------------------------------------------------------

shared_ptr<DataAcquisition> DataAcquisition_Factory::create(string testMode)
{
	printDebug(MID_DEBUG, "IN DataAcquisition_Factory::create");

	DataAcquisition *acqui = nullptr ;

	if ((testMode.compare(GRID_MODE_BASELINE) == 0) || (testMode.compare(GRID_MODE_CACHEBLK) == 0))
	{
		acqui = new DataAcquisition();
	}

#ifdef __NEC__
	else if ((testMode.compare(GRID_MODE_NEC) == 0))
	{
		acqui = new DataAcquisition_NEC(); 
	}
	else if ((testMode.compare(GRID_MODE_NEC_SCA) == 0))
	{
		acqui = new DataAcquisition_NEC_SCA(); 
	}
#endif

#ifdef __CUDA__
	else if ((testMode.compare(GRID_MODE_CUDA) == 0) || (testMode.compare(GRID_MODE_CUDA_OPTIM)) || (testMode.compare(GRID_MODE_CUDA_REF)))
	{
		acqui = new DataAcquisition_Cuda();
	}
#endif

// #ifdef __NEC_SCA__
// 	else if ((testMode.compare(GRID_MODE_NEC_SCA) == 0))
// 	{
// 		acqui = new DataAcquisition_NEC_SCA(); 
// 	}
// #endif

	else
	{
		printError("IN DataAcquisition_Factory::create, invalid testMode", testMode) ;
	}

	printDebug(MID_DEBUG, "OUT DataAcquisition_Factory::create");
	if (acqui != nullptr)
	{
		return shared_ptr<DataAcquisition>(acqui) ;
	}
	else
	{
		return nullptr ;
	}
}

} // namespace hpcscan
