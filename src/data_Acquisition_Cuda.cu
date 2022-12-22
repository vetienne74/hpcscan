
#include "data_Acquisition_Cuda.h"

#include <cassert>
#include <cmath>

#include "config.h"
#include "output_report.h"
#include "global.h" 

namespace hpcscan {

//-------------------------------------------------------------------------------------------------------

DataAcquisition_Cuda::DataAcquisition_Cuda(void)
{
    printDebug(MID_DEBUG, "IN DataAcquisition_Cuda::DataAcquisition_Cuda");

    printDebug(MID_DEBUG, "OUT DataAcquisition_Cuda::DataAcquisition_Cuda");
}

//-------------------------------------------------------------------------------------------------------

DataAcquisition_Cuda::~DataAcquisition_Cuda(void)
{
	printDebug(MID_DEBUG, "IN DataAcquisition_Cuda::~DataAcquisition_Cuda");	

	printDebug(MID_DEBUG, "OUT DataAcquisition_Cuda::~DataAcquisition_Cuda");
}

//-------------------------------------------------------------------------------------------------------

Rtn_code DataAcquisition_Cuda::initialize(DataAcqui_type acquiType, Grid& coefGrid, Myint nt)
{
	printDebug(MID_DEBUG, "In DataAcquisition_Cuda::initialize") ;     

    // allocate grid to store collected traces

    printDebug(MID_DEBUG, "Out DataAcquisition_Cuda::initialize") ;
    return(RTN_CODE_OK) ;
}

//-------------------------------------------------------------------------------------------------------

Rtn_code DataAcquisition_Cuda::appliSourceTerm(Grid &grid, Modeling_type modType, Myint it, Myfloat dt)
{
    printDebug(MID_DEBUG, "In DataAcquisition_Cuda::appliSourceTerm");    

    printDebug(MID_DEBUG, "Out DataAcquisition_Cuda::appliSourceTerm");
    return (RTN_CODE_OK);
}

//-------------------------------------------------------------------------------------------------------

Rtn_code DataAcquisition_Cuda::recordTrace(Grid &grid, Modeling_type modType, Myint it, Myfloat dt)
{
    printDebug(MID_DEBUG, "In DataAcquisition_Cuda::recordTrace");  

    printDebug(MID_DEBUG, "Out DataAcquisition_Cuda::recordTrace");
    return (RTN_CODE_OK);
}

} // namespace hpcscan
