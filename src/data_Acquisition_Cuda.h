
//-------------------------------------------------------------------------------------------------------
// This data acquisition is activated with command line option -testMode CUDA
// Derived class from DataAcquisition
// CUDA implementation (target GPU)
//-------------------------------------------------------------------------------------------------------


#ifndef HPCSCAN_DATA_ACQUISITION_CUDA_H_
#define HPCSCAN_DATA_ACQUISITION_CUDA_H_

#include "data_Acquisition.h"
#include "type_def.h"

using namespace std;

namespace hpcscan
{

    class DataAcquisition_Cuda : public DataAcquisition
    {
    public:
        // constructor
        DataAcquisition_Cuda(void);

        // destructor
        ~DataAcquisition_Cuda(void);

        // initialize
        Rtn_code initialize(DataAcqui_type acquiType, Grid &coefGrid, Myint nt);

        // apply source terms
        Rtn_code appliSourceTerm(Grid &grid, Modeling_type modType, Myint it, Myfloat dt) ;

        // record traces 
        Rtn_code recordTrace(Grid &grid, Modeling_type modType, Myint it, Myfloat dt) ;        

    private:
       
    };

} // namespace hpcscan
#endif