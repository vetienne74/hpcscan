
//-------------------------------------------------------------------------------------------------------
// This data acquisition is activated with command line option -testMode NEC_SCA
// Derived class from DataAcquisition
// NEC implementation (target SX-Aurora)
//-------------------------------------------------------------------------------------------------------


#ifndef HPCSCAN_DATA_ACQUISITION_NEC_SCA_H_
#define HPCSCAN_DATA_ACQUISITION_NEC_SCA_H_

#include "data_Acquisition_NEC.h"
#include "type_def.h"

using namespace std;

namespace hpcscan
{

    class DataAcquisition_NEC_SCA : public DataAcquisition_NEC
    {
    public:

        // write traces on disk
        //virtual Rtn_code writeTrace(string fileName) ;        

    private:
       
    };

} // namespace hpcscan
#endif