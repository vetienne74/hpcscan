
//-------------------------------------------------------------------------------------------------------
// This data acquisition is activated with command line option -testMode NEC
// Derived class from DataAcquisition
// NEC implementation (target SX-Aurora)
//-------------------------------------------------------------------------------------------------------


#ifndef HPCSCAN_DATA_ACQUISITION_NEC_H_
#define HPCSCAN_DATA_ACQUISITION_NEC_H_

#include "data_Acquisition.h"
#include "type_def.h"

using namespace std;

namespace hpcscan
{

    class DataAcquisition_NEC : public DataAcquisition
    {
    public:
        // constructor
        DataAcquisition_NEC(void);


    private:
   
    };

} // namespace hpcscan
#endif