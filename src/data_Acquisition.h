// ####################################################################################
// ##                                                                                ##
// ##              IMPORTANT: THIS FILE SHOULD NOT BE MODIFIED                       ##
// ## To implement specialization of some functions, you may create a new grid class ##
// ## that derives from this one. See for example: data_Acquisition_Cuda.cu                 ##
// ##                                                                                ##
// ####################################################################################

#ifndef HPCSCAN_DATA_ACQUISITION_H_
#define HPCSCAN_DATA_ACQUISITION_H_

#include "grid.h"
#include "type_def.h"

using namespace std;

namespace hpcscan
{

    class DataAcquisition
    {
    public:
        // constructor
        DataAcquisition(void);

        // destructor
        ~DataAcquisition(void);

        // initialize
        virtual Rtn_code initialize(DataAcqui_type acquiType, Grid &coefGrid, Myint nt);

        // apply source terms
        Rtn_code appliSourceTerm(Grid &grid, Modeling_type modType, Myint it, Myfloat dt) ;

        // record traces 
        virtual Rtn_code recordTrace(Grid &grid, Modeling_type modType, Myint it, Myfloat dt) ;

        // write traces on disk
        virtual Rtn_code writeTrace(string fileName) ;

        // time measures
        double trace_gather_time;
        double trace_write_time;
        double record_trace_time;
        double total_timeWriteTrace;

    protected:
        // number of sources
        Myint nSrc;

        Myint nt;

        // source index array
        Myint64 *srcIdx;

        // number of receivers globally
        Myint nRec;

        // number of recievers in this process
        Myint nRecLoc;

        // receiver index array
        Myint64 * tab_recIdxLoc = nullptr;

        Myint * tab_idxglob = nullptr;

        // array to store the traces
        Myfloat *trace;

         //number of sample in local trace
        Myint64 nSampleInTraceLoc;

        //number of sample in global trace
        Myint64 total_nSampleInTrace;
    };

} // namespace hpcscan
#endif