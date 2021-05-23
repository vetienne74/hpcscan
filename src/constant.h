#ifndef HPCSCAN_CONSTANT_H_
#define HPCSCAN_CONSTANT_H_

#include "type_def.h"

namespace hpcscan {

//-------------------------------------------------------------------------------------------------------
// Useful numerical values
//-------------------------------------------------------------------------------------------------------

const Myfloat64 PI = 3.1415926535897 ;

#ifdef _DOUBLE_PRECISION_
const Myfloat MAX_ERR_FLOAT = 1.e-11 ;
#else
const Myfloat MAX_ERR_FLOAT = 1.e-5 ;
#endif

const Myfloat TWO = 2.0 ;

//-------------------------------------------------------------------------------------------------------
// Grid modes
//-------------------------------------------------------------------------------------------------------
const string   GRID_MODE_BASELINE      = "Baseline" ;
const string   GRID_MODE_CACHEBLK      = "CacheBlk" ;
const string   GRID_MODE_CUDA          = "CUDA" ;
const string   GRID_MODE_DPCPP         = "DPC++" ;
const string   GRID_MODE_HIP           = "HIP" ;
const string   GRID_MODE_OPENACC       = "OpenAcc" ;
const string   GRID_MODE_NEC           = "NEC" ;
const string   GRID_MODE_NEC_SCA       = "NEC_SCA" ;

//-------------------------------------------------------------------------------------------------------
// Propagator types
//-------------------------------------------------------------------------------------------------------
const string   PROPA_TYPE_AC2STANDARD  = "Ac2Standard" ;
const string   PROPA_TYPE_AC2SPLITCOMP = "Ac2SplitComp" ;

} // namespace hpcscan

#endif
