#ifndef HPCSCAN_CONSTANT_H_
#define HPCSCAN_CONSTANT_H_

#include "type_def.h"

#include "mpi.h"

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

} // namespace hpcscan

#endif
