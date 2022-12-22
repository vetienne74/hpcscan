#ifndef HPCSCAN_GLOBAL_H_
#define HPCSCAN_GLOBAL_H_

#include <string>

#include "type_def.h"

#define BLOCK_WRITE_SIZE (1024*1024)

namespace hpcscan {

// debug level switch
extern Debug_level debug ;

// global number of MPI process
extern int nMpiProc ;

// global rank of MPI process
extern int myMpiRank ;

// total memory
extern Myint64 max_mem, current_mem ;

} // namespace hpcscan

#endif
