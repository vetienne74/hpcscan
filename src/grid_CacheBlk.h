
//-------------------------------------------------------------------------------------------------------
// This grid is activated with command line option -testMode CacheBlk
// Derived class from Grid
// Optimized with cache blocking techniques (target CPU)
//-------------------------------------------------------------------------------------------------------

#ifndef HPCSCAN_GRID_CACHEBLK_H_
#define HPCSCAN_GRID_CACHEBLK_H_

#include <string>

#include "mpi.h"

#include "grid.h"
#include "type_def.h"

namespace hpcscan {

class Grid_CacheBlk : public Grid
{
public:

	// constructor
	Grid_CacheBlk(Grid_type) ;

	// constructor
	Grid_CacheBlk(Grid_type, Dim_type, Myint64, Myint64, Myint64) ;

	// destructor
	~Grid_CacheBlk(void) ;

	// print info
	virtual void info(void) ;

	// compute FD_D2 along N1
	virtual Rtn_code FD_D2_N1(Point_type pType, const Grid&, Myint fdOrder) ;

	// compute FD_D2 along N2
	virtual Rtn_code FD_D2_N2(Point_type pType, const Grid&, Myint fdOrder) ;

	// compute FD_D2 along N3
	virtual Rtn_code FD_D2_N3(Point_type pType, const Grid&, Myint fdOrder) ;

	// compute FD_LAPLACIAN
	virtual Rtn_code FD_LAPLACIAN(Point_type pType, const Grid&, Myint fdOrder) ;

	// compute pressure with FD
	virtual Rtn_code computePressureWithFD(Grid& prcGridIn, Grid& coefGridIn, Myint fdOrder) ;

protected:

} ;

} // namespace hpcscan

#endif
