
//-------------------------------------------------------------------------------------------------------
// This grid is activated with command line option -testMode HIP_Optim
// Derived class from Grid_Hip
// HIP implementation (target GPU)
//-------------------------------------------------------------------------------------------------------

#ifndef HPCSCAN_GRID_HIP_OPTIM_H_
#define HPCSCAN_GRID_HIP_OPTIM_H_

#include "grid_Hip.h"

namespace hpcscan {

class Grid_Hip_Optim : public Grid_Hip
{
public:

	// constructor
	Grid_Hip_Optim(Grid_type) ;

	// constructor
	Grid_Hip_Optim(Grid_type, Dim_type, Myint64, Myint64, Myint64) ;

	// destructor
	~Grid_Hip_Optim(void) ;

	// initialise grid index and MPI data structure
	virtual Rtn_code initializeGrid(void) ;

	// print info
	virtual void info(void) ;

	// compute FD_D2 along N1
	virtual Rtn_code FD_D2_N1(Point_type pType, const Grid& Wgrid, Myint fdOrder) ;

	// compute FD_D2 along N2
	virtual Rtn_code FD_D2_N2(Point_type pType, const Grid& Wgrid, Myint fdOrder) ;

	// compute FD_D2 along N3
	virtual Rtn_code FD_D2_N3(Point_type pType, const Grid& Wgrid, Myint fdOrder) ;

	// compute FD_LAPLACIAN
	virtual Rtn_code FD_LAPLACIAN(Point_type pType, const Grid&, Myint fdOrder) ;

protected:

	int gpuBlkSize1 ;
	int gpuBlkSize2 ;
	int gpuBlkSize3 ;

} ;

} // namespace hpcscan

#endif
