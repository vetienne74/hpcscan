
//-------------------------------------------------------------------------------------------------------
// This grid is activated with command line option -testMode CUDA_Optim
// Derived class from Grid_Cuda
// CUDA implementation (target GPU)
//-------------------------------------------------------------------------------------------------------

#ifndef HPCSCAN_GRID_CUDA_OPTIM_H_
#define HPCSCAN_GRID_CUDA_OPTIM_H_

#include "grid_Cuda.h"

namespace hpcscan {

class Grid_Cuda_Optim : public Grid_Cuda
{
public:

	// constructor
	Grid_Cuda_Optim(Grid_type) ;

	// constructor
	Grid_Cuda_Optim(Grid_type, Dim_type, Myint64, Myint64, Myint64) ;

	// destructor
	~Grid_Cuda_Optim(void) ;

	// initialise grid index and MPI data structure
	virtual Rtn_code initializeGrid(void) ;

	// print info
	virtual void info(void) ;

	// compute FD_LAPLACIAN
	virtual Rtn_code FD_LAPLACIAN(Point_type pType, const Grid&, Myint fdOrder) ;

	// compute pressure with FD
	virtual Rtn_code computePressureWithFD(Grid& prcGridIn, Grid& coefGridIn, Myint fdOrder) ;

protected:

} ;

} // namespace hpcscan

#endif
