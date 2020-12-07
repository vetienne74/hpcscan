
#ifndef HPCSCAN_GRID_CUDA_H_
#define HPCSCAN_GRID_CUDA_H_

#include <string>

#include "mpi.h"

#include "grid.h"
#include "type_def.h"

namespace hpcscan {

class Grid_Cuda : public Grid
{
public:

	// constructor
	Grid_Cuda(Grid_type) ;

	// constructor
	Grid_Cuda(Grid_type, Dim_type, Myint64, Myint64, Myint64) ;

	// destructor
	~Grid_Cuda(void) ;

	// compute pressure with FD
	virtual Rtn_code computePressureWithFD(Grid& prcGridIn, Grid& coefGridIn, Myint fdOrder) ;

	// initialise grid index and MPI data structure
	virtual void initializeGrid(void) ;

	// print info
	virtual void info(void) ;

	// compute FD_LAPLACIAN
	virtual Rtn_code FD_LAPLACIAN(Point_type pType, const Grid&, Myint fdOrder) ;

	// fill grid with constant value
	virtual void fill(Point_type pointType, Myfloat val) ;

	// fill grid with predefined functions
	virtual void fill(Point_type pType, Func_type t1,  Func_type t2, Func_type t3,
			Myfloat64 param1, Myfloat64 param2, Myfloat64 param3, Myfloat64 amp) ;

	// get min and max in grid
	virtual Myfloat getMin(Point_type) ;
	virtual Myfloat getMax(Point_type) ;

	// L1 error between this grid and another
	virtual Myfloat L1Err(Point_type pointType, const Grid& gridIn) const ;

	// update pressure
	virtual Rtn_code updatePressure(Point_type pType, const Grid& prcGrid,
			const Grid& coefGrid, const Grid& laplaGrid) ;

	// apply boundary condition
	virtual Rtn_code applyBoundaryCondition(BoundCond_type boundCondType) ;

	// get sum of abs grid points
	virtual Myfloat getSumAbs(Point_type) const ;

	// get sum of abs diff of grid points
	virtual Myfloat getSumAbsDiff(Point_type, const Grid&) const;

protected:

} ;

} // namespace hpcscan

#endif
