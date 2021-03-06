
//-------------------------------------------------------------------------------------------------------
// This grid is activated with command line option -testMode NEC
// Derived class from Grid
// NEC compiler directives (target NEC Aurora Vector Engine)
//-------------------------------------------------------------------------------------------------------

#ifndef HPCSCAN_GRID_NEC_H_
#define HPCSCAN_GRID_NEC_H_

#include <string>

#include "mpi.h"

#include "grid.h"
#include "type_def.h"

namespace hpcscan {

class Grid_NEC : public Grid
{
public:

	// constructor
	Grid_NEC(Grid_type) ;

	// constructor
	Grid_NEC(Grid_type, Dim_type, Myint64, Myint64, Myint64) ;

	// destructor
	~Grid_NEC(void) ;

	// initialise grid index and MPI data structure
	virtual void initializeGrid(void) ;

	// print info
	virtual void info(void) ;

	// get min and max in grid
	virtual Myfloat getMin(Point_type) ;
	virtual Myfloat getMax(Point_type) ;

	// Max error between this grid and another (point wise)
	virtual Myfloat maxErr(Point_type, const Grid&) const ;

	// apply boundary condition
	virtual Rtn_code applyBoundaryCondition(BoundCond_type boundCondType) ;

	// compute FD_D2 along N1
	virtual Rtn_code FD_D2_N1(Point_type pType, const Grid& Wgrid, Myint fdOrder) ;

	// compute FD_D2 along N2
	virtual Rtn_code FD_D2_N2(Point_type pType, const Grid& Wgrid, Myint fdOrder) ;

	// compute FD_D2 along N3
	virtual Rtn_code FD_D2_N3(Point_type pType, const Grid& Wgrid, Myint fdOrder) ;

	// compute FD_LAPLACIAN
	virtual Rtn_code FD_LAPLACIAN(Point_type pType, const Grid& Wgrid, Myint fdOrder) ;

	// update pressure
	virtual Rtn_code updatePressure(Point_type pType, const Grid& prcGrid,
			const Grid& coefGrid, const Grid& laplaGrid) ;

	// compute pressure with FD
	virtual Rtn_code computePressureWithFD(Grid& prcGridIn, Grid& coefGridIn, Myint fdOrder) ;

protected:

#ifndef _DOUBLE_PRECISION_
	// flag to know packed stencil can be used
	bool flag_packed_stencil ;
#endif

	// grid padding
	virtual void padGridn1(void) ;
	virtual void padGridn2(void) ;
	virtual void padGridn3(void) ;

} ;

} // namespace hpcscan

#endif
