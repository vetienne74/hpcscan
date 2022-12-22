
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
	virtual Rtn_code initializeGrid(void) ;

	// print info
	virtual void info(void) ;

	// get min and max in grid
	virtual Myfloat getMin(Point_type) ;
	virtual Myfloat getMax(Point_type) ;

	// Max error between this grid and another (point wise)
	virtual Myfloat maxErr(Point_type, const Grid&) const ;

	// exchange halos with MPI
	virtual Rtn_code exchangeHalos(MPI_comm_mode_type) ;

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

	virtual Rtn_code read(string) ;
	
	virtual Rtn_code read(Point_type pointType, string) ;

	virtual Rtn_code write(string fileName);

	virtual Rtn_code write(Point_type pointType, string fileName);

protected:

	// Temporary 3D grid array
	Myfloat * tmp_grid_3d ;

	// flag to know packed stencil can be used
	bool flag_packed_stencil ;

	// grid padding
	virtual void padGridn1(void) ;
	virtual void padGridn2(void) ;
	virtual void padGridn3(void) ;

	//read all inside file overthrust
	Rtn_code read_all_at_once_then_loop(Point_type pointType, string fileName);

	//read inside file overthrust line by line
	Rtn_code read_line_by_line(Point_type pointType, string fileName);

	Rtn_code read_plane_by_plane(Point_type pointType, string fileName);


} ;

} // namespace hpcscan

#endif
