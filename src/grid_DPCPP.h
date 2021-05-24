
//-------------------------------------------------------------------------------------------------------
// This grid is activated with command line option -testMode DPC++
// Derived class from Grid
// DPC++ implementation (target CPU, GPU & FPGA)
//-------------------------------------------------------------------------------------------------------

#ifndef HPCSCAN_Grid_DPCPP_H_
#define HPCSCAN_Grid_DPCPP_H_

#include <string>

#include <CL/sycl.hpp>
#include "mpi.h"

#include "grid.h"
#include "type_def.h"

namespace hpcscan {

class Grid_DPCPP : public Grid
{
public:

	// constructor
	Grid_DPCPP(Grid_type) ;

	// constructor
	Grid_DPCPP(Grid_type, Dim_type, Myint64, Myint64, Myint64) ;

	// destructor
	~Grid_DPCPP(void) ;

	// print info
	virtual void info(void) ;

	// initialise grid index and MPI data structure
	virtual void initializeGrid(void) ;

	// fill array with constant value
	virtual void fillArray(Myfloat val) ;

protected:

	// pointer to device memory
	// equivalent to grid_3d on the host side
	Myfloat * d_grid_3d ;

	// pointer to device memory
	Myfloat * d_help_3d ;
	Myfloat * d_help_3d_2 ;

	// pointer to device queue
	// TO DO it should be initialized in initializeGrid
	sycl::queue Q ;

} ;

} // namespace hpcscan

#endif
