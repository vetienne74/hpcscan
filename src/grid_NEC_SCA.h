
#ifndef HPCSCAN_GRID_NEC_SCA_H_
#define HPCSCAN_GRID_NEC_SCA_H_

#include <string>

#include "mpi.h"

#ifdef __NEC__
#include <sca.h>
#endif

#include "grid_NEC.h"
#include "type_def.h"

namespace hpcscan {

class Grid_NEC_SCA : public Grid_NEC
{
public:

	// constructor
	Grid_NEC_SCA(Grid_type) ;

	// constructor
	Grid_NEC_SCA(Grid_type, Dim_type, Myint64, Myint64, Myint64) ;

	// destructor
	~Grid_NEC_SCA(void) ;

	// initialise grid index and MPI data structure
	virtual void initializeGrid(void) ;

	// print info
	virtual void info(void) ;

	// compute FD_D2 along N1
	virtual Rtn_code FD_D2_N1(Point_type pType, const Grid& Wgrid, Myint fdOrder) ;

	// compute FD_D2 along N2
	virtual Rtn_code FD_D2_N2(Point_type pType, const Grid& Wgrid, Myint fdOrder) ;

	// compute FD_D2 along N3
	virtual Rtn_code FD_D2_N3(Point_type pType, const Grid& Wgrid, Myint fdOrder) ;

	// compute FD_LAPLACIAN
	virtual Rtn_code FD_LAPLACIAN(Point_type pType, const Grid& Wgrid, Myint fdOrder) ;

protected:

#ifdef __NEC__
	// Handle for stencil code.
	sca_code_t code_FD_D2_N1 ;
	sca_code_t code_FD_D2_N2 ;
	sca_code_t code_FD_D2_N3 ;
	sca_code_t code_FD_LAPLACIAN ;
#endif

	// flag for SCA initialization
	bool flag_code_FD_D2_N1 ;
	bool flag_code_FD_D2_N2 ;
	bool flag_code_FD_D2_N3 ;
	bool flag_code_FD_LAPLACIAN ;

	// initialize SCA
	Rtn_code initialize_code_FD_D2_N1(Point_type pType, const Grid& Ugrid, Myint fdOrder) ;
	Rtn_code initialize_code_FD_D2_N2(Point_type pType, const Grid& Ugrid, Myint fdOrder) ;
	Rtn_code initialize_code_FD_D2_N3(Point_type pType, const Grid& Ugrid, Myint fdOrder) ;
	Rtn_code initialize_code_FD_LAPLACIAN(Point_type pType, const Grid& Ugrid, Myint fdOrder) ;

} ;

} // namespace hpcscan

#endif
