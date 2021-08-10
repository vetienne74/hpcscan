
//-------------------------------------------------------------------------------------------------------
// Grid factory used to create specialized Grid objects
// All grid implementations are referenced
//-------------------------------------------------------------------------------------------------------

#ifndef HPCSCAN_GRID_FACTORY_H_
#define HPCSCAN_GRID_FACTORY_H_

#include <memory>
#include <string>

#include "grid.h"
#include "type_def.h"

namespace hpcscan {

class Grid_Factory
{
public:

	// create Grid
	static shared_ptr<Grid> create(string gridMode, Grid_type gridType) ;
	static shared_ptr<Grid> create(string gridMode, Grid_type gridType, Dim_type dim,
			Myint64 n1, Myint64 n2, Myint64 n3) ;


} ;

} // namespace hpcscan

#endif
