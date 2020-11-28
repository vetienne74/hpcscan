
//-------------------------------------------------------------------------------------------------------
// Grid factory used to create specialized Grid objects
//-------------------------------------------------------------------------------------------------------

#include "grid_Factory.h"

#include "grid_CacheBlk.h"
#include "grid_GPU1.h"
#include "grid_NEC_SCA.h"
#include "output_report.h"

using namespace std;

namespace hpcscan {

//-------------------------------------------------------------------------------------------------------

shared_ptr<Grid> Grid_Factory::create(string gridMode, Grid_type gridType)
{
	printDebug(MID_DEBUG, "IN Grid_Factory::create");

	Grid *Rgrid = nullptr ;

	if (gridMode.compare("Baseline") == 0)
	{
		Rgrid = new Grid(gridType) ;
	}
	else if (gridMode.compare("CacheBlk") == 0)
	{
		Rgrid = new Grid_CacheBlk(gridType) ;
	}
	else if (gridMode.compare("GPU1") == 0)
	{
		Rgrid = new Grid_GPU1(gridType) ;
	}
	else if (gridMode.compare("NEC_SCA") == 0)
	{
		Rgrid = new Grid_NEC_SCA(gridType) ;
	}
	else
	{
		printError("IN Grid_Factory::create, invalid gridMode", gridMode) ;
	}

	printDebug(MID_DEBUG, "OUT Grid_Factory::create");
	if (Rgrid != nullptr)
	{
		return shared_ptr<Grid>(Rgrid) ;
	}
	else
	{
		printError("IN Grid_Factory::create, return nullptr", gridMode) ;
		return nullptr ;
	}
}

//-------------------------------------------------------------------------------------------------------

shared_ptr<Grid> Grid_Factory::create(string gridMode, Grid_type gridType, Dim_type dim,
		Myint64 n1, Myint64 n2, Myint64 n3)
{
	printDebug(MID_DEBUG, "IN Grid_Factory::create");

	Grid *Rgrid = nullptr ;

	if (gridMode.compare("Baseline") == 0)
	{
		Rgrid = new Grid(gridType, dim, n1, n2, n3) ;
	}
	else if (gridMode.compare("CacheBlk") == 0)
	{
		Rgrid = new Grid_CacheBlk(gridType, dim, n1, n2, n3) ;
	}
	else if (gridMode.compare("GPU1") == 0)
	{
		Rgrid = new Grid_GPU1(gridType, dim, n1, n2, n3) ;
	}
	else if (gridMode.compare("NEC_SCA") == 0)
	{
		Rgrid = new Grid_NEC_SCA(gridType, dim, n1, n2, n3) ;
	}
	else
	{
		printError("IN Grid_Factory::create, invalid gridMode") ;
	}

	printDebug(MID_DEBUG, "OUT Grid_Factory::create");
	if (Rgrid != nullptr)
	{
		return shared_ptr<Grid>(Rgrid) ;
	}
	else
	{
		return nullptr ;
	}
}


} // namespace hpcscan
