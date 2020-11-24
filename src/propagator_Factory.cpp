
//-------------------------------------------------------------------------------------------------------
// Propagator factory used to create specialized Propagator objects
//-------------------------------------------------------------------------------------------------------

#include "propagator_Factory.h"

#include "propagator_Ac2.h"
#include "propagator_Ac2SplitComp.h"
#include "output_report.h"

using namespace std;

namespace hpcscan {

//-------------------------------------------------------------------------------------------------------

shared_ptr<Propagator_Ac2> Propagator_Factory::create(string propaName)
{
	printDebug(MID_DEBUG, "IN Propagator_Factory::create");

	Propagator_Ac2 *propa = nullptr ;

	if (propaName.compare("Ac2Standard") == 0)
	{
		propa = new Propagator_Ac2() ;
	}
	else if (propaName.compare("Ac2SplitComp") == 0)
	{
		propa = new Propagator_Ac2SplitComp() ;
	}
	else
	{
		printError("IN Propagator_Factory::create, invalid propaName", propaName) ;
	}

	printDebug(MID_DEBUG, "OUT Propagator_Factory::create");
	if (propa != nullptr)
	{
		return shared_ptr<Propagator_Ac2>(propa) ;
	}
	else
	{
		return nullptr ;
	}
}

} // namespace hpcscan
