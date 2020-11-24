
#ifndef HPCSCAN_PROPAGATOR_FACTORY_H_
#define HPCSCAN_PROPAGATOR_FACTORY_H_

#include <memory>
#include <string>

#include "propagator_Ac2.h"
#include "type_def.h"

namespace hpcscan {

class Propagator_Factory
{
public:

	// create Propagator
	static shared_ptr<Propagator_Ac2> create(string propaType) ;

} ;

} // namespace hpcscan

#endif
