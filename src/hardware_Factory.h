
//-------------------------------------------------------------------------------------------------------
// Hardware factory used to create specialized Hardware objects
// All implementations are referenced
//-------------------------------------------------------------------------------------------------------

#ifndef HPCSCAN_HARDWARE_FACTORY_H_
#define HPCSCAN_HARDWARE_FACTORY_H_

#include <memory>
#include <string>

#include "hardware.h"
#include "type_def.h"

namespace hpcscan {

class Hardware_Factory
{
public:

	// create Hardware
	static shared_ptr<Hardware> create(string gridMode) ;
} ;

} // namespace hpcscan

#endif
