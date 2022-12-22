
#ifndef HPCSCAN_DATA_ACQUISITION_FACTORY_H_
#define HPCSCAN_DATA_ACQUISITION_FACTORY_H_

#include <memory>
#include <string>

#include "data_Acquisition.h"
#include "type_def.h"

namespace hpcscan {

class DataAcquisition_Factory
{
public:

	// create data_Acquisition
	static shared_ptr<DataAcquisition> create(string testMode) ;

} ;

} // namespace hpcscan

#endif
