#ifndef HPCSCAN_TESTCASE_MODELING_H_
#define HPCSCAN_TESTCASE_MODELING_H_

#include <string>

#include "testCase.h"

using namespace std;

namespace hpcscan {

class TestCase_Modeling : protected TestCase
{
public:

	TestCase_Modeling() ;
	Rtn_code run(void) ;
} ;

} // namespace hpcscan

#endif
