#ifndef HPCSCAN_TESTCASE_TEMPLATE_H_
#define HPCSCAN_TESTCASE_TEMPLATE_H_

#include <string>

#include "testCase.h"

using namespace std;

namespace hpcscan {

class TestCase_Template : protected TestCase
{
public:

	TestCase_Template() ;
	Rtn_code run(void) ;
} ;

} // namespace hpcscan

#endif
