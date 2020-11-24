#ifndef HPCSCAN_TESTCASE_UTIL_H_
#define HPCSCAN_TESTCASE_UTIL_H_

#include <string>

#include "testCase.h"

using namespace std;

namespace hpcscan {

class TestCase_Util : protected TestCase
{
public:

	TestCase_Util() ;
	Rtn_code run(void) ;
} ;

} // namespace hpcscan

#endif
