#ifndef HPCSCAN_TESTCASE_FD_D2_H_
#define HPCSCAN_TESTCASE_FD_D2_H_

#include <string>

#include "testCase.h"

using namespace std;

namespace hpcscan {

class TestCase_FD_D2 : protected TestCase
{
public:

	TestCase_FD_D2() ;
	virtual Rtn_code run(void) ;
} ;

} // namespace hpcscan

#endif
