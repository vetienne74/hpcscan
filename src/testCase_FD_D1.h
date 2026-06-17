#ifndef HPCSCAN_TESTCASE_FD_D1_H_
#define HPCSCAN_TESTCASE_FD_D1_H_

#include <string>

#include "testCase.h"

using namespace std;

namespace hpcscan {

class TestCase_FD_D1 : protected TestCase
{
public:

	TestCase_FD_D1() ;
	virtual Rtn_code run(void) ;
} ;

} // namespace hpcscan

#endif
