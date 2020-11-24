#ifndef HPCSCAN_TESTCASE_PROPA_H_
#define HPCSCAN_TESTCASE_PROPA_H_

#include <string>

#include "testCase.h"

using namespace std;

namespace hpcscan {

class TestCase_Propa : protected TestCase
{
public:

	TestCase_Propa() ;
	virtual Rtn_code run(void) ;
} ;

} // namespace hpcscan

#endif
