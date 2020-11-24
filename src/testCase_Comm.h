#ifndef HPCSCAN_TESTCASE_COMM_H_
#define HPCSCAN_TESTCASE_COMM_H_

#include <string>

#include "testCase.h"

using namespace std;

namespace hpcscan {

class TestCase_Comm : protected TestCase
{
public:

	TestCase_Comm() ;
	Rtn_code run(void) ;
} ;

} // namespace hpcscan

#endif
