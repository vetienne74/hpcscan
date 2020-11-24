#ifndef HPCSCAN_TESTCASE_MEMORY_H_
#define HPCSCAN_TESTCASE_MEMORY_H_

#include <string>

#include "testCase.h"

using namespace std;

namespace hpcscan {

class TestCase_Memory : protected TestCase
{
public:

	TestCase_Memory() ;
	Rtn_code run(void) ;
} ;

} // namespace hpcscan

#endif
