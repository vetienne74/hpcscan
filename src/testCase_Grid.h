#ifndef HPCSCAN_TESTCASE_GRID_H_
#define HPCSCAN_TESTCASE_GRID_H_

#include <string>

#include "testCase.h"

using namespace std;

namespace hpcscan {

class TestCase_Grid : protected TestCase
{
public:

	TestCase_Grid() ;
	Rtn_code run(void) ;
} ;

} // namespace hpcscan

#endif
