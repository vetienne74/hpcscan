#ifndef HPCSCAN_TESTCASE_MATRIX_H_
#define HPCSCAN_TESTCASE_MATRIX_H_

#include <string>

#include "testCase.h"

using namespace std;

namespace hpcscan
{

	class TestCase_Matrix : protected TestCase
	{
	public:
		TestCase_Matrix();
		Rtn_code run(void);
	};

} // namespace hpcscan

#endif
