#ifndef HPCSCAN_TESTCASE_H_
#define HPCSCAN_TESTCASE_H_

#include <fstream>
#include <string>

#include "grid.h"
#include "hardware_Factory.h"
#include "type_def.h"

using namespace std;

namespace hpcscan {

class TestCase
{
public:
	virtual Rtn_code run(void) = 0 ;

protected:

	// initialize test case
	virtual Rtn_code initialize(void) ;

	// finalize test case
	virtual void finalize(void) ;

	// check L1 error between 2 grids
	Rtn_code checkGridL1Err(Point_type pointType, const Grid&, const Grid&, Myfloat) ;

	// collective (all MPI proc) check L1 error between 2 grids
	Rtn_code checkAllProcGridL1Err(Point_type pointType, const Grid&, const Grid&, Myfloat) ;

	// check Max. error between 2 grids
	Rtn_code checkGridMaxErr(Point_type pointType, const Grid&, const Grid&, Myfloat) ;

	// collective (all MPI proc) check Max. error between 2 grids
	Rtn_code checkAllProcGridMaxErr(Point_type pointType, const Grid&, const Grid&, Myfloat) ;

	// check error between 2 floats
	Rtn_code checkFloatDiff(Myfloat, Myfloat, Myfloat) ;

	// check 2 integers are equal
	Rtn_code checkIntegerDiff(Myint, Myint) ;

	// check 2 boolean are equal
	Rtn_code checkBoolDiff(bool, bool) ;

	// relative error between 2 float
	Myfloat relErr(Myfloat float1, Myfloat float2) ;

	// performance log file
	// elapse time, GFlops, GB/s, etc
	ofstream perfLogFile ;

	// timer
	double testCaseStart, testCaseEnd ;

	// test case name
	string testCaseName ;

	// test case version
	string testCaseVersion ;

	// hardware
	shared_ptr<Hardware> hw ;
} ;

} // namespace hpcscan

#endif
