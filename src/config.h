
// this class is a singleton pattern

#ifndef HPCSCAN_CONFIG_H_
#define HPCSCAN_CONFIG_H_

#include <climits>
#include <fstream>
#include <string>

#include "type_def.h"

namespace hpcscan {

class Config
{
public:

	static Config* Instance() ;

	Rtn_code parse_argument(int argc, char* argv[]) ;

	Rtn_code info(void) ;

	Rtn_code initialize(void) ;

	Rtn_code finalize(Rtn_code) ;

	// auto padding (final grid is automatically augmented)
	bool autoPad ;

	// boundary name
	string boundary ;

	// cache block sizes
	Myint cb1, cb2, cb3 ;

	// debug log file
	ofstream debugLogFile ;

	// space dimension
	Dim_type dim ;

	// time step (for propagator)
	Myfloat64 dt ;

	// fd order
	Myint fdOrder ;

	// grid sampling
	Myfloat64 h ;

	// GPU block size (number of threads per block)
	Myint gpuBlkSize ;

	// GPU MPI aware flag
	bool gpuMpiAware ;

	// GPU grid size (number of blocks per grid)
	Myint gpuGridSize ;

	// inner domain size
	Myint n1, n2, n3 ;

	// grid padding (final grid is augmented by given number)
	Myint n1AddPad, n2AddPad, n3AddPad ;

	// grid padding (final grid is multiple of given number)
	Myint n1MulPad, n2MulPad, n3MulPad ;

	// layer width
	Myint nlayer ;

	// no. of time step (for propagator)
	Myint nt ;

	// no. of tries
	Myint ntry ;

	// subdomain decomposiiton
	Myint nsub1, nsub2, nsub3 ;

	// coef. parameters used for various purposes
	Myfloat64 param1, param2, param3, param4 ;

	// propagator name
	string propagator ;

	// ratio of stability time step
	Myfloat64 ratioCFL ;

	// time step increment (time in sec.) between snaphots for propagator
	Myfloat64 snapDt ;

	// time step increment (time step) between snaphots for propagator
	Myint snapInc ;

	// name of testCase to run ("All" to run all testCases)
	string testCaseName ;

	// name of testMode
	string testMode ;

	// maximum time (for propagator)
	Myfloat64 tmax ;

	// write grids
	bool writeGrid ;

	// host name
	char hostName[HOST_NAME_MAX];

	// user name
	char userName[LOGIN_NAME_MAX];

private:

	// private so that it can  not be called
	Config(void) ;

	// copy constructor is private
	Config(Config const&){} ;

	// assignment operator is private
	Config& operator=(Config const&){return *this;} ;

	// singleton
	static Config* m_pInstance ;

} ;

} // namespace hpcscan

#endif
