#ifndef TABLE_RESULTS_H_
#define TABLE_RESULTS_H_

#include <string>

#include "type_def.h"

using namespace std;

namespace hpcscan {

class Table_Results
{
public:

	// constructor
	Table_Results(string tableName, Myint nLine, Myint nCol) ;

	// destructor
	~Table_Results(void) ;

	// output table in terminal
	void display(void) ;

	// set on values in the table
	void seOneValue(Myint iLine, Myint iCol, Myfloat value) ;

	// table name
	string tableName ;

	// number of lines and columns
	Myint nLine, nCol ;

	// total number of values
	Myint nVal ;

	// table of values
	Myfloat * val ;

};

} // namespace hpcscan

#endif
