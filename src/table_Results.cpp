
//-------------------------------------------------------------------------------------------------------
// Management of results in a table format
//-------------------------------------------------------------------------------------------------------

#include "table_Results.h"

#include <iostream>

#include "global.h"
#include "output_report.h"

namespace hpcscan {

//-------------------------------------------------------------------------------------------------------

Table_Results::Table_Results(string tableNameIn, Myint nLineIn, Myint nColIn)
		{
	printDebug(MID_DEBUG, "IN Table_Results::Table_Results");

	tableName = tableNameIn ;
	nLine     = nLineIn ;
	nCol      = nColIn ;
	nVal      = nLine * nCol ;

	// allocate table
	val = new Myfloat[nLine*nCol] ;

	// initialize table
	for (Myint ii=0; ii<nVal; ii++)
	{
		val[ii] = UNSPECIFIED ;
	}

	printDebug(MID_DEBUG, "OUT Table_Results::Table_Results");
		}

//-------------------------------------------------------------------------------------------------------

Table_Results::~Table_Results(void)
		{
	printDebug(MID_DEBUG, "IN Table_Results::~Table_Results");

	delete[] val ;

	printDebug(MID_DEBUG, "OUT Table_Results::~Table_Results");
		}

//-------------------------------------------------------------------------------------------------------

void Table_Results::display(void)
{
	printDebug(MID_DEBUG, "IN Table_Results::display");

	if (myMpiRank == 0)
	{
		cout << "\n * Table " << tableName << " *\n" ;

		for (Myint iCol=0; iCol<nCol; iCol++)
		{
			cout << " | " ;
			for (Myint iLine=0; iLine<nLine; iLine++)
			{
				Myint idx = iLine + iCol * nLine ;
				cout << val[idx] << " | " ;
			}
			cout << "\n" ;
		}

	}

	printDebug(MID_DEBUG, "OUT Table_Results::display");
}

} // namespace hpcscan



