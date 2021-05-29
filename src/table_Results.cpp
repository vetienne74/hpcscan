
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

void Table_Results::seOneValue(Myint iLine, Myint iCol, Myfloat value)
{
	printDebug(FULL_DEBUG, "IN Table_Results::setOneValue");

	if ((iLine) < 0 || (iLine > nLine))
	{
		printError("Table_Results::seOneValue, invalid iLine", iLine) ;
		return ;
	}
	if ((iCol) < 0 || (iCol > nCol))
	{
		printError("Table_Results::seOneValue, invalid iCol", iCol) ;
		return ;
	}
	Myint idx = iLine + iCol * nLine ;
	val[idx] = value ;

	printDebug(FULL_DEBUG, "OUT Table_Results::setOneValue");
}

//-------------------------------------------------------------------------------------------------------

void Table_Results::display(void)
{
	printDebug(MID_DEBUG, "IN Table_Results::display");

	if (myMpiRank == 0)
	{
		cout << "\n * Table " << tableName << " *\n" ;

		// display column labels
		cout << "         | " ;
		for (Myint iCol = 0; iCol < nCol ; iCol++)
		{
			printf("%6d | ", iCol) ;
		}
		cout << "\n" ;

		// display line delimiter
		cout << "|--------|" ;
		for (Myint iCol = 0; iCol < nCol ; iCol++)
		{
			cout << "--------|" ;
		}
		cout << "\n" ;

		// display values
		for (Myint iLine = 0; iLine < nLine; iLine++)
		{
			printf("| %6d | ", iLine) ;
			for (Myint iCol = 0; iCol < nCol ; iCol++)
			{
				Myint idx = iLine + iCol * nLine ;
				if (val[idx] != UNSPECIFIED)
				{
					printf("%6.2f | ", val[idx]) ;
				}
				else
				{
					cout << "   -   | " ;
				}
			}
			cout << "\n" ;
		}
	}

	printDebug(MID_DEBUG, "OUT Table_Results::display");
}

} // namespace hpcscan



