
//-------------------------------------------------------------------------------------------------------
//
// This grid is activated with command line option -testMode Baseline
//
// ## THIS IS THE MAIN CLASS OF HPCSCAN AND THE REFERENCE IMPLEMENTATION ##
// ## ALL GRID IMPLEMENTATIONS SHOULD DERIVE FROM THIS CLASS             ##
//
// Handle all grid data in hpcscan and all operations performed on grids
//
// Details of the grid structure and indexes is provided in grid.cpp
//
//-------------------------------------------------------------------------------------------------------

// ####################################################################################
// ##                                                                                ##
// ##              IMPORTANT: THIS FILE SHOULD NOT BE MODIFIED                       ##
// ## To implement specialization of some functions, you may create a new grid class ##
// ## that derives from this one. See for example: grid_CacheBlk.cpp                 ##
// ##                                                                                ##
// ####################################################################################

#ifndef HPCSCAN_GRID_H_
#define HPCSCAN_GRID_H_

#include <string>

#include "mpi.h"

#include "type_def.h"

namespace hpcscan {

class Grid
{
public:

	// constructor
	Grid(Grid_type) ;

	// constructor
	Grid(Grid_type, Dim_type, Myint64, Myint64, Myint64) ;

	// destructor
	~Grid(void) ;

	// initialise grid index and MPI data structure
	virtual Rtn_code initializeGrid(void) ;

	// print info
	virtual void info(void) ;

	// write grid on disk
	virtual Rtn_code write(string) ;

	// write appropriate part of the local grid on disk
	// each proc writes its own part in a separated file
	virtual Rtn_code write(Point_type pointType, string) ;

	// write appropriate part of the global grid on disk 
	// each proc writes its own part in a global file
	// NOTE: only point type INNER_POINTS is supported in this version
	virtual Rtn_code writeGlobal(Point_type pointType, string file_name, Myint n_iter) ;	

	// read grid values from disk
	virtual Rtn_code read(string) ;

	// read grid values from disk and fill appropriate part of the grid
	virtual Rtn_code read(Point_type pointType, string) ;

	// L1 error between this grid and another
	virtual Myfloat L1Err(Point_type pointType, const Grid& gridIn) const ;

	// collective (all MPI process) L1 error between this grid and another
	virtual Myfloat allProcL1Err(Point_type, const Grid&) const ;

	// Max error between this grid and another (point wise)
	virtual Myfloat maxErr(Point_type, const Grid&) const ;

	// get grid indexes
	Rtn_code getGridIndex(Point_type, Myint64* i1Start, Myint64* i1End,
			Myint64* i2Start, Myint64* i2End, Myint64* i3Start, Myint64* i3End) const ;

	// get min and max in grid
	virtual Myfloat getMin(Point_type) ;
	virtual Myfloat getMax(Point_type) ;

	// collective (all MPI process) get min and max in grid
	virtual Myfloat allProcGetMin(Point_type) ;
	virtual Myfloat allProcGetMax(Point_type) ;

	// get sum of abs grid points
	virtual Myfloat getSumAbs(Point_type) const ;

	// get sum of abs diff of grid points
	virtual Myfloat getSumAbsDiff(Point_type, const Grid&) const;

	// get number of neighbours
	Myint getNumberOfNeighbour() ;

	// get number of grid points according to point type
	Myint64 getNumberOfGridPoint(Grid_type gridType, Point_type pointType) ;

	// get number of grid points involved in exchange halos
	Myint64 getNumberOfGridPointCommHalo(Grid_type gridType) ;

	// get number of grid points involved in boundary condition
	Myint64 getNumberOfGridPointBoundaryCondition(BoundCond_type boundCondType) ;

	// get min coordinates
	Myfloat getMinCoord(Axis_type) ;

	// get max coordinates
	Myfloat getMaxCoord(Axis_type) ;

	// get min spatial sampling
	Myfloat getMinSpaceSampling(void) ;

	// get neighbour proc
	Myint getNeighbourProc(Point_type) ;

	// check grids have same size
	bool sameSize(const Grid&) const ;

	// fill grid with constant value
	virtual void fill(Point_type pointType, Myfloat val) ;

	// fill grid with predefined functions
	virtual void fill(Point_type pType, Func_type t1,  Func_type t2, Func_type t3,
			Myfloat64 param1, Myfloat64 param2, Myfloat64 param3, Myfloat64 amp) ;

	// fill array with constant value
	virtual void fillArray(Myfloat val) ;

	// copy array
	virtual void copyArray(const Grid& gridIn) ;

	// add array
	virtual void addArray(const Grid& gridIn1, const Grid& gridIn2) ;

	// multiply array
	virtual void multiplyArray(const Grid& gridIn1, const Grid& gridIn2) ;

	// add and update array
	virtual void addUpdateArray(const Grid& gridIn) ;

	// define unit grid
	void defineUnitGrid() ;

	// exchange halos with MPI
	virtual Rtn_code exchangeHalos(MPI_comm_mode_type) ;

	// exchange one halo with MPI
	virtual Rtn_code exchangeHalo(MPI_comm_mode_type, Point_type pointType) ;

	// exchange all halos with MPI
	virtual Rtn_code myExchangeAll_halos();

	// apply boundary condition
	virtual Rtn_code applyBoundaryCondition(BoundCond_type boundCondType) ;

	// compute FD_D2 along N1
	virtual Rtn_code FD_D2_N1(Point_type pType, const Grid& Wgrid, Myint fdOrder) ;

	// compute FD_D2 along N2
	virtual Rtn_code FD_D2_N2(Point_type pType, const Grid& Wgrid, Myint fdOrder) ;

	// compute FD_D2 along N3
	virtual Rtn_code FD_D2_N3(Point_type pType, const Grid& Wgrid, Myint fdOrder) ;

	// compute FD_LAPLACIAN
	virtual Rtn_code FD_LAPLACIAN(Point_type pType, const Grid& Wgrid, Myint fdOrder) ;

	// update pressure
	virtual Rtn_code updatePressure(Point_type pType, const Grid& prcGrid,
			const Grid& coefGrid, const Grid& laplaGrid) ;

	// compute pressure with FD
	virtual Rtn_code computePressureWithFD(Grid& prcGridIn, Grid& coefGridIn, Myint fdOrder) ;

	// send grid with MPI_Send
	virtual Rtn_code sendWithMPI(Myint64 nGridPoint, Myint procDestId) ;

	// receive grid with MPI_Recv
	virtual Rtn_code recvWithMPI(Myint64 nGridPoint, Myint procSrcId) ;

	// send and receive grid with MPI_Sendrecv
	virtual Rtn_code sendRecvWithMPI(const Grid& gridDest, Myint idSend, Myint idRecv, Myint64 nGridPoint) ;

	// get number of FLOP per point for operator FD_D2
	Myint getFlopPerPtFD_D2(Myint fdOrder) ;

	// get number of FLOP per point for operator LAPLACIAN
	Myint getFlopPerPtFD_LAPLACIAN(Myint fdOrder) ;

	// get number of points in the stencil for operator FD_D2
	Myint getPtPerStencilFD_D2(Myint fdOrder) ;

	// get number of points in the stencil for operator LAPLACIAN
	Myint getPtPerStencilFD_LAPLACIAN(Myint fdOrder) ;

	// get Offset in regard of global grid
	Myint getOffsetGlobal(Axis_type axis) ;

	// check if global points given are part of local grid
	bool isInMyDomain(Myint64 i1,Myint64 i2, Myint64 i3) ;

	// get Global point from local point
	Myint64 localToGlobal(Axis_type axis, Myint64 i1) ;

	// get Local point from global point
	Myint64 globalToLocal(Axis_type axis, Myint64 i1) ;

	// get index from grid globally
	Rtn_code getGridIndexGlobal(Point_type PointType, Myint64* i1Start, Myint64* i1End,
				Myint64* i2Start, Myint64* i2End, Myint64* i3Start, Myint64* i3End);

	// get number of points in global grid
	Myint64 getGlobalnPoints(Axis_type axis);


	// Space dimension
	Dim_type dim ;

	// Grid type
	Grid_type gridType ;

	// 3D grid array
	Myfloat * grid_3d ;   // this will be allocated by Grid::initializeGrid on the CPU (we should remove it soon)

	// Grid size (local)
	Myint n1, n2, n3 ;

	// grid sampling
	Myfloat64 d1, d2, d3 ;

	// halo width (half of FD stencil)
	Myint haloWidth ;

	// grid mode
	string gridMode ;

	// inner points of Global grid
	Myint64 n1Inner, n2Inner, n3Inner ;

	// array size (points)
	Myint64 npoint ;

	//time measure
	double grid_writeGlobal_time;
	double grid_write_time;
	double i_o_pread_time;

protected:

	
	
	// grid indexes (local)
	Myint i1OffsetStart, i1OffsetEnd, i2OffsetStart, i2OffsetEnd, i3OffsetStart, i3OffsetEnd ;
	Myint i1Halo1Start , i1Halo1End , i2Halo1Start , i2Halo1End , i3Halo1Start , i3Halo1End ;
	Myint i1Halo2Start , i1Halo2End , i2Halo2Start , i2Halo2End , i3Halo2Start , i3Halo2End ;
	Myint i1Layer1Start, i1Layer1End, i2Layer1Start, i2Layer1End, i3Layer1Start, i3Layer1End ;
	Myint i1Layer2Start, i1Layer2End, i2Layer2Start, i2Layer2End, i3Layer2Start, i3Layer2End ;
	Myint i1InnerStart , i1InnerEnd , i2InnerStart , i2InnerEnd , i3InnerStart , i3InnerEnd ;
	Myint i1PadStart   , i1PadEnd   , i2PadStart   , i2PadEnd   , i3PadStart   , i3PadEnd ;

	// origin of local grid
	Myfloat Orig1, Orig2, Orig3 ;

	// grid offset using inner points (with respect to global grid)
	Myint i1OffsetGlobInner, i2OffsetGlobInner, i3OffsetGlobInner ;
	
	// grid offset using all points (with respect to global grid)
	Myint i1OffsetGlobAllPoints, i2OffsetGlobAllPoints, i3OffsetGlobAllPoints ;

	// proc id of neighbours (for subdomain decomposiiton)
	Myint i1ProcIdStart, i1ProcIdEnd, i2ProcIdStart, i2ProcIdEnd, i3ProcIdStart, i3ProcIdEnd ;

	// MPI data type used to exchanges halos (implementation 1)
	MPI_Datatype i1HaloDataType, i2HaloDataType, i3HaloDataType;

	// MPI data type used to exchanges halos (implementation 2)
	MPI_Datatype  i1Halo1DataTypeSend , i1Halo2DataTypeSend,
		i1Halo1DataTypeReceive , i1Halo2DataTypeReceive,

		i2Halo1DataTypeSend , i2Halo2DataTypeSend,
		i2Halo1DataTypeReceive , i2Halo2DataTypeReceive,

		i3Halo1DataTypeSend , i3Halo2DataTypeSend,
		i3Halo1DataTypeReceive , i3Halo2DataTypeReceive ;

	// MPI Data type used to define INNER_POINT in local & global grid 
	MPI_Datatype innerLocalGridType, innerGlobalGridType ;

	// create MPI Type to manipulate Inner points in Global Grid	
	MPI_Datatype createMpiTypeInnerGlobalGrid(void) ;

	// create MPI Type to manipulate Inner points in Local Grid	
	MPI_Datatype createMpiTypeInnerLocalGrid(void) ;

	// grid padding along axis 1, 2 and 3
	// add points at the end of each axis
	virtual void padGridn1(void) ;
	virtual void padGridn2(void) ;
	virtual void padGridn3(void) ;

	// grid offset along axis 1, 2 and 3
	// add points at the beginning of each axis
	virtual void offsetGridn1(void) ;
	virtual void offsetGridn2(void) ;
	virtual void offsetGridn3(void) ;

} ;

} // namespace hpcscan

#endif
